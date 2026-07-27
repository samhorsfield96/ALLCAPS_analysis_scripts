library(MLmetrics)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)

# point to data files
data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")

# ── ALLCAPS-typeable serotypes ─────────────────────────────────────────────────
allcaps_serotypes <- read_csv(
  file.path(data_root, "ALLCAPS_possible_serotypes.csv"),
  show_col_types = FALSE
)$Serotypes

# ── Load data ──────────────────────────────────────────────────────────────────
wide <- read_csv(file.path(data_root, "merged_benchmark_results_wide.csv"),
                 show_col_types = FALSE) %>%
  mutate(true_serogroup = as.character(true_serogroup)) %>%
  # remove ambiguous true serotype calls
  filter(true_serotype %in% allcaps_serotypes)

tools <- setdiff(names(wide), c("sample_id", "benchmark", "true_serotype", "true_serogroup"))

# ── Match functions ────────────────────────────────────────────────────────────
exact_match_serotype  <- function(vec, cls) !is.na(vec) & vec == cls
exact_match_serogroup <- function(vec, cls) {
  sub_cls <- sub("^([0-9]+).*", "\\1", cls)
  sub_vec = sub("^([0-9]+).*", "\\1", vec)
  !is.na(sub_vec) & sub_vec == sub_cls
} 
within_match <- function(vec, cls) !is.na(vec) & grepl(cls, vec, fixed = TRUE)

# ── Helper: compute per-class TP/FP/FN, sensitivity, specificity, F1 ──────────
# checked, all good
compute_metrics <- function(true_vec, pred_vec, match_fn) {
  classes <- sort(unique(true_vec[!is.na(true_vec)]))

  map_dfr(classes, function(cls) {
    actual_pos    <- match_fn(true_vec, cls)
    predicted_pos <- match_fn(pred_vec, cls)

    TP <- sum( actual_pos &  predicted_pos, na.rm = TRUE)
    FP <- sum(!actual_pos &  predicted_pos, na.rm = TRUE)
    FN <- sum( actual_pos & !predicted_pos, na.rm = TRUE)
    TN <- sum(!actual_pos & !predicted_pos, na.rm = TRUE)
    Total <- TP + FN
    
    f1 <- F1_Score(y_true = actual_pos, y_pred = predicted_pos, positive = TRUE)
    sensitivity <- Sensitivity(y_true = actual_pos, y_pred = predicted_pos, positive = TRUE)
    specificity <- Specificity(y_true = actual_pos, y_pred = predicted_pos, positive = TRUE)
    accuracy <- Accuracy(y_true = actual_pos, y_pred = predicted_pos)
    precision <- Precision(y_true = actual_pos, y_pred = predicted_pos, positive = TRUE)

    row <- tibble(class = as.character(cls), TP, FP, FN, TN, Total, sensitivity, specificity, precision, accuracy, f1)
    row <- row %>% replace(is.na(.), 0.0)
    row
  })
}

# ── Helper: build confusion matrix (true rows × predicted cols) ───────────────
# Rows = true classes (n), Cols = all unique predicted values (m).
build_confusion <- function(true_vec, pred_vec, true_classes) {
  pred_classes <- sort(unique(pred_vec[!is.na(pred_vec)]))
  all_cols     <- sort(unique(c(true_classes, pred_classes)))

  mat <- matrix(0L,
                nrow = length(true_classes),
                ncol = length(all_cols),
                dimnames = list(true = true_classes, predicted = all_cols))

  for (i in seq_along(true_vec)) {
    tc <- true_vec[i]
    pc <- pred_vec[i]
    if (is.na(tc) || !(tc %in% true_classes)) next
    if (is.na(pc)) next
    if (!(pc %in% all_cols)) next
    mat[as.character(tc), as.character(pc)] <- mat[as.character(tc), as.character(pc)] + 1L
  }
  as.data.frame(mat)
}

# ── Run analysis for one benchmark × level × match-type ───────────────────────
run_analysis <- function(data, true_col, match_type) {
  match_fn <- if (match_type == "exact" & true_col == "true_serotype") {
    exact_match_serotype
  } else if (match_type == "exact" & true_col == "true_serogroup") {
    exact_match_serogroup
  } else {
    within_match
  } 
  classes  <- sort(unique(data[[true_col]][!is.na(data[[true_col]])]))

  metrics_list   <- list()
  confusion_list <- list()

  for (tool in tools) {
    true_vec <- data[[true_col]]
    pred_vec <- data[[tool]]

    metrics_list[[tool]] <- compute_metrics(true_vec, pred_vec, match_fn) %>%
      mutate(tool = tool, .before = 1)

    confusion_list[[tool]] <- build_confusion(true_vec, pred_vec, classes)
  }

  list(
    metrics   = bind_rows(metrics_list),
    confusion = confusion_list
  )
}

# ── Execute per benchmark ──────────────────────────────────────────────────────
benchmarks <- unique(wide$benchmark)

analysis_types <- list(
  serotype_exact   = list(col = "true_serotype",  match = "exact"),
  serotype_within  = list(col = "true_serotype",  match = "within"),
  serogroup_exact = list(col = "true_serogroup", match = "exact")
)

all_metrics <- list()

for (bm in benchmarks) {
  bm_data  <- filter(wide, benchmark == bm)
  bm_label <- gsub("[^A-Za-z0-9]", "_", bm)

  for (nm in names(analysis_types)) {
    cfg      <- analysis_types[[nm]]
    result   <- run_analysis(bm_data, cfg$col, cfg$match)

    # Collect metrics
    all_metrics[[paste(bm_label, nm, sep = "__")]] <- result$metrics %>%
      mutate(benchmark = bm, analysis = nm, .before = 1)

    # Save per-tool confusion matrices
    for (tool in tools) {
      conf     <- result$confusion[[tool]]
      tool_lbl <- gsub("[^A-Za-z0-9]", "_", tool)
      out      <- file.path(data_root, "confusion_matrices",
                            paste0("confusion_", bm_label, "_", nm, "_", tool_lbl, ".csv"))
      write_csv(tibble(true_class = rownames(conf), conf), out)
      message("Saved: ", out)
    }
  }
}

# ── Save combined metrics ──────────────────────────────────────────────────────
metrics_combined <- bind_rows(all_metrics)
out_metrics <- file.path(data_root, "metrics_all.csv")
write_csv(metrics_combined, out_metrics)
message("Saved: ", out_metrics)

message("\nDone.")
