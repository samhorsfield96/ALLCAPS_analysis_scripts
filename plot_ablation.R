library(MLmetrics)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(patchwork)
library(stringr)

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")

# ── ALLCAPS-typeable serotypes ─────────────────────────────────────────────────
allcaps_serotypes <- read_csv(
  file.path(data_root, "ALLCAPS_possible_serotypes.csv"),
  show_col_types = FALSE
)$Serotypes

# Process a single benchmark directory and return a merged data frame
process_benchmark_dir <- function(dir_path) {
  message("Processing: ", dir_path)
  
  # --- Ground truth ---
  gt_file <- file.path(dir_path, "ground_truth.csv")
  if (!file.exists(gt_file)) {
    warning("No ground_truth.csv in ", dir_path, " — skipping.")
    return(NULL)
  }
  ground_truth <- read_csv(gt_file, show_col_types = FALSE) %>%
    select(sample_id, true_serotype = Serotype) %>%
    distinct(sample_id, .keep_all = TRUE) %>%
    mutate(
      true_serotype = trimws(sub("(?i)serogroup\\s*", "", true_serotype, perl = TRUE)),
      true_serotype = sub("^0+([0-9])", "\\1", true_serotype),
      true_serogroup = sub("^([0-9]+).*", "\\1", true_serotype)
    )
    # # remove ambiguous true serotype calls
    # filter(true_serotype %in% allcaps_serotypes)
  
  # label training or testing genomes
  ALLCAPS_pred_file <- file.path(dir_path, "allcaps_predictions.csv")
  ALLCAPS_pred_df <- read_csv(ALLCAPS_pred_file, show_col_types = FALSE) %>%
    filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
    mutate(sample_id = sub("#.*$", "", sample_id)) %>%
    group_by(sample_id) %>%
    slice_max(order_by = serotype_confidence, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    #mutate(pred_argmax = sub("^0+([0-9])", "\\1", pred_argmax)) %>%
    select(sample_id, ALLCAPS = pred_argmax)
  
  prokbert_pred_file <- file.path(dir_path, "prokbert_predictions.csv")
  prokbert_pred_df <- read_csv(prokbert_pred_file, show_col_types = FALSE) %>%
    filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
    mutate(sample_id = sub("#.*$", "", sample_id)) %>%
    group_by(sample_id) %>%
    slice_max(order_by = serotype_confidence, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    #mutate(pred_argmax = sub("^0+([0-9])", "\\1", pred_argmax)) %>%
    select(sample_id, ProkBERT = pred_argmax)

  merged <- ground_truth %>%
    right_join(ALLCAPS_pred_df, by = "sample_id", relationship = "one-to-many")
  
  merged <- merged %>% 
    left_join(prokbert_pred_df, by = "sample_id")
  
  merged$benchmark <- str_replace(basename(dir_path), "_", " ")
  
  merged
}

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
  
  for (tool in tools) {
    true_vec <- data[[true_col]]
    pred_vec <- data[[tool]]
    
    metrics_list[[tool]] <- compute_metrics(true_vec, pred_vec, match_fn) %>%
      mutate(tool = tool, .before = 1)
  }
  
  bind_rows(metrics_list)
}

dir_path <- file.path(data_root, "GPS_benchmark")
combined <- process_benchmark_dir(dir_path)
# remove ambiguous true serotype calls
combined <- combined %>%
  filter(true_serotype %in% allcaps_serotypes)

out_file <- file.path(data_root, "merged_ablation.csv")
write_csv(combined, out_file)

analysis_types <- list(
  serotype_exact   = list(col = "true_serotype",  match = "exact"),
  serogroup_exact = list(col = "true_serogroup", match = "exact")
)

tools <- c("ALLCAPS", "ProkBERT")
all_metrics <- list()
bm <- "GPS benchmark"
bm_label <- gsub("[^A-Za-z0-9]", "_", bm)

for (nm in names(analysis_types)) {
  cfg      <- analysis_types[[nm]]
  result   <- run_analysis(combined, cfg$col, cfg$match)
  
  # Collect metrics
  all_metrics[[paste(bm_label, nm, sep = "__")]] <- result %>%
    mutate(benchmark = bm, analysis = nm, .before = 1)
}

# ── Save combined metrics ──────────────────────────────────────────────────────
metrics <- bind_rows(all_metrics)

out_metrics <- file.path(data_root, "metrics_ablation_all.csv")
write_csv(metrics, out_metrics)
message("Saved: ", out_metrics)

# ── Plot combined metrics ──────────────────────────────────────────────────────
plot_dir <- file.path(data_root, "plots")
dir.create(plot_dir, showWarnings = FALSE)

# ── Helper: build one boxplot for a single benchmark ─────────────────────────
make_boxplot <- function(df, bm_label) {
  df_long <- df %>%
    select(tool, class, all_of(metric_cols)) %>%
    pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = metric_cols))
  
  balanced <- df_long %>%
    group_by(tool, metric) %>%
    summarise(balanced = mean(value, na.rm = TRUE), .groups = "drop")
  
  ymax_per <- df_long %>%
    group_by(tool, metric) %>%
    summarise(y_pos = max(value, na.rm = TRUE), .groups = "drop")
  
  annot <- balanced %>%
    left_join(ymax_per, by = c("tool", "metric")) %>%
    mutate(label = sprintf("%.3f", balanced), y_annot = pmin(y_pos + 0.04, 1.05))
  
  tool_levels <- sort(unique(df_long$tool))
  greys <- grDevices::gray.colors(sum(tool_levels != "ALLCAPS"), start = 0.4, end = 0.85)
  clrs  <- setNames(character(length(tool_levels)), tool_levels)
  clrs[tool_levels != "ALLCAPS"] <- greys
  clrs["ALLCAPS"] <- "#CC0000"
  
  ggplot(df_long, aes(x = tool, y = value, fill = tool)) +
    geom_boxplot(outlier.size = 0.8, na.rm = TRUE) +
    geom_point(data = annot, aes(x = tool, y = balanced),
               shape = 23, size = 3, fill = "white", colour = "black",
               inherit.aes = FALSE) +
    geom_text(data = annot, aes(x = tool, y = y_annot, label = label),
              size = 4, vjust = 0, inherit.aes = FALSE) +
    facet_wrap(~metric, ncol = 2, labeller = labeller(metric = metric_labels)) +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = clrs) +
    labs(x = NULL, y = "Score") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1, size = 12),
      axis.text.y  = element_text(size = 12),
      strip.text   = element_text(face = "bold", size = 14),
      axis.title.y = element_text(size = 16),
      plot.title   = element_text(face = "bold", size = 14),
      legend.position = "none"
    )
}

combos <- metrics %>%
  distinct(benchmark, analysis)

# Metrics to plot
metric_cols <- c("sensitivity", "specificity", "precision", "f1")
metric_labels <- c(
  sensitivity = "Sensitivity",
  specificity = "Specificity",
  precision   = "Precision",
  f1          = "F1 Score"
)

for (ana in unique(combos$analysis)) {
  plots <- lapply(unique(combos$benchmark), function(bm) {
    df <- metrics %>% filter(benchmark == bm, analysis == ana)
    if (nrow(df) == 0) return(NULL)
    make_boxplot(df, bm)
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) next
  
  combined_plot <- wrap_plots(plots, ncol = 1)
  ana_lbl  <- gsub("[^A-Za-z0-9]", "_", ana)
  out_file <- file.path(plot_dir, paste0("ablation_boxplot_", ana_lbl, ".pdf"))
  ggsave(out_file, combined_plot, width = 14, height = 10 * length(plots))
  message("Saved: ", out_file)
}

message("\nDone.")



