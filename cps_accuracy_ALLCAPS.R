library(MLmetrics)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(patchwork)

# ── ALLCAPS-typeable serotypes ─────────────────────────────────────────────────
data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
allcaps_serotypes <- read_csv(
  file.path(data_root, "ALLCAPS_possible_serotypes.csv"),
  show_col_types = FALSE
)$Serotypes

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
  
  # --- Prediction files ---
  pred_files <- setdiff(
    list.files(dir_path, pattern = "\\.csv$", full.names = TRUE),
    gt_file
  )
  
  pred_list <- lapply(pred_files, function(f) {
    tool_name <- tools::file_path_sans_ext(basename(f))
    message("  Reading: ", basename(f))
    
    if (grepl("allcaps", basename(f), ignore.case = TRUE)) {
      df <- read_csv(f, show_col_types = FALSE) %>%
        mutate(
          True_cps = !grepl("NONCBL#", sample_id, fixed = TRUE),
          sample_id = sample_id |>
            sub("^NONCBL#", "", x = _) |>
            sub("#.*$", "", x = _)
        ) %>%
        group_by(sample_id) %>%
        slice_max(order_by = serotype_confidence, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(
          sample_id,
          predicted_serotype = pred_argmax,
          is_cbl,
          is_novel_energy,
          True_cps
        ) %>%
        mutate(tool = "ALLCAPS")
    } else {
      return(NULL)
    }
    
    df
  })
  
  preds <- bind_rows(compact(pred_list)) %>%
    mutate(predicted_serotype = sub("^0+([0-9])", "\\1", predicted_serotype))
  
  preds <- preds %>%
    rename(Predicted_cps = is_cbl)
  
  # --- Merge predictions with ground truth ---
  merged <- ground_truth %>%
    right_join(preds, by = "sample_id", relationship = "one-to-many")
  
  merged$benchmark <- basename(dir_path)
  
  # # only include samples with sample calls
  # merged <- merged %>%
  #   filter(true_serotype %in% allcaps_serotypes)
  
  merged
}

exact_match  <- function(vec, classes_vec, cls) vec[classes_vec == cls]

compute_metrics <- function(true_vec, pred_vec, classes_vec) {
  classes  <- sort(unique(classes_vec[!is.na(classes_vec)]))
  
  map_dfr(classes, function(cls) {
    actual_pos    <- exact_match(true_vec, classes_vec, cls)
    predicted_pos <- exact_match(pred_vec, classes_vec, cls)
    
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
    
    row <- tibble(Serotype = as.character(cls), TP, FP, FN, TN, Total, sensitivity, specificity, precision, accuracy, f1)
    row <- row %>% replace(is.na(.), 0.0)
    row
  })
}

run_analysis <- function(data) {
  true_col <- "True_cps"
  predicted_col <- "Predicted_cps"
  classes_col <- "true_serotype"
  
  true_vec <- data[[true_col]]
  pred_vec <- data[[predicted_col]]
  classes_vec <- data[[classes_col]]
  tool <- unique(data$tool)
  
  metrics <- compute_metrics(true_vec, pred_vec, classes_vec)
  metrics$tool <- tool
  metrics
}

benchmark_dirs <- list.dirs(data_root, recursive = FALSE, full.names = TRUE)

all_results <- lapply(benchmark_dirs, process_benchmark_dir)
all_results <- Filter(Negate(is.null), all_results)

combined <- bind_rows(all_results)
benchmarks <- unique(combined$benchmark)

all_metrics <- list()

for (bm in benchmarks) {
  bm_data  <- filter(combined, benchmark == bm)
  bm_label <- gsub("[^A-Za-z0-9]", "_", bm)
  
  result  <- run_analysis(bm_data)
  
  # Collect metrics
  all_metrics[[bm_label]] <- result %>%
    mutate(benchmark = bm, .before = 1)
}

# ── Save combined metrics ──────────────────────────────────────────────────────
metrics_combined <- bind_rows(all_metrics)
out_metrics <- file.path(data_root, "cps_accuracy_ALLCAPS.csv")
write_csv(metrics_combined, out_metrics)
message("Saved: ", out_metrics)

# Plot results
plot_dir <- file.path(data_root, "plots")
dir.create(plot_dir, showWarnings = FALSE)

# Metrics to plot
metric_cols <- c("sensitivity", "specificity", "precision", "f1")
metric_labels <- c(
  sensitivity = "Sensitivity",
  specificity = "Specificity",
  precision   = "Precision",
  f1          = "F1 Score"
)

# ── Helper: build one boxplot for a single benchmark ─────────────────────────
make_boxplot <- function(df, bm_label) {
  df_long <- df %>%
    select(tool, Serotype, all_of(metric_cols)) %>%
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
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),
      strip.text   = element_text(face = "bold", size = 14),
      axis.title.y = element_text(size = 16),
      plot.title   = element_text(face = "bold", size = 14),
      legend.position = "none"
    )
}

# only us GPS benchmark
benchmarks <- c("GPS_benchmark")
plots <- lapply(benchmarks, function(bm) {
  df <- metrics_combined %>% filter(benchmark == bm)
  if (nrow(df) == 0) return(NULL)
  make_boxplot(df, bm)
})
plots <- Filter(Negate(is.null), plots)
if (length(plots) == 0) next

# plots <- lapply(benchmarks, function(bm) {
#   df <- metrics_combined %>% filter(benchmark == bm)
#   if (nrow(df) == 0) return(NULL)
#   make_boxplot(df, bm)
# })
# plots <- Filter(Negate(is.null), plots)
# if (length(plots) == 0) next


combined_plot <- wrap_plots(plots, ncol = 1) +
  #plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 24, face = 'bold'))
out_file <- file.path(plot_dir, paste0("cps_accuracy_ALLCAPS.pdf"))
ggsave(out_file, combined_plot, width = 10, height = 6)
out_file <- file.path(plot_dir, paste0("cps_accuracy_ALLCAPS.png"))
ggsave(out_file, combined_plot, width = 10, height = 6)
message("Saved: ", out_file)

message("\nDone.")

