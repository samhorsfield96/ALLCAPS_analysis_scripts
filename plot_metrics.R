library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggtext)

# ── Load data ──────────────────────────────────────────────────────────────────
data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
metrics   <- read_csv(file.path(data_root, "metrics_all.csv"), show_col_types = FALSE)

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

# ── Loop over each benchmark × analysis combination ───────────────────────────
combos <- metrics %>%
  distinct(benchmark, analysis)

for (i in seq_len(nrow(combos))) {
  bm  <- combos$benchmark[i]
  ana <- combos$analysis[i]

  df <- metrics %>%
    filter(benchmark == bm, analysis == ana)

  # Long format for faceting over metric
  df_long <- df %>%
    select(tool, class, all_of(metric_cols)) %>%
    pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = metric_cols))

  # Balanced (macro-average) per tool × metric, ignoring NA
  balanced <- df_long %>%
    group_by(tool, metric) %>%
    summarise(balanced = mean(value, na.rm = TRUE), .groups = "drop")

  # Merge for annotation positioning: place label just above upper whisker
  ymax_per <- df_long %>%
    group_by(tool, metric) %>%
    summarise(y_pos = max(value, na.rm = TRUE), .groups = "drop")

  annot <- balanced %>%
    left_join(ymax_per, by = c("tool", "metric")) %>%
    mutate(
      label    = sprintf("Bal=%.3f", balanced),
      y_annot  = pmin(y_pos + 0.04, 1.05)
    )

  p <- ggplot(df_long, aes(x = tool, y = value, fill = tool)) +
    geom_boxplot(outlier.size = 0.8, na.rm = TRUE) +
    geom_text(
      data    = annot,
      aes(x = tool, y = y_annot, label = label),
      size    = 2.5,
      vjust   = 0,
      inherit.aes = FALSE
    ) +
    facet_wrap(~metric, ncol = 2, labeller = labeller(metric = metric_labels)) +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.2)) +
    scale_fill_manual(values = {
      tool_levels <- sort(unique(df_long$tool))
      greys <- grDevices::gray.colors(sum(tool_levels != "ALLCAPS"), start = 0.4, end = 0.85)
      clrs  <- setNames(character(length(tool_levels)), tool_levels)
      clrs[tool_levels != "ALLCAPS"] <- greys
      clrs["ALLCAPS"] <- "#CC0000"
      clrs
    }) +
    labs(
      x     = NULL,
      y     = "Score",
      fill  = "Tool"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1, size = 8),
      strip.text   = element_text(face = "bold"),
      legend.position = "none"
    )

  bm_lbl  <- gsub("[^A-Za-z0-9]", "_", bm)
  ana_lbl <- gsub("[^A-Za-z0-9]", "_", ana)
  out_file <- file.path(plot_dir, paste0("boxplot_", bm_lbl, "_", ana_lbl, ".pdf"))
  ggsave(out_file, p, width = 12, height = 9)
  message("Saved: ", out_file)
}

message("\nDone.")

# ── ALLCAPS accuracy vs training set size ──────────────────────────────────────
training_counts <- read_csv(
  file.path(data_root, "GPS_training_allcaps_counts.csv"),
  show_col_types = FALSE
)

allcaps_metrics <- metrics %>%
  filter(tool == "ALLCAPS") %>%
  inner_join(training_counts, by = c("class" = "Serotype")) %>%
  select(benchmark, analysis, class, n_genomes_without_allcaps,
         all_of(metric_cols)) %>%
  pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metric_cols))

for (bm in unique(allcaps_metrics$benchmark)) {
  for (ana in unique(allcaps_metrics$analysis)) {
    df_plot <- allcaps_metrics %>%
      filter(benchmark == bm, analysis == ana, !is.na(value))

    if (nrow(df_plot) == 0) next

    p <- ggplot(df_plot,
                aes(x = n_genomes_without_allcaps, y = value, label = class)) +
      geom_point(colour = "#CC0000", size = 1.8, alpha = 0.7) +
      geom_smooth(method = "loess", se = TRUE, colour = "grey30",
                  linewidth = 0.7, fill = "grey80") +
      geom_text(size = 2, vjust = -0.6, colour = "grey20") +
      facet_wrap(~metric, ncol = 2, labeller = labeller(metric = metric_labels)) +
      scale_x_log10(labels = scales::comma) +
      scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
      labs(
        x        = "Training genomes per serotype",
        y        = "Score"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.text = element_text(face = "bold"))

    bm_lbl  <- gsub("[^A-Za-z0-9]", "_", bm)
    ana_lbl <- gsub("[^A-Za-z0-9]", "_", ana)
    out_file <- file.path(plot_dir,
                          paste0("allcaps_vs_training_", bm_lbl, "_", ana_lbl, ".pdf"))
    ggsave(out_file, p, width = 12, height = 9)
    message("Saved: ", out_file)
  }
}

message("\nTraining-size plots done.")
