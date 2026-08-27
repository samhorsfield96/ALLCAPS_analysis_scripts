library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(patchwork)
library(ggrepel)

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

# ── Loop over each analysis, combining all benchmarks into one PDF ────────────
combos <- metrics %>%
  distinct(benchmark, analysis)

for (ana in unique(combos$analysis)) {
  plots <- lapply(unique(combos$benchmark), function(bm) {
    df <- metrics %>% filter(benchmark == bm, analysis == ana)
    if (nrow(df) == 0) return(NULL)
    make_boxplot(df, bm)
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) next

  combined_plot <- wrap_plots(plots, ncol = 1) +
    plot_annotation(tag_levels = "A") & 
    theme(plot.tag = element_text(size = 24, face = 'bold'))
  ana_lbl  <- gsub("[^A-Za-z0-9]", "_", ana)
  out_file <- file.path(plot_dir, paste0("boxplot_", ana_lbl, ".pdf"))
  ggsave(out_file, combined_plot, width = 14, height = 10 * length(plots))
  out_file <- file.path(plot_dir, paste0("boxplot_", ana_lbl, ".png"))
  ggsave(out_file, combined_plot, width = 14, height = 10 * length(plots))
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
      theme(
        axis.text.y  = element_text(size = 12),
        strip.text   = element_text(face = "bold", size = 14),
        axis.title.y  = element_text(size = 16),
        axis.title.x  = element_text(size = 16))

    bm_lbl  <- gsub("[^A-Za-z0-9]", "_", bm)
    ana_lbl <- gsub("[^A-Za-z0-9]", "_", ana)
    out_file <- file.path(plot_dir,
                          paste0("allcaps_vs_training_", bm_lbl, "_", ana_lbl, ".pdf"))
    ggsave(out_file, p, width = 12, height = 9)
    message("Saved: ", out_file)
  }
}

message("\nTraining-size plots done.")

# ── Per-serotype comparison: ALLCAPS vs SeroBA and ALLCAPS vs Pneumo_Typer ────
comparison_pairs <- list(
  list(tool_b = "SeroBA",       label = "SeroBA"),
  list(tool_b = "Pneumo-Typer", label = "Pneumo-Typer")
)

for (pair in comparison_pairs) {
  tool_b <- pair$tool_b
  pair_label <- pair$label

  df_pair <- metrics %>%
    filter(tool %in% c("ALLCAPS", tool_b)) %>%
    select(benchmark, analysis, tool, class, all_of(metric_cols)) %>%
    pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = metric_cols)) %>%
    pivot_wider(names_from = tool, values_from = value) %>%
    rename(ALLCAPS_val = ALLCAPS, other_val = all_of(tool_b)) %>%
    filter(!is.na(ALLCAPS_val) & !is.na(other_val))

  for (bm in unique(df_pair$benchmark)) {
    for (ana in unique(df_pair$analysis)) {
      df_plot <- df_pair %>% filter(benchmark == bm, analysis == ana)
      if (nrow(df_plot) == 0) next

      p <- ggplot(df_plot, aes(x = ALLCAPS_val, y = other_val, label = class)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
        geom_point(colour = "#CC0000", size = 1.8, alpha = 0.7) +
        geom_text_repel(size = 2, colour = "grey20") +
        facet_wrap(~metric, ncol = 2, labeller = labeller(metric = metric_labels)) +
        scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
        scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
        labs(
          x = "ALLCAPS",
          y = pair_label
        ) +
        theme_bw(base_size = 11) +
        theme(
          axis.text.x  = element_text(size = 12),
          axis.text.y  = element_text(size = 12),
          strip.text   = element_text(face = "bold", size = 14),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16)
        )

      bm_lbl  <- gsub("[^A-Za-z0-9]", "_", bm)
      ana_lbl <- gsub("[^A-Za-z0-9]", "_", ana)
      pair_lbl <- gsub("[^A-Za-z0-9]", "_", pair_label)
      out_file <- file.path(plot_dir,
                            paste0("serotype_compare_ALLCAPS_vs_", pair_lbl, "_", bm_lbl, "_", ana_lbl, ".pdf"))
      ggsave(out_file, p, width = 12, height = 9)
      out_file <- file.path(plot_dir,
                            paste0("serotype_compare_ALLCAPS_vs_", pair_lbl, "_", bm_lbl, "_", ana_lbl, ".png"))
      ggsave(out_file, p, width = 12, height = 9)
      message("Saved: ", out_file)
    }
  }
}

message("\nPer-serotype comparison plots done.")
