library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(Polychrome)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggsci)
library(ggpubr)
library(MLmetrics)
library(pROC)
library(purrr)
library(ggrepel)
library(tidyverse)

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
# assign out put directories
plot_dir <- file.path(data_root, "plots")
dir.create(plot_dir, showWarnings = FALSE)

#TODO plot AUROCs across k-nn values
files_threshold_per_serotype <- list.files(file.path(data_root), pattern = "*_nn_best_thresholds_per_serotype.csv", recursive = TRUE, full.names = TRUE)
files_threshold_all <- list.files(file.path(data_root), pattern = "*_nn_best_thresholds_all.csv", recursive = TRUE, full.names = TRUE)

file_training_testing <- file.path(data_root, "nn_training_testing_data.csv")
df_training_testing <- read_csv(file_training_testing, show_col_types = FALSE)
training_counts <- df_training_testing %>%
  filter(dataset == "Training" & benchmark == "GPS benchmark") %>%
  group_by(nn_serotype) %>%
  summarise(sum = n())

training_distances <- df_training_testing %>%
  filter(dataset == "Training" & benchmark == "GPS benchmark") %>%
  group_by(nn_serotype, is_held_out) %>%
  summarise(average_nn_dist = mean(knn_distance))

i <- 1
for (file in files_threshold_per_serotype) {
  df <- read_csv(file, show_col_types = FALSE)
  
  if (i == 1) {
    df_threshold_per_serotype <- df
  } else {
    df_threshold_per_serotype <- rbind(df_threshold_per_serotype, df)
  }
  i <- i + 1
}

i <- 1
for (file in files_threshold_all) {
  df <- read_csv(file, show_col_types = FALSE)
  
  if (i == 1) {
    df_threshold_all <- df
  } else {
    df_threshold_all <- rbind(df_threshold_all, df)
  }
  i <- i + 1
}

# generate plot showing best K for each serotype
p.box <- ggplot(
  df_threshold_per_serotype,
  aes(
    x = as.factor(k),
    y = auc,
    colour = as.factor(k)
  )
) +
  geom_boxplot(outliers = TRUE) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme_light() +
  labs(
    x = "K-value",
    y = "AUROC"
  ) +
  # Mean point
  stat_summary(
    fun = median,
    geom = "point",
    shape = 18,
    size = 3,
    colour = "black"
  ) +
  
  # Mean value as text
  stat_summary(
    fun = median,
    geom = "text",
    aes(label = sprintf("%.3f", after_stat(y))),
    vjust = 2.0,
    colour = "black",
    size = 3.5
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_color_npg() +
  scale_y_continuous(limits = c(0, 1.0)) +
  geom_point(
    data = df_threshold_all,
    aes(x = factor(k), y = auc),
    shape = 4,
    size = 3,
    stroke = 1.2
  )
p.box

ggsave(file.path(plot_dir, "k_vs_AUROC_boxplot.png"), plot=p.box, width = 9, height = 6)
ggsave(file.path(plot_dir, "k_vs_AUROC_boxplot.pdf"), plot=p.box, width = 9, height = 6)


#TODO plot accuracy values across k-nn values
files_metrics_per_serotype <- list.files(file.path(data_root), pattern = "*_nn_metrics_threshold_per_serotype.csv", recursive = TRUE, full.names = TRUE)
files_metrics_all <- list.files(file.path(data_root), pattern = "*_nn_metrics_threshold_all.csv", recursive = TRUE, full.names = TRUE)

i <- 1
for (file in files_metrics_per_serotype) {
  df <- read_csv(file, show_col_types = FALSE)
  k_val <- unique(df$k)
  
  if (i == 1) {
    df_metrics_per_serotype <- df
  } else {
    df_metrics_per_serotype <- rbind(df_metrics_per_serotype, df)
  }
  i <- i + 1
}

i <- 1
for (file in files_metrics_all) {
  df <- read_csv(file, show_col_types = FALSE)
  
  if (i == 1) {
    df_metrics_all <- df
  } else {
    df_metrics_all <- rbind(df_metrics_all, df)
  }
  i <- i + 1
}

# Metrics to plot
metric_cols <- c("sensitivity", "specificity", "precision", "f1")
metric_labels <- c(
  sensitivity = "Sensitivity",
  specificity = "Specificity",
  precision   = "Precision",
  f1          = "F1 Score"
)

# per serotype results
df_metrics_per_serotype_long <- df_metrics_per_serotype %>%
  select(Serotype, k, all_of(metric_cols)) %>%
  pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metric_cols))

df_metrics_per_serotype_long$k <- factor(
  df_metrics_per_serotype_long$k,
  levels = sort(unique(df_metrics_per_serotype_long$k))
)
balanced <- df_metrics_per_serotype_long %>%
  group_by(k, metric) %>%
  summarise(balanced = mean(value, na.rm = TRUE), .groups = "drop")

ymax_per <- df_metrics_per_serotype_long %>%
  group_by(k, metric) %>%
  summarise(y_pos = max(value, na.rm = TRUE), .groups = "drop")

annot <- balanced %>%
  left_join(ymax_per, by = c("k", "metric")) %>%
  mutate(label = sprintf("%.3f", balanced), y_annot = pmin(y_pos + 0.04, 1.05))

p.df_metrics_per_serotype <- ggplot(df_metrics_per_serotype_long, aes(x = k, y = value, fill = k)) +
  geom_boxplot(outlier.size = 0.8, na.rm = TRUE) +
  geom_point(data = annot, aes(x = k, y = balanced),
             shape = 23, size = 3, fill = "white", colour = "black",
             inherit.aes = FALSE) +
  geom_text(data = annot, aes(x = k, y = y_annot, label = label),
            size = 4, vjust = 0, inherit.aes = FALSE) +
  facet_wrap(~metric, ncol = 2, labeller = labeller(metric = metric_labels)) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.2)) +
  scale_fill_npg() +
  labs(x = "K-value", y = "Score") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    strip.text   = element_text(face = "bold", size = 14),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.title   = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

ggsave(file.path(plot_dir, "nn_sensitivity_threshold_per_serotype.pdf"), plot=p.df_metrics_per_serotype, width = 9, height = 11)
ggsave(file.path(plot_dir, "nn_sensitivity_threshold_per_serotype.png"), plot=p.df_metrics_per_serotype, width = 9, height = 11)

# per serotype results
df_metrics_all_long <- df_metrics_all %>%
  select(Serotype, k, all_of(metric_cols)) %>%
  pivot_longer(all_of(metric_cols), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = metric_cols))

df_metrics_all_long$k <- factor(
  df_metrics_all_long$k,
  levels = sort(unique(df_metrics_all_long$k))
)

balanced <- df_metrics_all_long %>%
  group_by(k, metric) %>%
  summarise(balanced = mean(value, na.rm = TRUE), .groups = "drop")

ymax_per <- df_metrics_all_long %>%
  group_by(k, metric) %>%
  summarise(y_pos = max(value, na.rm = TRUE), .groups = "drop")

annot <- balanced %>%
  left_join(ymax_per, by = c("k", "metric")) %>%
  mutate(label = sprintf("%.3f", balanced), y_annot = pmin(y_pos + 0.04, 1.05))

p.df_metrics_all <- ggplot(df_metrics_all_long, aes(x = k, y = value, fill = k)) +
  geom_boxplot(outlier.size = 0.8, na.rm = TRUE) +
  geom_point(data = annot, aes(x = k, y = balanced),
             shape = 23, size = 3, fill = "white", colour = "black",
             inherit.aes = FALSE) +
  geom_text(data = annot, aes(x = k, y = y_annot, label = label),
            size = 4, vjust = 0, inherit.aes = FALSE) +
  facet_wrap(~metric, ncol = 2, labeller = labeller(metric = metric_labels)) +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.2)) +
  scale_fill_npg() +
  labs(x = "K-value", y = "Score") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 40, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    strip.text   = element_text(face = "bold", size = 14),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    plot.title   = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

ggsave(file.path(plot_dir, "nn_sensitivity_threshold_all.pdf"), plot=p.df_metrics_all, width = 9, height = 11)
ggsave(file.path(plot_dir, "nn_sensitivity_threshold_all.png"), plot=p.df_metrics_all, width = 9, height = 11)

# TODO plot sensitivity vs specificity for each LOO serotype scatter for each k-value
# also include distributions of distances for each serotype, 

k_vals <- sort(unique(df_metrics_per_serotype_long$k))

for (k_val in k_vals) {
  k_df <- subset(df_metrics_per_serotype, k == k_val & Total > 0)
  # Order serotypes by sensitivity (highest -> lowest)
  serotype_order <- k_df %>%
    arrange(sensitivity) %>%
    pull(Serotype)
  
  plot_data <- k_df %>%
    select(Serotype, sensitivity, specificity) %>%
    pivot_longer(
      cols = c(sensitivity, specificity),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Serotype = factor(Serotype, levels = serotype_order)
    )
  
  plot_data$Metric[plot_data$Metric == "sensitivity"] <- "Sensitivity"
  plot_data$Metric[plot_data$Metric == "specificity"] <- "Specificity"
  
  p.acc <- ggplot(plot_data, aes(x = Serotype, y = Value, group = Serotype)) +
    geom_line(
      aes(group = Serotype),
      colour = "grey60",
      linewidth = 0.6
    ) +
    geom_point(
      aes(colour = Metric),
      size = 3
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1)
    ) +
    labs(
      x = "Serotype",
      y = "Statistic value",
      colour = NULL
    ) +
    theme_light() +
    coord_flip() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 16),
      axis.title = element_text(face = "bold", size = 16),
      legend.position = "right"
    ) + scale_colour_npg()
  p.acc
  #TODO add number of genomes in training to each point, next to the serotype labels on the y axis
  
  ggsave(file.path(plot_dir, paste0(k_val, "_accuracy_threshold_per_serotype.pdf")), plot=p.acc, width = 7, height = 12)
  ggsave(file.path(plot_dir, paste0(k_val, "_accuracy_threshold_per_serotype.png")), plot=p.acc, width = 7, height = 12)
  
  # plot sensitivity, specificity vs sample size
  plot_data <- merge(k_df, training_counts, by.x = "Serotype", by.y = "nn_serotype")
  
  plot_data <- plot_data %>%
    select(Serotype, sensitivity, specificity, sum) %>%
    pivot_longer(
      cols = c(sensitivity, specificity),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Serotype = factor(Serotype, levels = serotype_order)
    )
  
  plot_data$Metric[plot_data$Metric == "sensitivity"] <- "Sensitivity"
  plot_data$Metric[plot_data$Metric == "specificity"] <- "Specificity"
  
  p.acc.sum <- ggplot(plot_data, aes(x = sum, y = Value, group = Serotype)) +
    geom_point(
      size = 3
    ) +
    labs(
      x = "Total training samples",
      y = "Statistic value",
      colour = NULL
    ) +
    facet_grid(. ~ Metric) +
    theme_light() +
    scale_x_log10() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 16),
      axis.title = element_text(face = "bold", size = 16),
      legend.position = "right"
    ) + scale_colour_npg()
    p.acc.sum
  ggsave(file.path(plot_dir, paste0(k_val, "_accuracy_vs_training_size_per_serotype.pdf")), plot=p.acc.sum, width = 6, height = 4)
  ggsave(file.path(plot_dir, paste0(k_val, "_accuracy_vs_training_size_per_serotype.png")), plot=p.acc.sum, width = 6, height = 4)
    
  # plot sensitivity, specificity vs average held in/held out distance
  plot_data <- merge(k_df, training_distances, by.x = "Serotype", by.y = "nn_serotype")
  
  plot_data <- plot_data %>%
    select(Serotype, sensitivity, specificity, average_nn_dist, is_held_out) %>%
    pivot_longer(
      cols = c(sensitivity, specificity),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Serotype = factor(Serotype, levels = serotype_order)
    )
  
  plot_data$Metric[plot_data$Metric == "sensitivity"] <- "Sensitivity"
  plot_data$Metric[plot_data$Metric == "specificity"] <- "Specificity"
  plot_data$is_held_out[plot_data$is_held_out == TRUE] <- "Held-out"
  plot_data$is_held_out[plot_data$is_held_out == FALSE] <- "Held-in"
  
  p.acc.dist <- ggplot(plot_data, aes(x = average_nn_dist, y = Value, group = Serotype)) +
    geom_point(
      size = 3
    ) +
    labs(
      x = "Average NN distance",
      y = "Statistic value",
      colour = NULL
    ) +
    facet_grid(is_held_out ~ Metric) +
    theme_light() +
    scale_x_log10() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 16),
      axis.title = element_text(face = "bold", size = 16),
      legend.position = "right"
    ) + scale_colour_npg()
  p.acc.dist
  
  ggsave(file.path(plot_dir, paste0(k_val, "_accuracy_vs_distance_per_serotype.pdf")), plot=p.acc.dist, width = 6, height = 6)
  ggsave(file.path(plot_dir, paste0(k_val, "_accuracy_vs_distance_per_serotype.png")), plot=p.acc.dist, width = 6, height = 6)
  
  # plot sensitivity, specificity vs average held in/held out distance ratio
  training_distances_ratio <- training_distances %>%
    pivot_wider(
      names_from = is_held_out,
      values_from = average_nn_dist,
      names_prefix = "held_out_"
    ) %>%
    mutate(
      ratio = held_out_TRUE / held_out_FALSE
    )
    
  plot_data <- merge(k_df, training_distances_ratio, by.x = "Serotype", by.y = "nn_serotype")
  
  plot_data <- plot_data %>%
    select(Serotype, sensitivity, specificity, ratio) %>%
    pivot_longer(
      cols = c(sensitivity, specificity),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(
      Serotype = factor(Serotype, levels = serotype_order)
    )
  
  plot_data$Metric[plot_data$Metric == "sensitivity"] <- "Sensitivity"
  plot_data$Metric[plot_data$Metric == "specificity"] <- "Specificity"
  
  p.acc.dist.ratio <- ggplot(plot_data, aes(x = ratio, y = Value, group = Serotype)) +
    geom_point(
      size = 3
    ) +
    labs(
      x = "Average NN distance ratio",
      y = "Statistic value",
      colour = NULL
    ) +
    facet_grid(.~ Metric) +
    theme_light() +
    scale_x_log10() +
    theme(
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 16),
      axis.title = element_text(face = "bold", size = 16),
      legend.position = "right"
    ) + scale_colour_npg()
  p.acc.dist.ratio
  
  ggsave(file.path(plot_dir, paste0(k_val, "_accuracy_vs_distance_ratio_per_serotype.pdf")), plot=p.acc.dist.ratio, width = 6, height = 4)
  ggsave(file.path(plot_dir, paste0(k_val, "_accuracy_vs_distance_ratio_size_per_serotype.png")), plot=p.acc.dist.ratio, width = 6, height = 4)
    
}

#TODO reclassify ATB data, need to determine which is nearest neighbour and quote accuracy of nearest neighbour assignment per serotype