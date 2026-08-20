library(dplyr)
library(readr)
library(tidyr)
library(Polychrome)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggsci)
library(ggpubr)

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
metrics   <- read.csv(file.path(data_root, "knn-sweep", "k_sweep_per_fold.csv"))
plot_dir <- file.path(data_root, "plots")
dir.create(plot_dir, showWarnings = FALSE)

allcaps_serotypes <- read.csv(
  file.path(dirname(rstudioapi::getSourceEditorContext()$path),
            "data", "ALLCAPS_possible_serotypes.csv")
)$Serotypes


metrics <- metrics %>%
  rename(serotype = fold) %>%
  mutate(
    serotype = trimws(sub("(?i)serogroup\\s*", "", serotype, perl = TRUE)),
    serotype = sub("^0+([0-9])", "\\1", serotype),
    serogroup = sub("^([0-9]+).*", "\\1", serotype)
  )  %>%
  # remove ambiguous true serotype calls
  filter(serotype %in% allcaps_serotypes)

# Define ordered serogroups
serogroups <- sort(as.numeric(unique(as.character(metrics$serogroup))))

# Generate colours
cols <- createPalette(
  length(serogroups),
  seedcolors = c("#000000", "#E41A1C", "#377EB8", "#4DAF4A")
)
# Name colours by factor levels
names(cols) <- serogroups

# pick which k to plot
k_vals <- sort(unique(metrics$k))

p_list <- list()

for (i in seq_along(k_vals)) {
  
  k_val <- k_vals[i]
  k_df <- subset(metrics, k == k_val)
  
  # order by auroc
  k_df <- k_df[order(-k_df$auroc), ]
  k_df$serotype <- factor(k_df$serotype, levels = k_df$serotype)
  
  k_df$bar_colour <- ifelse(seq_len(nrow(k_df)) %% 2 == 0, "darkgrey", "lightgrey")
  
  p <- ggplot(k_df, aes(x = serotype, y = auroc, fill = bar_colour)) +
    geom_col() +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    scale_fill_identity() +
    labs(
      x = "Serotype",
      y = "AUROC",
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 16, face="bold"),
      legend.position = "right"
    ) + coord_cartesian(ylim = c(0.4, 1.0))
  
  ggsave(paste0(file.path(plot_dir, "knn_"), k_val, ".pdf"), plot=p, width = 18, height = 6)
  ggsave(paste0(file.path(plot_dir, "knn_"), k_val, ".png"), plot=p, width = 18, height = 6)
}

# Ensure k is ordered numerically
metrics <- metrics %>%
  arrange(
    as.numeric(as.character(serogroup)),
    as.numeric(str_extract(serotype, "^\\d+")),
    str_extract(serotype, "[A-Za-z]+$")
  )

metrics$serotype <- factor(
  metrics$serotype,
  levels = unique(metrics$serotype)
)

p.all <- ggplot(
  metrics,
  aes(
    x = k,
    y = auroc,
    colour = serogroup
  )
) +
  geom_point(size = 3) +
  geom_line() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  facet_wrap(~ serotype) +   # replace facet_variable with your column
  theme_light() +
  labs(
    x = "K-value",
    y = "AUROC"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) + 
  scale_x_log10() +
  scale_colour_manual(
    values = cols,
    breaks = as.character(serogroups)
  ) +
  scale_y_continuous(limits = c(0, 1.0))
p.all

ggsave(file.path(plot_dir, "knn_all.png"), plot=p.all, width = 18, height = 12)
ggsave(file.path(plot_dir, "knn_all.pdf"), plot=p.all, width = 18, height = 12)

# generate plot showing best K for each serotype
p.box <- ggplot(
  metrics,
  aes(
    x = as.factor(k),
    y = auroc,
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
  scale_y_continuous(limits = c(0, 1.0))
p.box

ggsave(file.path(plot_dir, "k_vs_AUROC_boxplot.png"), plot=p.box, width = 9, height = 6)
ggsave(file.path(plot_dir, "k_vs_AUROC_boxplot.pdf"), plot=p.box, width = 9, height = 6)

# TODO combine all knn-raw data into single table, then determine assignment to model serogroups rather than serotypes
# analyse closest nearest neighbour
nn_metrics   <- read.csv(file.path(data_root, "knn-sweep", "knn_nn_summary.csv"))

nn_metrics <- nn_metrics %>%
  rename(serotype = held_out_serotype) %>%
  mutate(
    serotype = trimws(sub("(?i)serogroup\\s*", "", serotype, perl = TRUE)),
    serotype = sub("^0+([0-9])", "\\1", serotype),
    serogroup = sub("^([0-9]+).*", "\\1", serotype),
    modal_nn_serotype = trimws(sub("(?i)serogroup\\s*", "", modal_nn_serotype, perl = TRUE)),
    modal_nn_serotype =sub("^0+([0-9])", "\\1", modal_nn_serotype),
    modal_nn_serogroup = sub("^([0-9]+).*", "\\1", modal_nn_serotype),
  )  %>%
  # remove ambiguous true serotype calls
  filter(serotype %in% allcaps_serotypes) |>
  select(serotype, serogroup, n_query, modal_nn_serotype, modal_nn_serogroup, modal_nn_serotype_frac,
         median_nn_distance)
nn_metrics$within_serogroup = ifelse((nn_metrics$serotype != nn_metrics$serogroup), TRUE, FALSE)

# order by modal frequency
nn_metrics <- nn_metrics[order(-nn_metrics$modal_nn_serotype_frac), ]
nn_metrics$serotype <- factor(nn_metrics$serotype, levels = nn_metrics$serotype)

nn_metrics$bar_colour <- ifelse(nn_metrics$within_serogroup, "#CC0000", "lightgrey")

p <- ggplot(nn_metrics, aes(x = serotype, y = modal_nn_serotype_frac, fill = bar_colour)) +
  geom_col() +
  scale_fill_identity() +
  labs(
    x = "Serotype",
    y = "Prop. serotypes assigned",
  ) +
  geom_text(
    aes(label = modal_nn_serogroup),
    hjust = -0.1,
    angle = 90,
    size = 3
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 16, face="bold"),
    legend.position = "right"
  ) + coord_cartesian(ylim = c(0.0, 1.0))
p

ggsave(paste0(file.path(plot_dir, "novel_nn_assignment", ".pdf")), plot=p, width = 18, height = 6)
ggsave(paste0(file.path(plot_dir, "novel_nn_assignment", ".png")), plot=p, width = 18, height = 6)
