library(dplyr)
library(readr)
library(tidyr)
library(Polychrome)
library(ggplot2)

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

# iterate over k values
for (k_val in k_vals) {
  k_df <- subset(metrics, k == k_val)
  
  k_df <- k_df[order(-k_df$auroc), ]
  
  k_df$serotype <- factor(
    k_df$serotype,
    levels = k_df$serotype
  )
  
  p <-ggplot(k_df, aes(x =reorder(serotype, -auroc), serotype, y = auroc, fill = serogroup)) +
    geom_col() +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x  = element_text(angle = 40, hjust = 1, size = 12),
      axis.text.y  = element_text(size = 12),
      strip.text   = element_text(face = "bold", size = 14),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.position = "right"
    ) +
    scale_fill_manual(
      values = cols,
      breaks = as.character(serogroups)  # legend in ascending order
    ) +
    xlab("Serotype") +
    ylab("AUROC") +
    labs(fill = "Serogroup") +
    geom_hline(yintercept = 0.5, linetype="dashed")
  p
  
  ggsave(paste0(file.path(plot_dir, "knn_"), k_val, ".pdf"), plot=p, width = 18, height = 6)
  ggsave(paste0(file.path(plot_dir, "knn_"), k_val, ".png"), plot=p, width = 18, height = 6)
}
