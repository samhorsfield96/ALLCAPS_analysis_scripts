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

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
metrics   <- read.csv(file.path(data_root, "knn-sweep", "k_sweep_per_fold.csv"))
ground_truth <- read.csv(file.path(data_root, "GPS_benchmark", "ground_truth.csv"))
plot_dir <- file.path(data_root, "plots")
dir.create(plot_dir, showWarnings = FALSE)
create_merged_knn_file <- FALSE

allcaps_serotypes <- read.csv(
  file.path(dirname(rstudioapi::getSourceEditorContext()$path),
            "data", "ALLCAPS_possible_serotypes.csv")
)$Serotypes

# read in AUROC metrics
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
palette_seed <- 42
set.seed(palette_seed)
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

#create figure showing classification accuracy of different novel and non-novel samples as well as average distance for novel and non-novel
process_knn_raw <- function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  
  file_str_list <- strsplit(file, "/")[[1]]
  loo_serotype <- file_str_list[length(file_str_list) - 1]
  
  sample_id_new <- str_split_fixed(df$sample_id, "#", 2)
  df$sample_id <- sample_id_new[,1]
  df$Contig_ID <- sample_id_new[,2]
  
  df$Contig_ID <- as.character(df$Contig_ID)
  df$nn_genogroup <- as.character(df$nn_genogroup)
  
  df$loo_serotype <- loo_serotype
  df$nn_serotype <- as.character(df$nn_serotype)
  
  df <- df %>%
    mutate(
      loo_serotype = trimws(sub("(?i)serogroup\\s*", "", loo_serotype, perl = TRUE)),
      loo_serotype = sub("^0+([0-9])", "\\1", loo_serotype),
      loo_serogroup = sub("^([0-9]+).*", "\\1", loo_serotype),
      nn_serotype = trimws(sub("(?i)serogroup\\s*", "", nn_serotype, perl = TRUE)),
      nn_serotype =sub("^0+([0-9])", "\\1", nn_serotype),
      nn_serogroup = sub("^([0-9]+).*", "\\1", nn_serotype),
      nn_genogroup = trimws(sub("(?i)serogroup\\s*", "", nn_genogroup, perl = TRUE))
    )
  
  df <- df |>
    select(sample_id, Contig_ID, loo_serotype, loo_serogroup, knn_distance, nn_distance, nn_serotype, nn_serogroup, nn_genogroup, is_novel)
  
  df
}

# create merged file of knn distances at K=1
if (create_merged_knn_file) {
  files <- list.files(file.path(data_root, "knn-sweep", "knn_raw"), pattern = "knn_query_distances.csv", recursive = TRUE, full.names = TRUE)
  all_results <- lapply(files, process_knn_raw)
  combined_query <- bind_rows(all_results)
  combined_query$is_held_out <- TRUE
  
  files <- list.files(file.path(data_root, "knn-sweep", "knn_raw"), pattern = "knn_id_distances.csv", recursive = TRUE, full.names = TRUE)
  all_results <- lapply(files, process_knn_raw)
  combined_training <- bind_rows(all_results)
  combined_training$is_held_out <- FALSE
  
  combined <- rbind(combined_query, combined_training)
  
  ground_truth$Contig_ID <- as.character(ground_truth$Contig_ID)

  # merge with serotype information
  combined_merged <- left_join(combined, ground_truth, by = c("sample_id", "Contig_ID"))
  combined_merged$Is_capsule <- ifelse(combined_merged$Is_capsule == 1, TRUE, FALSE)
  combined_merged <- combined_merged %>%
    mutate(
      Serotype = trimws(sub("(?i)serogroup\\s*", "", Serotype, perl = TRUE)),
      Serotype = sub("^0+([0-9])", "\\1", Serotype),
      Serogroup = sub("^([0-9]+).*", "\\1", Serotype),
    )
  
  # determine with genogroup assignments
  genogroups <- combined %>%
    distinct(nn_serotype, nn_serogroup, nn_genogroup)
  colnames(genogroups) <- c("Serotype", "Serogroup", "Genogroup")
  
  # merge genogroups
  final_knn_df <- combined_merged %>%
    left_join(
      genogroups %>%
        group_by(Serotype) %>%
        slice(1) %>%
        ungroup() %>%
        select(Serotype, Genogroup),
      by = "Serotype"
    )
  
  write_csv(final_knn_df, file.path(data_root, "all_1_nn_results.csv"))
  
} else {
  final_knn_df <- read_csv(file.path(data_root, "all_1_nn_results.csv"), show_col_types = FALSE)
}

# filter combined merged
final_knn_df <- final_knn_df %>% filter(loo_serotype %in% allcaps_serotypes)

# helper for calculating per serotype metrics
match_serotype  <- function(vec, cls) !is.na(vec) & vec == cls
compute_metrics <- function(final_knn_df) {
  classes <- sort(unique(final_knn_df$loo_serotype))
  map_dfr(classes, function(cls) {
    
    class_rows <- final_knn_df[match_serotype(final_knn_df$loo_serotype, cls),]
    
    TP <- sum(class_rows$is_held_out & class_rows$is_novel, na.rm = TRUE)
    FP <- sum(!class_rows$is_held_out & class_rows$is_novel, na.rm = TRUE)
    FN <- sum( class_rows$is_held_out & !class_rows$is_novel, na.rm = TRUE)
    TN <- sum(!class_rows$is_held_out & !class_rows$is_novel, na.rm = TRUE)
    Total <- TP + FN
    
    f1 <- F1_Score(y_true = class_rows$is_held_out, y_pred = class_rows$is_novel, positive = TRUE)
    sensitivity <- Sensitivity(y_true = class_rows$is_held_out, y_pred = class_rows$is_novel, positive = TRUE)
    specificity <- Specificity(y_true = class_rows$is_held_out, y_pred = class_rows$is_novel, positive = TRUE)
    accuracy <- Accuracy(y_true = class_rows$is_held_out, y_pred = class_rows$is_novel)
    precision <- Precision(y_true = class_rows$is_held_out, y_pred = class_rows$is_novel, positive = TRUE)
    
    row <- tibble(Serotype = as.character(cls), TP, FP, FN, TN, Total, sensitivity, specificity, precision, accuracy, f1)
    row <- row %>% replace(is.na(.), 0.0)
    row
  })
}

knn_metrics <- compute_metrics(final_knn_df)

write_csv(knn_metrics, file.path(data_root, "all_1_nn_metrics.csv"))




