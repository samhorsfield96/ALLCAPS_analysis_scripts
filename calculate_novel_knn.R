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

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
ground_truth <- read.csv(file.path(data_root, "merged_ground_truth.csv"))
colnames(ground_truth) <- c("sample_id", "Serotype", "Serogroup", "dataset", "benchmark")

# sample serotypes predicted by ALLCAPS
allcaps_serotypes <- read.csv(
  file.path(dirname(rstudioapi::getSourceEditorContext()$path),
            "data", "ALLCAPS_possible_serotypes.csv")
)$Serotypes

# assign out put directories
plot_dir <- file.path(data_root, "plots")
dir.create(plot_dir, showWarnings = FALSE)
create_merged_knn_file <- FALSE

#create file containing all knns
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
    select(sample_id, Contig_ID, loo_serotype, loo_serogroup, nn_distance, nn_serotype, nn_serogroup, nn_genogroup, is_novel)
  
  df
}

# create merged file of knn distances at K=1
files <- list.files(file.path(data_root, "knn-sweep", "knn_raw"), pattern = "knn_query_distances.csv", recursive = TRUE, full.names = TRUE)
all_results <- lapply(files, process_knn_raw)
combined_query <- bind_rows(all_results)
combined_query$is_held_out <- TRUE

files <- list.files(file.path(data_root, "knn-sweep", "knn_raw"), pattern = "knn_id_distances.csv", recursive = TRUE, full.names = TRUE)
all_results <- lapply(files, process_knn_raw)
combined_training <- bind_rows(all_results)
combined_training$is_held_out <- FALSE

combined <- rbind(combined_query, combined_training)

# merge with serotype information
combined_merged <- left_join(combined, ground_truth, by = c("sample_id"))
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

# reset is novel based distances
final_knn_df <- final_knn_df %>%
  select(-is_novel)

# filter combined merged
final_knn_df <- final_knn_df %>% filter(loo_serotype %in% allcaps_serotypes)

# generate ROC curves per serotype
roc_results <- final_knn_df %>%
  group_by(nn_serotype) %>%
  filter(
    n_distinct(is_held_out) == 2,
    !is.na(nn_distance)
  ) %>%
  nest() %>%
  mutate(
    roc = map(
      data,
      ~ roc(
        response = .x$is_held_out,
        predictor = .x$nn_distance,
        levels = c(FALSE, TRUE),
        direction = "<",
        quiet = TRUE
      )
    ),
    auc = map_dbl(roc, auc)
  )

# pick best threshold per serotype
best_thresholds <- roc_results %>%
  mutate(
    best = purrr::map(
      roc,
      ~ pROC::coords(
        .x,
        x = "best",
        best.method = "youden",
        ret = c("threshold", "sensitivity", "specificity"),
        transpose = FALSE
      )
    )
  ) %>%
  select(nn_serotype, auc, best) %>%
  tidyr::unnest(best)

write_csv(best_thresholds, file.path(data_root, "1_nn_best_thresholds_per_serotype.csv"))

# generate ROC curves for all
roc_result_all <- final_knn_df %>%
  nest() %>%
  mutate(
    roc = map(
      data,
      ~ roc(
        response = .x$is_held_out,
        predictor = .x$nn_distance,
        levels = c(FALSE, TRUE),
        direction = "<",
        quiet = TRUE
      )
    ),
    auc = map_dbl(roc, auc)
  )

# pick best threshold for all
best_threshold_all <- roc_result_all %>%
  mutate(
    best = purrr::map(
      roc,
      ~ pROC::coords(
        .x,
        x = "best",
        best.method = "youden",
        ret = c("threshold", "sensitivity", "specificity"),
        transpose = FALSE
      )
    )
  ) %>%
  select(auc, best) %>%
  tidyr::unnest(best)

write_csv(best_thresholds, file.path(data_root, "1_nn_best_thresholds_all.csv"))

# classify each serotype based on ROC
final_knn_df <- final_knn_df %>%
  left_join(
    best_thresholds %>%
      select(nn_serotype, threshold),
    by = "nn_serotype"
  ) %>%
  mutate(
    is_novel_threshold_per_serotype = if_else(
      is.na(threshold),
      FALSE,
      nn_distance > threshold
    )
  )

# repeat using the cut-off for all
final_knn_df$is_novel_threshold_all <- ifelse(final_knn_df$nn_distance > best_threshold_all$threshold, TRUE, FALSE)

write_csv(final_knn_df, file.path(data_root, "all_1_nn_results.csv"))

# calculating per serotype novel classification accuracy
match_serotype  <- function(vec, cls) !is.na(vec) & vec == cls
compute_metrics <- function(final_knn_df, pred_col_name, f1_beta) {
  classes <- sort(unique(final_knn_df$loo_serotype))
  map_dfr(classes, function(cls) {
    class_rows <- final_knn_df[match_serotype(final_knn_df$loo_serotype, cls),]
    true_col <- as.logical(class_rows$is_held_out)
    pred_col <- as.logical(class_rows[[pred_col_name]])
    
    TP <- sum(true_col & pred_col, na.rm = TRUE)
    FP <- sum(!true_col & pred_col, na.rm = TRUE)
    FN <- sum( true_col & !pred_col, na.rm = TRUE)
    TN <- sum(!true_col & !pred_col, na.rm = TRUE)
    Total <- TP + FN
    
    f1 <- F1_Score(y_true = true_col, y_pred = pred_col, positive = TRUE)
    sensitivity <- Sensitivity(y_true = true_col, y_pred = pred_col, positive = TRUE)
    specificity <- Specificity(y_true = true_col, y_pred = pred_col, positive = TRUE)
    accuracy <- Accuracy(y_true = true_col, y_pred = pred_col)
    precision <- Precision(y_true = true_col, y_pred = pred_col, positive = TRUE)
    
    weighted_f1 <- (1 + f1_beta^2) * (precision * sensitivity) /
      (f1_beta^2 * precision + sensitivity)
    
    row <- tibble(Serotype = as.character(cls), TP, FP, FN, TN, Total, sensitivity, specificity, precision, accuracy, f1, weighted_f1)
    row <- row %>% replace(is.na(.), 0.0)
    row
  })
}

compute_distances <- function(final_knn_df, pred_col_name) {
  classes <- sort(unique(final_knn_df$loo_serotype))
  map_dfr(classes, function(cls) {
    
    class_rows <- final_knn_df[match_serotype(final_knn_df$loo_serotype, cls),]
    true_col <- as.logical(class_rows$is_held_out)
    pred_col <- as.logical(class_rows[[pred_col_name]])
    
    distances_TP <- class_rows[true_col & pred_col,]
    distances_FP <- class_rows[!true_col & pred_col,]
    distances_FN <- class_rows[true_col & !pred_col,]
    distances_TN <- class_rows[!true_col & !pred_col,]
    
    distances_held_out <- class_rows[true_col,]
    distances_training <- class_rows[!true_col,]
    
    quantiles_distances_TP <-  quantile(distances_TP$nn_distance, c(0.25, 0.5, 0.75))
    quantiles_distances_FP <-  quantile(distances_FP$nn_distance, c(0.25, 0.5, 0.75))
    quantiles_distances_FN <-  quantile(distances_FN$nn_distance, c(0.25, 0.5, 0.75))
    quantiles_distances_TN <-  quantile(distances_TN$nn_distance, c(0.25, 0.5, 0.75))
    
    quantiles_distances_held_out <- quantile(distances_held_out$nn_distance, c(0.25, 0.5, 0.75))
    quantiles_distances_training <- quantile(distances_training$nn_distance, c(0.25, 0.5, 0.75))
    
    row <- tibble(Serotype = as.character(cls),  
                  LQ_distances_TP = quantiles_distances_TP[1], median_distances_TP = quantiles_distances_TP[2], UQ_distances_TP = quantiles_distances_TP[3],
                  LQ_distances_FP = quantiles_distances_FP[1], median_distances_FP = quantiles_distances_FP[2], UQ_distances_FP = quantiles_distances_FP[3],
                  LQ_distances_TN = quantiles_distances_TN[1], median_distances_TN = quantiles_distances_TN[2], UQ_distances_TN = quantiles_distances_TN[3],
                  LQ_distances_FN = quantiles_distances_FN[1], median_distances_FN = quantiles_distances_FN[2], UQ_distances_FN = quantiles_distances_FN[3],
                  LQ_distances_held_out = quantiles_distances_held_out[1], median_distances_held_out = quantiles_distances_held_out[2], UQ_distances_held_out = quantiles_distances_held_out[3],
                  LQ_distances_training = quantiles_distances_training[1], median_distances_training = quantiles_distances_training[2], UQ_distances_training = quantiles_distances_training[3])
    row <- row %>% replace(is.na(.), 0.0)
    row
  })
}

# determine classification accuracy
f1_beta <- 4
knn_metrics_per_serotype <- compute_metrics(final_knn_df, "is_novel_threshold_per_serotype", f1_beta)
knn_metrics_all <- compute_metrics(final_knn_df, "is_novel_threshold_all", f1_beta)
write_csv(knn_metrics_per_serotype, file.path(data_root, "all_1_nn_metrics_threshold_per_serotype.csv"))
write_csv(knn_metrics_all, file.path(data_root, "all_1_nn_metrics_threshold_all.csv"))

# determine distances of TPs etc
knn_distances_per_serotype <- compute_distances(final_knn_df, "is_novel_threshold_per_serotype")
knn_distances_all <- compute_distances(final_knn_df, "is_novel_threshold_all")

# per serotype
plot_data <- knn_distances_per_serotype %>%
  select(Serotype,
         LQ_distances_TP,
         median_distances_TP,
         UQ_distances_TP,
         LQ_distances_FP,
         median_distances_FP,
         UQ_distances_FP,
         LQ_distances_TN,
         median_distances_TN,
         UQ_distances_TN,
         LQ_distances_FN,
         median_distances_FN,
         UQ_distances_FN) %>%
  pivot_longer(
    cols = c(
      LQ_distances_TP,
      median_distances_TP,
      UQ_distances_TP,
      LQ_distances_FP,
      median_distances_FP,
      UQ_distances_FP,
      LQ_distances_TN,
      median_distances_TN,
      UQ_distances_TN,
      LQ_distances_FN,
      median_distances_FN,
      UQ_distances_FN
    ),
    names_to = c("statistic", "group"),
    names_pattern = "^(LQ|median|UQ)_distances_(TP|FP|TN|FN)$"
  ) %>%
  pivot_wider(
    names_from = statistic,
    values_from = value
  ) %>%
  rename(
    ymin = LQ,
    middle = median,
    ymax = UQ
  )

#plot_data$group[plot_data$group == "held_out"] <- "Novel"
#plot_data$group[plot_data$group == "training"] <- "Training"

plot_data <- plot_data %>%
  mutate(
    serotype_number = as.numeric(sub("^([0-9]+).*", "\\1", Serotype)),
    serotype_suffix = sub("^[0-9]+", "", Serotype)
  ) %>%
  arrange(serotype_number, serotype_suffix) %>%
  mutate(
    Serotype = factor(Serotype, levels = unique(Serotype))
  ) %>%
  select(-serotype_number, -serotype_suffix)

plot_data$group <- factor(plot_data$group, levels = c("TP", "FP", "TN", "FN"))

# generate distance plots showing distances between NNs for held-out and known serotypes
p.distance <- ggplot(
  plot_data,
  aes(
    x = Serotype,
    y = middle,
    ymin = ymin,
    ymax = ymax,
    group = group
  )
) +
  geom_errorbar(
    aes(color = group),
    position = position_dodge(width = 0.5),
    width = 0.15
  ) +
  geom_point(
    aes(color = group),
    position = position_dodge(width = 0.5),
    size = 2
  ) +
  coord_flip() +
  labs(
    x = "Serotype",
    y = "1st NN Cosine Distance",
    color = "Classification"
  ) +
  theme_light() +
  scale_color_npg()
p.distance

ggsave(file.path(plot_dir, "knn_distance_quartiles_threshold_per_serotype.png"), plot=p.distance, width = 9, height = 11)
ggsave(file.path(plot_dir, "knn_distance_quartiles_threshold_per_serotype.pdf"), plot=p.distance, width = 9, height = 11)

# all threshold
plot_data <- knn_distances_all %>%
  select(Serotype,
         LQ_distances_TP,
         median_distances_TP,
         UQ_distances_TP,
         LQ_distances_FP,
         median_distances_FP,
         UQ_distances_FP,
         LQ_distances_TN,
         median_distances_TN,
         UQ_distances_TN,
         LQ_distances_FN,
         median_distances_FN,
         UQ_distances_FN) %>%
  pivot_longer(
    cols = c(
      LQ_distances_TP,
      median_distances_TP,
      UQ_distances_TP,
      LQ_distances_FP,
      median_distances_FP,
      UQ_distances_FP,
      LQ_distances_TN,
      median_distances_TN,
      UQ_distances_TN,
      LQ_distances_FN,
      median_distances_FN,
      UQ_distances_FN
    ),
    names_to = c("statistic", "group"),
    names_pattern = "^(LQ|median|UQ)_distances_(TP|FP|TN|FN)$"
  ) %>%
  pivot_wider(
    names_from = statistic,
    values_from = value
  ) %>%
  rename(
    ymin = LQ,
    middle = median,
    ymax = UQ
  )

#plot_data$group[plot_data$group == "held_out"] <- "Novel"
#plot_data$group[plot_data$group == "training"] <- "Training"

plot_data <- plot_data %>%
  mutate(
    serotype_number = as.numeric(sub("^([0-9]+).*", "\\1", Serotype)),
    serotype_suffix = sub("^[0-9]+", "", Serotype)
  ) %>%
  arrange(serotype_number, serotype_suffix) %>%
  mutate(
    Serotype = factor(Serotype, levels = unique(Serotype))
  ) %>%
  select(-serotype_number, -serotype_suffix)

plot_data$group <- factor(plot_data$group, levels = c("TP", "FP", "TN", "FN"))

# generate distance plots showing distances between NNs for held-out and known serotypes
p.distance <- ggplot(
  plot_data,
  aes(
    x = Serotype,
    y = middle,
    ymin = ymin,
    ymax = ymax,
    group = group
  )
) +
  geom_errorbar(
    aes(color = group),
    position = position_dodge(width = 0.5),
    width = 0.15
  ) +
  geom_point(
    aes(color = group),
    position = position_dodge(width = 0.5),
    size = 2
  ) +
  coord_flip() +
  labs(
    x = "Serotype",
    y = "1st NN Cosine Distance",
    color = "Sample Type"
  ) +
  theme_light() +
  scale_color_npg()
p.distance

ggsave(file.path(plot_dir, "knn_distance_quartiles_threshold_all.png"), plot=p.distance, width = 9, height = 11)
ggsave(file.path(plot_dir, "knn_distance_quartiles_threshold_all.pdf"), plot=p.distance, width = 9, height = 11)

