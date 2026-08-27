library(phytools)
library(ape)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
knn.file <- file.path(data_root, "ATB_knn_query_distances.csv")
sample.file <- file.path(data_root, "ATB_query_results.csv")

# combine nn and sample prediction dfs
sample.df <- read_csv(sample.file, show_col_types = FALSE) %>%
  filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
  mutate(Contig_ID = sub("^.*#", "", sample_id),
        sample_id = sub("#.*$", "", sample_id)) %>%
  select(
    sample_id,
    Contig_ID,
    predicted_serotype = pred_argmax,
    serotype_confidence,
    is_cbl,
    is_novel_serogroup,
    is_novel_energy
  ) %>%
  mutate(
    predicted_serogroup = sub("^([0-9]+).*", "\\1", predicted_serotype)
  )

knn.df <- read_csv(knn.file, show_col_types = FALSE) %>%
  filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
  mutate(Contig_ID = sub("^.*#", "", sample_id),
         sample_id = sub("#.*$", "", sample_id)) %>%
  select(
    sample_id,
    Contig_ID,
    nn_distance,
    nn_serotype,
    nn_genogroup
  ) %>%
  mutate(
    nn_serogroup = sub("^([0-9]+).*", "\\1", nn_serotype)
  )

combined.df <- inner_join(sample.df, knn.df, by = c("sample_id", "Contig_ID"))

# read in best thresholds
best_thresholds <- read_csv(file.path(data_root, "1_nn_best_thresholds_per_serotype.csv"), show_col_types = FALSE)

# classify each serotype based on ROC
combined.df <- combined.df %>%
  left_join(
    best_thresholds %>%
      select(nn_serotype, threshold),
    by = "nn_serotype"
  ) %>%
  mutate(
    is_novel_serogroup = if_else(
      is.na(threshold),
      FALSE,
      nn_distance > threshold
    )
  )

combined.df$matching_nn_pred <- combined.df$predicted_serotype == combined.df$predicted_serotype

write_csv(combined.df, file.path(data_root, "ATB_query_results_knn.csv"))
