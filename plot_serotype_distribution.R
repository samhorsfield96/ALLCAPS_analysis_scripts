library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)

# ── ALLCAPS-typeable serotypes ────────────────────────────────────────────────
allcaps_serotypes <- read_csv(
  file.path(dirname(rstudioapi::getSourceEditorContext()$path),
            "data", "ALLCAPS_possible_serotypes.csv"),
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
    ) %>% 
    # remove ambiguous true serotype calls
    filter(true_serotype %in% allcaps_serotypes)
  
  # label training or testing genomes
  pred_file <- file.path(dir_path, "allcaps_predictions.csv")
  pred_df <- read_csv(pred_file, show_col_types = FALSE) %>%
    filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
    mutate(sample_id = sub("#.*$", "", sample_id)) %>%
    group_by(sample_id) %>%
    slice_max(order_by = serotype_confidence, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(sample_id, predicted_serotype = pred_argmax, is_cbl = is_cbl, is_novel_energy = is_novel_energy) %>%
    mutate(tool = "ALLCAPS")
  
  ground_truth$dataset <- ifelse(ground_truth$sample_id %in% pred_df$sample_id, "Testing", "Training")
  
  ground_truth$benchmark <- str_replace(basename(dir_path), "_", " ")
  ground_truth
}

# --- Main ---
data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")

benchmark_dirs <- list.dirs(data_root, recursive = FALSE, full.names = TRUE)

all_results <- lapply(benchmark_dirs, process_benchmark_dir)
all_results <- Filter(Negate(is.null), all_results)

combined <- bind_rows(all_results)

out_file <- file.path(data_root, "merged_ground_truth.csv")
write_csv(combined, out_file)

plot_df <- combined %>%
  count(true_serotype, dataset, benchmark, name = "count") %>%
  mutate(true_serotype = fct_reorder(true_serotype, count, .fun = sum, .desc = TRUE))

p <- ggplot(plot_df, aes(x = true_serotype, y = count)) +
  geom_col(fill = "steelblue") +
  geom_text(
    aes(label = count),
    hjust = -0.1,
    angle = 90,
    size = 3
  ) +
  facet_wrap(dataset ~ benchmark, ncol=1, scales = "free_y") +
  labs(
    x = "Serotype",
    y = "Count"
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.3))
  ) 
p

out_file <- file.path(data_root, "plots", "serotype_distributions")
ggsave(paste0(out_file, ".pdf"), p, width = 12, height = 14)
ggsave(paste0(out_file, ".png"), p, width = 12, height = 14)
