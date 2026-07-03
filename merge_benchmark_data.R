library(dplyr)
library(readr)
library(tidyr)

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
      # ALLCAPS: sample_id contains "#contig" — strip the hash suffix
      df <- read_csv(f, show_col_types = FALSE) %>%
        filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
        mutate(sample_id = sub("#.*$", "", sample_id)) %>%
        group_by(sample_id) %>%
        slice_max(order_by = serotype_confidence, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(sample_id, predicted_serotype = pred_argmax, is_cbl = is_cbl, is_novel_energy = is_novel_energy) %>%
        mutate(tool = "ALLCAPS")
    } else {
      # Standard format: sample_id, tool, predicted_serotype
      df <- read_csv(f, show_col_types = FALSE) %>%
        select(sample_id, tool, predicted_serotype)
      df$is_cbl = NA
      df$is_novel_energy = NA
    }
    df
  })

  preds <- bind_rows(pred_list) %>%
    mutate(predicted_serotype = sub("^0+([0-9])", "\\1", predicted_serotype))

  # --- Merge predictions with ground truth ---
  merged <- ground_truth %>%
    right_join(preds, by = "sample_id", relationship = "one-to-many")

  merged$benchmark <- basename(dir_path)
  merged
}

# --- Main ---
data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")

# Alternatively, set the path manually:
# data_root <- "/Users/samhorsfield/Software/ALLCAPS_analysis_scripts/data"

benchmark_dirs <- list.dirs(data_root, recursive = FALSE, full.names = TRUE)

all_results <- lapply(benchmark_dirs, process_benchmark_dir)
all_results <- Filter(Negate(is.null), all_results)

combined <- bind_rows(all_results)

# Keep only samples present in all tools and ground truth within each benchmark
all_tools <- unique(combined$tool)
complete_samples <- combined %>%
  group_by(benchmark, sample_id) %>%
  filter(all(all_tools %in% tool)) %>%
  ungroup()

combined <- combined %>%
  semi_join(complete_samples, by = c("benchmark", "sample_id"))

message("\nMerge complete. Dimensions: ", nrow(combined), " rows × ", ncol(combined), " columns")
print(head(combined))

# Save long-format output
out_file <- file.path(data_root, "merged_benchmark_results.csv")
write_csv(combined, out_file)
message("Saved to: ", out_file)

# --- Wide format ---
combined_wide <- combined %>%
  select(sample_id, benchmark, true_serotype, true_serogroup, tool, predicted_serotype) %>%
  pivot_wider(names_from = tool, values_from = predicted_serotype)

out_file_wide <- file.path(data_root, "merged_benchmark_results_wide.csv")
write_csv(combined_wide, out_file_wide)
message("Wide format saved to: ", out_file_wide)
