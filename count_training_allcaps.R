library(dplyr)
library(readr)

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")

# ── Load and process ground truth (one row per sample) ────────────────────────
gt <- read_csv(file.path(data_root, "GPS_benchmark", "ground_truth.csv"),
               show_col_types = FALSE) %>%
  select(sample_id, Serotype) %>%
  mutate(
    Serotype = trimws(sub("(?i)serogroup\\s*", "", Serotype, perl = TRUE)),
    Serotype = sub("^0+([0-9])", "\\1", Serotype)
  ) %>%
  filter(!grepl("/", Serotype, fixed = TRUE)) %>%
  distinct(sample_id, .keep_all = TRUE)

# ── Load ALLCAPS predictions and get the predicted sample_ids ─────────────────
allcaps <- read_csv(file.path(data_root, "GPS_benchmark", "allcaps_predictions.csv"),
                    show_col_types = FALSE) %>%
  filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
  filter(!grepl("/", sample_id, fixed = TRUE)) %>%
  filter(!grepl("serotype", sample_id, ignore.case = TRUE)) %>%
  mutate(sample_id = sub("#.*$", "", sample_id)) %>%
  group_by(sample_id) %>%
  slice_max(order_by = serotype_confidence, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(sample_id)

# ── Samples in ground truth but NOT in ALLCAPS predictions ────────────────────
missing <- gt %>%
  anti_join(allcaps, by = "sample_id")

# ── Count per serotype ─────────────────────────────────────────────────────────
counts <- missing %>%
  count(Serotype, name = "n_genomes_without_allcaps") %>%
  arrange(desc(n_genomes_without_allcaps)) %>%
  filter(Serotype != "NON-CBL")

message("Genomes in ground truth without ALLCAPS predictions: ", nrow(missing))
print(counts, n = Inf)

# ── Save ───────────────────────────────────────────────────────────────────────
out <- file.path(data_root, "GPS_benchmark_training_allcaps_counts.csv")
write_csv(counts, out)
message("Saved to: ", out)
