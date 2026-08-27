library(phytools)
library(ape)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ggtree", ask = FALSE)
BiocManager::install("ggtreeExtra", ask = FALSE)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(Polychrome)
library(ggpubr)
library(patchwork)

# set to false if not generating new intermediate files
WRITE_NEW_INTERMEDIATE <- FALSE

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")
sample.file <- file.path(data_root, "ATB_query_results.csv")

generate_tree <- function(tree, dataframe){
  # remove tips with no data
  tips_with_data <- dataframe$tip   # or meta$tip
  new.tree <- drop.tip(
    tree,
    setdiff(tree$tip.label, tips_with_data)
  )
  
  # --- midpoint root the tree ---
  new.tree <- midpoint.root(new.tree)
  
  # -------------------------------
  # Reorder taxonomy to tree tips
  # -------------------------------
  new.dataframe <- dataframe[dataframe$tip %in% new.tree$tip.label,]
  
  # -------------------------------
  # HARD SAFETY CHECKS
  # -------------------------------
  stopifnot(
    nrow(new.dataframe) == length(new.tree$tip.label),
    all(rownames(new.dataframe$tip) == new.tree$tip.label)
  )
  
  list(new.tree, new.dataframe)
}

tree.file <- file.path(data_root, "pneumo_ATB_tree.nwk")
tree <- read.tree(tree.file)

if (WRITE_NEW_INTERMEDIATE == TRUE) {
  # ensure pick contig with that is non-novel, is cbl and then has highest serotype confidence
  sample.df <- read_csv(sample.file, show_col_types = FALSE) %>%
    filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
    mutate(sample_id = sub("#.*$", "", sample_id)) %>%
    group_by(sample_id) %>%
    arrange(
      desc(is_cbl & !is_novel_serogroup),
      desc(serotype_confidence)
    ) %>%
    slice(1) %>%
    ungroup() %>%
    select(
      sample_id,
      predicted_serotype = pred_argmax,
      is_cbl,
      is_novel_serogroup,
      is_novel_energy
    ) %>%
    mutate(
      tool = "ALLCAPS",
      predicted_serogroup = sub("^([0-9]+).*", "\\1", predicted_serotype)
    )
  colnames(sample.df) <- c("tip", "predicted_serotype", "is_cbl", "is_novel_serogroup", "is_novel_energy", "tool", "predicted_serogroup")

  
  #only run if downsampling tree
  tree_list <- generate_tree(tree, sample.df)
  tree  <- tree_list[[1]]
  sample.df    <- tree_list[[2]]
  write.tree(tree, file = "pneumo_ATB_tree.nwk")
  
  poppunk.file <- file.path(data_root, "ATB_PopPUNK_results.csv")
  poppunk.df <- read_csv(poppunk.file, show_col_types = FALSE)
  
  merged.df <- merge(sample.df, poppunk.df, by.x = "tip", by.y = "sample", all.x = TRUE)
  
  merged.df$final_serotype_prediction <- merged.df$predicted_serotype
  merged.df$final_serotype_prediction[merged.df$is_cbl == "TRUE" & merged.df$is_novel_serogroup == "TRUE"] <- "Novel"
  merged.df$final_serotype_prediction[merged.df$is_cbl == "FALSE"] <- "Non-cps"
  
  merged.df$final_serogroup_prediction <- merged.df$predicted_serogroup
  merged.df$final_serogroup_prediction[merged.df$is_cbl == "TRUE" & merged.df$is_novel_serogroup == "TRUE"] <- "Novel"
  merged.df$final_serogroup_prediction[merged.df$is_cbl == "FALSE"] <- "Non-cps"
  
  # order and factor
  merged.df <- merged.df[order(merged.df$predicted_serotype),]
  merged.df$predicted_serotype <- factor(merged.df$predicted_serotype)
  merged.df <- merged.df[order(merged.df$predicted_serogroup),]
  merged.df$predicted_serogroup <- factor(merged.df$predicted_serogroup)
  merged.df <- merged.df[order(merged.df$final_serotype_prediction),]
  merged.df$final_serotype_prediction <- factor(merged.df$final_serotype_prediction)
  merged.df <- merged.df[order(merged.df$final_serogroup_prediction),]
  merged.df$final_serogroup_prediction <- factor(merged.df$final_serogroup_prediction)
  
  seroba_meta <- read.csv(file.path(data_root, "seroba_meta.tsv"), sep = "\t", header = FALSE)
  colnames(seroba_meta) <- c("Name", "Genogroup", "Associated_Serotype", "SeroBA_Unknown", "SeroBA_Type")
  seroba_meta[] <- lapply(seroba_meta, function(x) {
    if (is.character(x)) sub("^0+", "", x) else x
  })
  
  merged.df <- merge(merged.df, seroba_meta, by.x = "predicted_serotype", by.y = "Associated_Serotype", all.x = TRUE, all.y = FALSE) 
  
  write.csv(merged.df, file.path(data_root, "parsed_ATB_query_results.csv"), row.names = FALSE, quote = FALSE)
} else {
  merged.file <- file.path(data_root, "parsed_ATB_query_results.csv")
  merged.df <- read_csv(merged.file, show_col_types = FALSE)
}

# generate GPSC proportion data
# minimum clade to plot
min.GPSC.size <- 300

proportion_data <- merged.df %>%
  group_by(GPSC, predicted_serogroup, final_serogroup_prediction) %>%
  tally() %>%
  group_by(GPSC) %>%
  mutate(Proportion = n / sum(n)) %>%
  filter(!is.na(GPSC),  GPSC != "nan", !str_detect(GPSC, ";"), final_serogroup_prediction != "Non-cps") %>%
  ungroup() %>%
  group_by(GPSC) %>%
  mutate(GPSC_total = sum(n)) %>%
  filter(GPSC_total > min.GPSC.size) %>%
  ungroup() %>%
  mutate(GPSC = forcats::fct_reorder(GPSC, GPSC_total, .desc = TRUE)) %>%
  arrange(desc(GPSC_total), GPSC)

proportion_data$predicted_serogroup <- as.character(proportion_data$predicted_serogroup)

# Define ordered serogroups and colours for plots
serogroups <- as.character(sort(as.numeric(unique(as.character(proportion_data$predicted_serogroup)))))
palette_seed <- 42

set.seed(palette_seed)
# Generate colours
cols <- createPalette(
  length(serogroups),
  seedcolors = c("#000000", "#E41A1C", "#377EB8", "#4DAF4A")
)

# Name colours by factor levels
names(cols) <- as.character(serogroups)
cols <- c(cols, Novel="grey")

top_proportion_data <- proportion_data %>%
  filter(final_serogroup_prediction != "Novel") %>%
  group_by(GPSC) %>%
  # Keeps the entire row where Proportion is highest for that GPSC
  #slice_max(order_by = Proportion, n = 1, with_ties = FALSE) %>%
  ungroup() 

# Find GPSCs present in top_proportion_data but missing from top_proportion_data_novel
missing_GPSC <- proportion_data %>%
  anti_join(top_proportion_data, by = "GPSC")

# Add placeholder rows
placeholder_rows <- missing_GPSC %>%
  mutate(
    final_serogroup_prediction = "NA",
    Proportion = 0,
    n = NA
  )

# Add them to the novel data
top_proportion_data <- bind_rows(
  top_proportion_data,
  placeholder_rows
)

top_proportion_data_novel <- proportion_data %>%
  filter(final_serogroup_prediction == "Novel") %>%
  group_by(GPSC) %>%
  ungroup() 

# Find GPSCs present in top_proportion_data but missing from top_proportion_data_novel
missing_GPSC <- proportion_data %>%
  anti_join(top_proportion_data_novel, by = "GPSC")

# Add placeholder rows
placeholder_rows <- missing_GPSC %>%
  mutate(
    final_serogroup_prediction = "Novel",
    Proportion = 0,
    n = NA
  )

# Add them to the novel data
top_proportion_data_novel <- bind_rows(
  top_proportion_data_novel,
  placeholder_rows
)

# --- 1. YOUR DATA DOWNSAMPLING LOGIC ---
# Calculate serotype proportions per GPSC

# Identify one representative tip per GPSC
representatives <- inner_join(top_proportion_data, merged.df, by = "GPSC") %>%
  arrange(GPSC, tip) %>% 
  group_by(GPSC) %>%
  summarise(tip = first(tip), .groups = "drop")

# Downsample the tree using your custom function
tree_list <- generate_tree(tree, representatives)
collapsed_tree  <- tree_list[[1]]

# Rename the remaining tips to their corresponding GPSC name
name_bridge <- match(collapsed_tree$tip.label, representatives$tip)
collapsed_tree$tip.label <- representatives$GPSC[name_bridge]

write.tree(collapsed_tree, file.path(data_root, "GPSC_tree.nwk"))

# Format data frame to match the new tip labels (GPSC names) for the facet plot
plot_data <- top_proportion_data %>%
  rename(label = GPSC)
write.csv(plot_data, file.path(data_root, "per_GPSC_data_observed.csv"), row.names = FALSE, quote = FALSE)
#add empty rows to enable correct colour plotting
plot_data <- plot_data %>%
  tidyr::complete(
    label,
    final_serogroup_prediction = names(cols),
    predicted_serogroup = names(cols),
    fill = list(n = 0)
  )

# Plot tree for known serotypes
p.observed <- ggtree(
  collapsed_tree,
  layout = "circular",
  branch.length = "none"
)

# Set factor order
plot_data$predicted_serogroup <- factor(
  plot_data$predicted_serogroup,
  levels = names(cols)
)

# Proportion bars
p.observed <- p.observed +
  geom_fruit(
    data = plot_data,
    geom = geom_col,
    mapping = aes(
      y = label,
      x = Proportion,
      fill = predicted_serogroup
    ),
    orientation = "y",
    pwidth = 0.18,
    offset = 0.02
  ) +
  scale_fill_manual(values = cols,
                    limits = names(cols),
                    breaks = names(cols),
                    drop = FALSE)

# Reference ring at proportion = 1
p.observed <- p.observed +
  geom_fruit(
    data = plot_data,
    geom = geom_tile,
    mapping = aes(
      y = label,
      x = 0,
      width = 0.01,
      height = 1
    ),
    fill = "black",
    inherit.aes = FALSE,
    pwidth = 0.18,
    offset = 0.00
  ) +   
  geom_tiplab(
    aes(label = label),
    offset = 4.0,
    align = FALSE
  ) + labs(fill = "Predicted Serogroup")

p.observed
ggsave(file.path(data_root, "plots", "ATB_tree_plot.pdf"), width = 14, height = 10, plot = p.observed)

# generate same tree for novel genomes
plot_data <- top_proportion_data_novel %>%
  rename(label = GPSC)
write.csv(plot_data, file.path(data_root, "per_GPSC_data_novel.csv"), row.names = FALSE, quote = FALSE)
#add empty rows to enable correct colour plotting
plot_data <- plot_data %>%
  tidyr::complete(
    label,
    final_serogroup_prediction = names(cols),
    predicted_serogroup = names(cols),
    fill = list(n = 0)
  )

# --- 2. THE LINEAR FACET PLOT (TDbook Architecture) ---

p.novel <- ggtree(
  collapsed_tree,
  layout = "circular",
  branch.length = "none"
)

# Set factor order
plot_data$predicted_serogroup <- factor(
  plot_data$predicted_serogroup,
  levels = names(cols)
)

# Proportion bars
p.novel <- p.novel +
  geom_fruit(
    data = plot_data,
    geom = geom_col,
    mapping = aes(
      y = label,
      x = Proportion,
      fill = predicted_serogroup
    ),
    orientation = "y",
    pwidth = 0.18,
    offset = 0.02
  ) +
  scale_fill_manual(values = cols,
                    limits = names(cols),
                    breaks = names(cols),
                    drop = FALSE)

# Reference ring at proportion = 1
p.novel <- p.novel +
  geom_fruit(
    data = plot_data,
    geom = geom_tile,
    mapping = aes(
      y = label,
      x = 0,
      width = 0.01,
      height = 1
    ),
    fill = "black",
    inherit.aes = FALSE,
    pwidth = 0.18,
    offset = 0.00
  ) +   
  geom_tiplab(
    aes(label = label),
    offset = 4.0,
    align = FALSE
  ) + labs(fill = "Predicted Serogroup")

p.novel
ggsave(file.path(data_root, "plots", "ATB_tree_plot_novel.pdf"), width = 14, height = 10, plot = p.novel)

# generate plot for counts
plot_data <- proportion_data %>%
  group_by(GPSC, final_serogroup_prediction) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  tidyr::complete(
    GPSC,
    final_serogroup_prediction = names(cols),
    fill = list(n = 0)
  ) %>%
  mutate(
    final_serogroup_prediction = factor(
      final_serogroup_prediction,
      levels = names(cols)
    )
  )

p.dist <- ggplot(plot_data, aes(x=GPSC, y=n, group=GPSC, fill = final_serogroup_prediction)) +
  geom_col() +
  scale_fill_manual(values = cols,
                    limits = names(cols),
                    breaks = names(cols),
                    drop = FALSE) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
        axis.title = element_text(size = 14)) +
  ylab("Total genomes") +
  xlab("GPSC") +
  labs(fill = "Predicted Serogroup")
p.dist

p.observed <- p.observed +
  theme(legend.position = "none")

p.novel <- p.novel +
  theme(legend.position = "none")

p.observed <- p.observed +
  theme(plot.margin = margin(5.5, 50, 5.5, 5.5))

p.novel <- p.novel +
  theme(plot.margin = margin(5.5, 5.5, 5.5, 50))

combined_plot <- (
  (p.observed | p.novel) /
    p.dist
) +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

combined_plot
ggsave(file.path(data_root, "plots", "ATB_tree_hist.pdf"), width = 15, height = 10, plot = combined_plot)
ggsave(file.path(data_root, "plots", "ATB_tree_hist.png"), width = 15, height = 10, plot = combined_plot)

# For distirbution of annotation proportions
proportion_data <- merged.df %>%
  group_by(GPSC, predicted_serogroup, final_serogroup_prediction) %>%
  tally() %>%
  group_by(GPSC) %>%
  mutate(Proportion = n / sum(n)) %>%
  filter(!is.na(GPSC), GPSC != "nan", !str_detect(GPSC, ";"), final_serogroup_prediction != "Novel", final_serogroup_prediction != "Non-cps") %>%
  ungroup() %>%
  group_by(GPSC) %>%
  mutate(GPSC_total = sum(n)) %>%
  filter(GPSC_total > min.GPSC.size) %>%
  ungroup() %>%
  mutate(GPSC = forcats::fct_reorder(GPSC, GPSC_total, .desc = TRUE)) %>%
  arrange(desc(GPSC_total), GPSC)

top_proportion_data <- proportion_data %>%
  filter(final_serogroup_prediction != "Novel") %>%
  group_by(GPSC) %>%
  # Keeps the entire row where Proportion is highest for that GPSC
  slice_max(order_by = Proportion, n = 1, with_ties = FALSE) %>%
  ungroup() 

box_stats <- top_proportion_data %>%
  summarise(
    median = median(Proportion, na.rm = TRUE),
    Q1 = quantile(Proportion, 0.25, na.rm = TRUE),
    Q3 = quantile(Proportion, 0.75, na.rm = TRUE)
  )

p.box.large <- ggplot(top_proportion_data, aes(y = Proportion)) +
  geom_boxplot() +
  geom_text(
    data = box_stats,
    aes(
      x = 0.25,
      y = 0.25,
      label = paste0(
        "Q3: ", round(Q3, 2),
        "\nMedian: ", round(median, 2),
        "\nQ1: ", round(Q1, 2)
      )
    ),
  ) +
  theme_light() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab("Majority Serotype Proportion") +
  xlab(NULL)

p.box.large

min.GPSC.size <- 1

proportion_data <- merged.df %>%
  group_by(GPSC, predicted_serogroup, final_serogroup_prediction) %>%
  tally() %>%
  group_by(GPSC) %>%
  mutate(Proportion = n / sum(n)) %>%
  filter(!is.na(GPSC), GPSC != "nan", !str_detect(GPSC, ";"), final_serogroup_prediction != "Novel", final_serogroup_prediction != "Non-cps") %>%
  ungroup() %>%
  group_by(GPSC) %>%
  mutate(GPSC_total = sum(n)) %>%
  filter(GPSC_total > min.GPSC.size) %>%
  ungroup() %>%
  mutate(GPSC = forcats::fct_reorder(GPSC, GPSC_total, .desc = TRUE)) %>%
  arrange(desc(GPSC_total), GPSC)

top_proportion_data <- proportion_data %>%
  group_by(GPSC) %>%
  # Keeps the entire row where Proportion is highest for that GPSC
  slice_max(order_by = Proportion, n = 1, with_ties = FALSE) %>%
  filter(n > min.GPSC.size) %>%
  ungroup() 

box_stats <- top_proportion_data %>%
  summarise(
    median = median(Proportion, na.rm = TRUE),
    Q1 = quantile(Proportion, 0.25, na.rm = TRUE),
    Q3 = quantile(Proportion, 0.75, na.rm = TRUE)
  )

p.box1 <- ggplot(top_proportion_data, aes(y = Proportion)) +
  geom_boxplot() +
  geom_text(
    data = box_stats,
    aes(
      x = 0.25,
      y = 0.25,
      label = paste0(
        "Q3: ", round(Q3, 2),
        "\nMedian: ", round(median, 2),
        "\nQ1: ", round(Q1, 2)
      )
    ),
  ) +
  theme_light() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab("Majority Serotype Proportion") +
  xlab(NULL)

p.box1

combined_box <- ggarrange(
  p.box.large,
  p.box1,
  labels = c("A", "B"),
  ncol = 2,
  heights = c(1, 1),
  widths = c(1, 1)
)
combined_box

ggsave(file.path(data_root, "plots", "ATB_GPSC_Serotype_proportions.png"), width = 8, height = 4, plot = combined_box)

# number of difference serotypes within each GPSC
within_GPSC_difference <- proportion_data %>%
  group_by(GPSC) %>%
  summarise(n_classes = n_distinct(predicted_serogroup))

box_stats <- within_GPSC_difference %>%
  summarise(
    median = median(n_classes, na.rm = TRUE),
    Q1 = quantile(n_classes, 0.25, na.rm = TRUE),
    Q3 = quantile(n_classes, 0.75, na.rm = TRUE)
  )

p.box.GPSC <- ggplot(within_GPSC_difference, aes(y = n_classes)) +
  geom_boxplot() +
  geom_text(
    data = box_stats,
    aes(
      x = 0.25,
      y = 10,
      label = paste0(
        "Q3: ", round(Q3, 2),
        "\nMedian: ", round(median, 2),
        "\nQ1: ", round(Q1, 2)
      )
    ),
  ) +
  theme_light() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylab("Number of serotypes per GPSC") +
  xlab(NULL)

p.box.GPSC

ggsave(file.path(data_root, "plots", "ATB_GPSC_Serotype_number.png"), width = 4, height = 4, plot = p.box.GPSC)
