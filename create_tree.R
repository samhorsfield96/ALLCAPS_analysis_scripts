library(phytools)
library(ape)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("ggtree")
  BiocManager::install("ggtreeExtra")
}
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(Polychrome)
library(ggpubr)

# set to false if not generating new intermediate files
WRITE_NEW_INTERMEDIATE <- FALSE
DOWNSAMPLE <- TRUE

data_root <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "data")

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
  sample.file <- file.path(data_root, "ATB_query_results.csv")
  sample.df <- read_csv(sample.file, show_col_types = FALSE) %>%
    filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
    mutate(sample_id = sub("#.*$", "", sample_id)) %>%
    group_by(sample_id) %>%
    slice_max(order_by = serotype_confidence, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(sample_id, predicted_serotype = pred_argmax, is_cbl = is_cbl, is_novel_energy = is_novel_energy) %>%
    mutate(tool = "ALLCAPS",
           predicted_serogroup = sub("^([0-9]+).*", "\\1", predicted_serotype))
  colnames(sample.df) <- c("tip", "predicted_serotype", "is_cbl", "is_novel_energy", "tool", "predicted_serogroup")

  
  #only run if downsampling tree
  tree_list <- generate_tree(tree, sample.df)
  tree  <- tree_list[[1]]
  sample.df    <- tree_list[[2]]
  write.tree(tree, file = "pneumo_ATB_tree.nwk")
  
  poppunk.file <- file.path(data_root, "ATB_PopPUNK_results.csv")
  poppunk.df <- read_csv(poppunk.file, show_col_types = FALSE)
  
  merged.df <- merge(sample.df, poppunk.df, by.x = "tip", by.y = "sample", all.x = TRUE)
  
  write.csv(merged.df, file.path(data_root, "parsed_ATB_query_results.csv"), row.names = FALSE, quote = FALSE)
} else {
  merged.file <- file.path(data_root, "parsed_ATB_query_results.csv")
  merged.df <- read_csv(merged.file, show_col_types = FALSE)
}

# order and factor
merged.df <- merged.df[order(merged.df$predicted_serotype),]
merged.df$predicted_serotype <- factor(merged.df$predicted_serotype)
merged.df <- merged.df[order(merged.df$predicted_serogroup),]
merged.df$predicted_serogroup <- factor(merged.df$predicted_serogroup)

seroba_meta <- read.csv(file.path(data_root, "seroba_meta.tsv"), sep = "\t", header = FALSE)
colnames(seroba_meta) <- c("Name", "Genogroup", "Associated_Serotype", "Unknown", "Type")
seroba_meta[] <- lapply(seroba_meta, function(x) {
  if (is.character(x)) sub("^0+", "", x) else x
})

merged.df <- merge(merged.df, seroba_meta, by.x = "predicted_serotype", by.y = "Associated_Serotype", all.x = TRUE, all.y = FALSE) 

# generate GPSC proportion data
# minimum clade to plot
min.GPSC.size <- 100

proportion_data <- merged.df %>%
  group_by(GPSC, predicted_serogroup) %>%
  tally() %>%
  group_by(GPSC) %>%
  mutate(Proportion = n / sum(n)) %>%
  filter(!is.na(GPSC)) %>%
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


if (DOWNSAMPLE == TRUE)
{
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
  write.csv(plot_data, file.path(data_root, "per_GPSC_data.csv"), row.names = FALSE, quote = FALSE)
  
  # --- 2. THE LINEAR FACET PLOT (TDbook Architecture) ---
  
  p <- ggtree(
    collapsed_tree,
    layout = "circular",
    branch.length = "none"
  )
  
  # Define ordered serogroups
  serogroups <- sort(as.numeric(unique(as.character(plot_data$predicted_serogroup))))
  
  # Set factor order
  plot_data$predicted_serogroup <- factor(
    plot_data$predicted_serogroup,
    levels = serogroups
  )
  
  # Generate colours
  cols <- createPalette(
    length(serogroups),
    seedcolors = c("#000000", "#E41A1C", "#377EB8", "#4DAF4A")
  )
  
  # Name colours by factor levels
  names(cols) <- serogroups
  
  # Proportion bars
  p <- p +
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
    scale_fill_manual(values = cols)
  
  # Reference ring at proportion = 1
  p <- p +
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
    ) + labs(fill = "Predicted Serogroup")
  
  p
  ggsave(file.path(data_root, "plots", "ATB_tree_plot.pdf"), width = 14, height = 10, plot = p)

} else {
  
  # Data to attach to tree
  plot_data <- merged.df %>%
    rename(label = tip)   # geom_fruit matches on the tree tip labels
  
  p <- ggtree(
    tree,
    layout = "circular",
    branch.length = "none"
  )
  
  # Get the serogroups in a fixed order
  serogroups <- sort(unique(plot_data$predicted_serogroup))
  
  # Generate 48 (or however many are needed) colours
  cols <- createPalette(
    length(serogroups),
    seedcolors = c("#000000", "#E41A1C", "#377EB8", "#4DAF4A")
  )
  
  # Name the colours so ggplot matches them correctly
  names(cols) <- serogroups
  
  p <- p +
    geom_fruit(
      data = plot_data,
      geom = geom_point,
      mapping = aes(
        y = label,
        fill = predicted_serogroup,
      ),
      shape = 21,
      colour = "black",
      stroke = 0.2,
      offset = 0.02,
      pwidth = 0.08
    ) +
    scale_fill_manual(values = cols) +
    theme(
      legend.position = "right"
    )
  
  p
  
}

proportion_data$predicted_serogroup <- factor(
  proportion_data$predicted_serogroup,
  levels = serogroups
)

p.dist <- ggplot(proportion_data, aes(x=GPSC, y=n, group=GPSC, fill = predicted_serogroup)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) +
  ylab("Total genomes") +
  xlab("GPSC") 
p.dist

combined_plot <- ggarrange(
  p,
  p.dist,
  labels = c("A", "B"),
  ncol = 1,
  common.legend = TRUE,
  legend = "right",
  heights = c(1.5, 1),
  widths = c(1.5, 1)
)
combined_plot
ggsave(file.path(data_root, "plots", "ATB_tree_hist.pdf"), width = 15, height = 10, plot = combined_plot)
ggsave(file.path(data_root, "plots", "ATB_tree_hist.png"), width = 15, height = 10, plot = combined_plot)


# For distirbution of annotation proportions
box_stats <- top_proportion_data %>%
  summarise(
    median = median(Proportion, na.rm = TRUE),
    Q1 = quantile(Proportion, 0.25, na.rm = TRUE),
    Q3 = quantile(Proportion, 0.75, na.rm = TRUE)
  )

p.box100 <- ggplot(top_proportion_data, aes(y = Proportion)) +
  geom_boxplot() +
  geom_text(
    data = box_stats,
    aes(
      x = 0.25,
      y = 0.5,
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

p.box100

min.GPSC.size <- 1

proportion_data <- merged.df %>%
  group_by(GPSC, predicted_serogroup) %>%
  tally() %>%
  group_by(GPSC) %>%
  mutate(Proportion = n / sum(n)) %>%
  filter(!is.na(GPSC)) %>%
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
      y = 0.5,
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
  p.box100,
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
