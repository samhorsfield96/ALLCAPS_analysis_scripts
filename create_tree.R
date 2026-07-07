library(phytools)
library(ape)
library(dplyr)
library(stringr)
library(tidyr)
library(readr)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("ggtreeExtra")
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

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

if (DOWNSAMPLE == TRUE)
{
  # --- 1. YOUR DATA DOWNSAMPLING LOGIC ---
  # Calculate serotype proportions per GPSC
  
  # minimum clade to plot
  min.GPSC.serotype.size <- 10
  
  proportion_data <- merged.df %>%
    group_by(GPSC, predicted_serogroup) %>%
    tally() %>%
    group_by(GPSC) %>%
    mutate(Proportion = n / sum(n)) %>%
    ungroup()
  
  top_proportion_data <- proportion_data %>%
    group_by(GPSC) %>%
    # Keeps the entire row where Proportion is highest for that GPSC
    slice_max(order_by = Proportion, n = 1, with_ties = FALSE) %>%
    filter(n > min.GPSC.serotype.size) %>%
    ungroup() 
  
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
  
  # Format data frame to match the new tip labels (GPSC names) for the facet plot
  plot_data <- top_proportion_data %>%
    rename(tip = GPSC)
  
  # --- 2. THE LINEAR FACET PLOT (TDbook Architecture) ---
  
  # Initialize a standard horizontal linear tree layout
  p <- ggtree(collapsed_tree) 
  
  # # Attach your basic strain info dataset and color the tip points
  # p <- p %<+% plot_data + 
  #   geom_tippoint(aes(color = tip)) +
  #   geom_tiplab(size = 2, offset = 0.02) # Optional: Displays the GPSC name at the tip
  
  # Reset the fill scale so your bar plot colors don't conflict with tip colors
  p <- p + new_scale_fill()
  
  # Add the horizontal stacked bar charts inside a panel facet
  p <- p + geom_facet(
    panel = "Serogroup Proportions", 
    data = plot_data, 
    geom = geom_col,               # geom_col is preferred over geom_bar(stat="identity") here
    mapping = aes(
      x = Proportion,              # The length of the bar sections (adds up to 1.0)
      fill = predicted_serogroup    # Splits the bar into stacked colors based on Serotype
    ), 
    position = position_stack(),   # Stacks the serotypes together horizontally
    orientation = 'y',             # Keeps the bars aligned horizontally with tree rows
    width = 0.6
  ) +
    scale_fill_viridis_d(option = "turbo", name = "Serogroup") +
    theme_tree2(legend.position = "right")
  
  # plot_data_test <- subset(plot_data, tip %in% unique(plot_data$tip)[1:10])
  # p.test <- ggplot(plot_data_test, aes(x = tip, y = Proportion, fill= predicted_serotype)) +
  #   geom_col()
  # p.test
  
  # Render the final horizontal plot
  p
  
} else {
  # 1. Plot the base circular tree (hide tip labels because of the large scale)
  p <- ggtree(tree, layout = "circular", linewidth = 0.1)
  
  # 2. Add the First Ring: Serotype (as a tile/heatmap ring)
  p <- p + geom_fruit(
    data = trait_data,
    geom = geom_tile,
    mapping = aes(y = TaxonID, fill = Serotype),
    offset = 0.05,        # Distance from the tree tips
    pwidth = 0.05         # Width of the ring
  ) +
    scale_fill_viridis_d(option = "plasma", name = "Serotype") # Discrete color palette
  
  # 3. Add the Second Ring: GPSC (stacked outside the first ring)
  p <- p + geom_fruit(
    data = trait_data,
    geom = geom_tile,
    mapping = aes(y = TaxonID, fill = GPSC),
    offset = 0.05,        # Distance from the Serotype ring
    pwidth = 0.05
  ) +
    # Using a new fill scale helper for a second discrete trait
    ggnewscale::new_scale_fill() + 
    scale_fill_discrete(name = "GPSC")
  
  # 4. Render the plot
  print(p)
  
}






