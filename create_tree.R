library(phytools)
library(dplyr)
library(stringr)
library(tidyr)

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
  new.dataframe <- dataframe[new.tree$tip.label, , drop = FALSE]
  
  # -------------------------------
  # HARD SAFETY CHECKS
  # -------------------------------
  stopifnot(
    nrow(new.dataframe) == length(new.tree$tip.label),
    all(rownames(new.dataframe) == new.tree$tip.label)
  )
  
  list(new.tree, new.dataframe)
}

tree.file <- file.path(data_root, "ATB_tree_whole.nwk")
tree <- read.tree(tree.file)

sample.file <- file.path(data_root, "ATB_query_results.csv")
sample.df <- read_csv(sample.file, show_col_types = FALSE) %>%
  filter(!grepl("NONCBL#", sample_id, fixed = TRUE)) %>%
  mutate(sample_id = sub("#.*$", "", sample_id)) %>%
  group_by(sample_id) %>%
  slice_max(order_by = serotype_confidence, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(sample_id, predicted_serotype = pred_argmax, is_cbl = is_cbl, is_novel_energy = is_novel_energy) %>%
  mutate(tool = "ALLCAPS")
sample.df$tip <- sample.df$sample_id

tree_list <- generate_tree(tree, sample.df)
sub_tree  <- tree_list[[1]]
sub_df    <- tree_list[[2]]

#output new subtree
write.tree(tree, file = "my_tree.nwk")
