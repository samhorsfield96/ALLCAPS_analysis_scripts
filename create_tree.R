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
  new.dataframe <- dataframe[dataframe$sample_id %in% new.tree$tip.label,]
  
  # -------------------------------
  # HARD SAFETY CHECKS
  # -------------------------------
  stopifnot(
    nrow(new.dataframe) == length(new.tree$tip.label),
    all(rownames(new.dataframe$sample_id) == new.tree$tip.label)
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
colnames(sample.df) <- c("tip", "predicted_serotype", "is_cbl", "is_novel_energy","tool")

write.csv(sample.df, file.path(data_root, "parsed_ATB_query_results.csv"), row.names = FALSE, quote = FALSE)

#only run if downsampling tree
#tree_list <- generate_tree(tree, sample.df)
#tree  <- tree_list[[1]]
#sample.df    <- tree_list[[2]]
#write.tree(tree, file = "pneumo_ATB_tree.nwk")
