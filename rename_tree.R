library(ape)
library(getopt)
library(tidyverse)

# Define the command-line arguments
spec <- matrix(c(
  'input',    'i', 1, 'character', 'Input tree file',
  'nexus',    'n', 2, 'logical',   'Input tree is nexus (default newick)',
  'csv',      'c', 2, 'character', 'Metadata CSV. Specify custom labels with -l flag, or default assumes new names in first column, old names in second',
  'tips',     't', 2, 'character', 'Column name with original tip names (default second column)',
  'label',    'l', 2, 'character', 'Custom labels: specify columns to use in labels: comma separated column names in order
                                    Default without this option is new names in first column, old names in second column',
  'output',   'o', 1, 'character', 'Output tree file',
  'renamed',  'r', 2, 'logical',   'CSV output with old and new tips names',
  'drop_old', 'd', 2, 'logical',   'Drop original tip names (default keep old name at start of new name)',
  'prefix',   'p', 2, 'character', 'Add specified string to the start of all labels'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

# setwd('/home/aileen/onedrive/treebuilding/241115carabidae/bayes/')
# opt <- data.frame(
#   input = c('backup/vas_mito_ry.nexus.con.tre'),
#   csv = c('/home/aileen/onedrive/treebuilding/241115carabidae/metadata/taxonomy.csv'),
#   label = ('tree_family,tree_subfamily,tree_tribe,tree_subtribe,tree_genus,tree_species,tree_subspecies'),
#   output = c('backup/tax.vas_mito_ry.nexus.con.tre'),
#   tips = c('tree_id'),
#   nexus = T
# )

# Function to strip quotes if present
strip_quotes <- function(x) {
  gsub("(^['\"]|['\"]$)", '', x)
}

# Load the tree and CSV data
df <- read.csv(opt$csv)
df[is.na(df)] <- ''
df[] <- lapply(df, as.character)

# Change tip name column in label string if necessary:
if (!is.null(opt$label)) {
  cat("Tree will be renamed with the following column(s) from", opt$csv, ": ", opt$label, '\n')
  if (!is.null(opt$tips)) {
    opt$label <- str_replace(opt$label, opt$tips, 'old_id')
  }
} else {
  cat("Tree will be renamed with the", colnames(df[0]), "column from", opt$csv, '\n')
}


if (is.null (opt$nexus)) {
  tree <- read.tree(opt$input)
} else {
  tree <- read.nexus(opt$input)
}

# tree <- read.nexus(opt$input)
# cat("Tree has", length(tree$tip.label), "tips and", tree$Nnode, "internal nodes\n")

#1 Get tip names
tip_names <- strip_quotes(tree$tip.label)
# Strip quotes from tree tip labels
tree$tip.label <- strip_quotes(tree$tip.label)

# CSV with pre-defined labels
if (is.null (opt$label)) {

  # Load CSV 
  df <- read.csv(opt$csv, header=FALSE)

  # Rename tips using the mapping from the CSV file
  for (i in 1:nrow(df)) {
    old_name <- strip_quotes(as.character(df[i, 2]))
    new_name <- strip_quotes(as.character(df[i, 1]))
    
    # Perform the renaming
    if (new_name != "") {
      tree$tip.label[tree$tip.label == old_name] <- new_name}
  }
} else {
    column_names <- strsplit(opt$label, ",")[[1]]
    # Assume tip names are in first colum
    if (is.null(opt$tips)) {
      df <- df %>%
        rename(old_id = names(.)[1]) %>%
        mutate(old_id = strip_quotes(old_id))

    } else {
      df <- df %>%
        rename(old_id = opt$tips) %>%
        mutate(old_id = strip_quotes(old_id))
    }
    df <- df %>%
      filter(old_id %in% tip_names)
    # Write fasta IDs
    df <- df %>%
      rowwise() %>% 
      mutate(fasta_id = ifelse(is.null(opt$drop_old),
        paste(c(old_id, c_across(all_of(column_names))), collapse = "_"), 
        paste(c(c_across(all_of(column_names))), collapse = "_"))) %>%
      ungroup()
  
    # Add prefix
    if (!is.null(opt$prefix)) {
      df$fasta_id <- paste0(opt$prefix, "_", df$fasta_id)
    }
    df <- df %>%
      mutate(fasta_id = gsub("_+$", "", fasta_id)) %>%
      select(old_id, fasta_id) %>%
      unique()
  matches <- match(tree$tip.label, df$old_id)
  tree$tip.label[!is.na(matches)] <- df$fasta_id[matches[!is.na(matches)]]
}

if (!is.null(tree$node.label)) {
  expected <- tree$Nnode
  actual <- length(tree$node.label)

  if (actual != expected) {
    warning(paste("PastML node.label length (", actual, ") does not match Nnode (", expected, "). Writing tree with original labels anyway.", sep = ""))
    
    # Optionally: truncate to expected length if needed
    tree$node.label <- tree$node.label[1:min(expected, actual)]
  }
}



# Write tree
if (is.null (opt$nexus)) {
  write.tree(tree, file=opt$output)
  cat("Tree written to", opt$output, '\n')
} else {
  write.nexus(tree, file=opt$output)
  cat("Tree written to", opt$output, '\n')
}

