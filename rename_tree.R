library(ape)
library(getopt)
library(dplyr)

# Define the command-line arguments
spec <- matrix(c(
  'input',    'i', 1, 'character', 'Input tree file',
  'nexus',    'n', 2, 'logical',   'Input tree is nexus (default newick)',
  'csv',      'c', 2, 'character', 'Metadata CSV. Specify custom labels with -l flag, or default assumes new names in first column, old names in second',
  'tips',     't', 2, 'character', 'Column name with original tip names (if not first column)',
  'label',    'l', 2, 'character', 'Custom labels: specify columns to use in labels: comma separated column names in order
                                    Default without this option is new names in first column, old names in second column',
  'output',   'o', 1, 'character', 'Output tree file',
  'drop_old', 'd', 2, 'logical',   'Drop original tip names (default keep old name at start of new name)'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

# setwd('/home/aileen/onedrive/treebuilding/13Ktree/gmt/raxml')
# opt <- data.frame(
#   input = c('gmt13_coipart_nt.raxml_treeshrink/output.bestTree'),
#   csv = c('~/onedrive/treebuilding/metadata/mmg_test.csv'),
#   label = c('superfamily,family,subfamily,tribe,genus,species,morphospecies'),
#   output = c('test.tree'),
#   tips = c('lab_id')
# )

# Function to strip quotes if present
strip_quotes <- function(x) {
  gsub("(^['\"]|['\"]$)", '', x)
}

# Load the tree and CSV data
df <- read.csv(opt$csv)
df[is.na(df)] <- ''

if (is.null (opt$nexus)) {
  tree <- read.tree(opt$input)
} else {
  tree <- read.nexus(opt$input)
}

# Get tip names
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
}


# CSV with taxonomy columns
if (!is.null (opt$label)) {
  # Get specified metadata for labels
  if (!is.null (opt$label)) {
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
    if (is.null(opt$drop_old)) {
      df <- df %>%
        rowwise() %>% 
        mutate(fasta_id = paste(c(old_id, c_across(all_of(column_names))), collapse = "_")) %>% 
        ungroup()
    } else {
      df <- df %>%
        rowwise() %>% 
        mutate(fasta_id = paste(c(c_across(all_of(column_names))), collapse = "_")) %>% 
        ungroup()
    }
    df <- df %>%
      mutate(fasta_id = gsub("_+$", "", fasta_id)) %>%
      select(old_id, fasta_id) %>%
      unique()
    # Get tip names from specified column
  }
  matches <- match(tree$tip.label, df$old_id)
  tree$tip.label[!is.na(matches)] <- df$fasta_id[matches[!is.na(matches)]]
}



# Default taxonomy labels from master_metadata.csv
    # } else {
  #   meta <- meta %>%
  #     mutate(fasta_id = case_when(
  #       suborder == 'Myxophaga' | suborder == 'Archostemata' ~ paste(ncbi_taxid, suborder, family, genus, species, subspecies, sep = '_'), # Add suborder for small suborders
  #       family != "" ~ paste(ncbi_taxid, family, subfamily, tribe, genus, species, subspecies, sep = '_'),  # If family
  #       !is.na(superfamily) & superfamily != "" ~ paste(ncbi_taxid, superfamily, sep = '_'), # Otherwise use lowest available rank of higher taxonomy
  #       !is.na(infraorder) & infraorder != "" ~ paste(ncbi_taxid, infraorder, sep = '_'),
  #       !is.na(suborder) & suborder != "" ~ paste(ncbi_taxid, suborder, sep = '_'),
  #       TRUE ~ order  # If both family and suborder are blank
  #     )) %>%
  #     mutate(fasta_id = gsub("_+$", "", fasta_id)) %>%
  #     unique() %>%
  #     select(ncbi_taxid, rec_id, genbank_accession, bold_id, bold_bin, lab_id, order, suborder, infraorder,
  #           superfamily,  family, subfamily, tribe, genus, subgenus, species, subspecies, species_id, fasta_id)
  # }
# }

# Write tree
if (is.null (opt$nexus)) {
  write.tree(tree, file=opt$output)
} else {
  write.nexus(tree, file=opt$output)
}
