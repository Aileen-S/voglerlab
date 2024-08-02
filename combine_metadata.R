# Load necessary libraries
library(dplyr)
library(getopt)

# Define the command-line arguments
spec <- matrix(c(
  'bold',      'b', 1, 'character', 'BOLD metadata file',
  'genbank',   'g', 1, 'character', 'GenBank metadata file',
  'lab',       'l', 1, 'character', 'lab MMG ids file',
  'outgroup',  'o', 1, 'character', 'lab outgroup MMG ids file',
  'mmg',       'm', 1, 'character', 'lab MMG metadata file',
  'output',    'c', 1, 'character', 'output CSV file'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

# Print the new working directory
print(paste("Working directory set to:", getwd()))

# GenBank
if (!is.null(opt$genbank)) {
  print('Getting GenBank metadata')
  #gen <- read.csv('genbank/metadata.csv')  
  gen <- read.csv(opt$genbank)
  print(paste(nrow(gen), 'rows in', opt$genbank))
  gen <- gen %>% 
  select(ncbi_taxid, genbank_accession, bold_id, bold_bin, lab_id, suborder, infraorder, superfamily, 
                      family, subfamily, tribe, species) %>%
  mutate(rec_id = genbank_accession)
  # Make NCBI taxon ID to add to other files
  ncbi <- gen %>% select(ncbi_taxid, species) %>% filter(!grepl("sp$", species)) %>% distinct()
  names(ncbi)[names(ncbi) == 'ncbi_taxid'] <- 'ncbi'
}

# BOLD
if (!is.null(opt$bold)) {
  print('Getting BOLD metadata')
  #bold <- read.csv('bold/metadata.csv')
  bold <- read.csv(opt$bold)
  print(paste(nrow(bold), 'rows in', opt$bold))
  bold <- bold %>% 
  select(ncbi_taxid, genbank_accession, bold_id, bold_bin, lab_id, suborder, infraorder, superfamily, 
                      family, subfamily, tribe, species) %>%
  mutate(rec_id = bold_id)
}

lab_ids <- character()

# Get lab MMG ids
if (!is.null(opt$lab)) {
  print('Getting lab IDs')
  lab <- readLines(opt$lab)
  print(paste(length(lab), 'ids in', opt$lab))
  lab_ids <- c(lab_ids, lab)
}

# Get outgroup MMG ids
if (!is.null(opt$outgroup)) {
  print('Getting outgroup IDs')
  out <- readLines(opt$outgroup)
  print(paste(length(out), 'ids in', opt$outgroup))
  lab_ids <- c(lab_ids, out)
}

# Remove duplicates from the list
lab_ids <- unique(lab_ids)
# Get lab MMG metadata
if (!is.null(opt$mmg)) {
  print('Getting MMG metadata')
  mmg <- read.csv(opt$mmg)
  # Filter rows with ID list
  mmg <- mmg %>% filter(mt_id %in% lab_ids)
  print(paste(nrow(mmg), 'of', length(lab_ids), 'requested ids founds in', opt$mmg))
  # Fix column names
  empty <- c('bold_id',	'bold_bin')
  mmg[ , empty] <- ''
  names(mmg)[names(mmg) == 'mt_id'] <- 'lab_id'
  # Remove non-species level ncbi taxon IDs
  mmg <- mmg %>%
    mutate(ncbi_taxid = if_else(ncbi_id_rank != 'species', NA_integer_, ncbi_taxid)) %>% 
    select(ncbi_taxid, genbank_accession, bold_id, bold_bin, lab_id, suborder, infraorder, superfamily, 
           family, subfamily, tribe, species)  %>%
    mutate(rec_id = lab_id)
}

# # Check if each dataframe exists
# dfs <- list(
#   gen = if (exists("gen", envir = .GlobalEnv)) gen else NULL,
#   bold = if (exists("bold", envir = .GlobalEnv)) bold else NULL,
#   mmg = if (exists("lab", envir = .GlobalEnv)) mmg else NULL,
# )
# 
# # Remove NULL elements from the list
# dfs <- Filter(Negate(is.null), data_frames)

# Combine metadata
all <- bind_rows(gen, bold, mmg)
all <- all %>%
  distinct() %>%
  # Clean species names
  mutate(species = gsub(" sp[.,].*$", " sp", species))


# Merge NCBI TXIDs
add <- merge(all, ncbi, by = 'species', all.x = TRUE)

add <- add %>%
  mutate(ncbi_taxid = if_else(is.na(ncbi_taxid), ncbi, ncbi_taxid))
all <- add %>% select(ncbi_taxid, rec_id, genbank_accession, bold_id, bold_bin, lab_id, suborder, infraorder, superfamily, 
                                            family, subfamily, tribe, species)
all[is.na(all)] <- ''

# Print metadata to rename with NCBI taxids
#write.csv(all, 'supermatrix/rename_1.csv', row.names = FALSE)

# Get metadata to rename with taxonomy
all <- all %>%
  mutate(spec_id = if_else(ncbi_taxid == '', rec_id, ncbi_taxid)) %>%
  #rename(fasta_id_2 = ncbi_taxid, fasta_id_1 = rec_id) %>%
  mutate(species = replace(species, is.na(species) | species == "", "sp")) %>%
  mutate(species = gsub(" ", "_", species)) %>%
  mutate(taxonomy = paste(spec_id, family, subfamily, tribe, species, sep = '_')) %>%
  select(taxonomy, rec_id, ncbi_taxid, genbank_accession, bold_id, bold_bin, lab_id, suborder, infraorder, superfamily, 
         family, subfamily, tribe, species)
  
write.csv(all, 'supermatrix/rename.csv', row.names = FALSE)
paste('Combined metadata written to', opt$output)
