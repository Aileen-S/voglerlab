# Packages
library('bold')
library('tidyverse')
library('seqinr')
library('getopt')

# Get command line arguments
# column 1 = long flag
# column 2 = short flag
# column 3 = argument mask (0=no argument, 1=required argument, 2=optional argument)
# column 4 = argument type (logical, integer, double, complex, character)
# column 5 = optional, description/help message

#####################################
# Get COX1 sequences from BOLD. Choose either taxon to download from
# web or csv and sequences to search existing files.
# Filters one per BIN.

spec <- matrix(c(
  'taxon',      't', 2, 'character', 'BOLD search term',
  'csv',        'c', 2, 'character', 'Path to existing BOLD metadata to be filtered',
  'sequences',  's', 2, 'character', 'Path to existing fasta to be filtered',
  'genbank',    'g', 2, 'logical',   'Remove sequences also on GenBank',
  'fasta',      'f', 1, 'character', 'Output fasta file',
  'metadata',   'm', 1, 'character', 'Output metadata csv'
), byrow = T, ncol = 5)

opt <- getopt(spec)

##############################
# Get fastas and metadata
##############################

# Check either taxon or input file path have been 
#if ( !is.null(opt$taxon) & (!is.null(opt&csv) | (!is.null(opt$sequences) ) ) ) {
#  stop('Choose either taxon or csv&sequences, not both.')
#}

if ( !is.null(opt$taxon) ) {
  out <- bold_seqspec(taxon=opt$taxon, sepfasta = TRUE)
  meta <- out[['data']]
  fasta <- out[['fasta']]
  print(paste(nrow(meta), 'records found for', opt$taxon))
  write.csv(meta, paste('raw', opt$metadata, sep = '_'), row.names = FALSE)
  print(paste('Saved metadata to raw_', opt$metadata, sep = ''))
  write.fasta(fasta, names(fasta), paste('raw', opt$fasta, sep = '_'))
  print(paste('Saved sequences to raw_', opt$fasta, sep = ''))
  } else {
  meta <- read.csv(opt$csv)
  fasta <- read.fasta(opt$sequences)
  } 

####################################
# Filter dataframe and write to file
####################################

# Add sp for empty species values
meta$species_name[which(meta$species_name == "")] <- paste(meta$genus_name[which(meta$species_name == "")], 'sp', sep = "_")

# Add column for fasta IDs
meta$fasta_id <- paste(meta$processid, meta$family_name, meta$subfamily_name, meta$species_name, sep = '_')
meta$fasta_id <- gsub(" ", "_", meta$fasta_id)

# Filter out GenBank sequences
if ( !is.null(opt$genbank) ) {
  f_meta <- meta %>% drop_na(genbank_accession)
  print(paste(nrow(f_meta), 'records remaining after those also on GenBank removed'))
} else {
    f_meta <- meta
  }
# Keep only COI sequences
f_meta <- f_meta %>% filter(markercode=="COI-5P"|markercode=="COI-3P")
print(paste(nrow(f_meta), 'records with COI-5P sequences'))

# Keep one record per bin
f_meta <- f_meta %>% distinct(bin_uri, .keep_all = TRUE)
print(paste(nrow(f_meta), 'unique BINs. Saved one sequence for each BIN.'))

# Write metadata to CSV
write.csv(f_meta, opt$metadata, row.names = FALSE)
print(paste('Metadata saved to', opt$meta))

# Get in same format as genbank metadata/mitogenome metadata?


#######################
# Filter and save fasta
#######################

# Save process IDs
ids <- f_meta$processid

# Filter fasta
fasta <- fasta[ids]

# Add taxonomy to fasta
# Update the identifiers in your list
for (i in seq_along(fasta)) {
  old_identifier <- names(fasta)[i]
  new_identifier <- f_meta$fasta_id[f_meta$processid == old_identifier]
  # If a corresponding new identifier is found, update the name in the list
  if (!is.na(new_identifier) && length(new_identifier) == 1) {
    names(fasta)[i] <- new_identifier
  }
}

# Write the named list to a FASTA file
write.fasta(fasta, names(fasta), opt$fasta)
print(paste('Sequences saved to', opt$fasta))


