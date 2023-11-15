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
spec <- matrix(c(
  'taxon', 't', 1, 'character', 'BOLD search term',
  'fasta', 'f', 1, 'character', 'Output fasta file',
  'metadata', 'm', 1, 'character', 'Output metadata csv'
), byrow = T, ncol = 5)

opt <- getopt(spec)

##############################
# Download fastas and metadata
##############################

# Download metadata and fastas together

out <- bold_seqspec(taxon=opt$taxon, sepfasta = TRUE)
meta <- out[['data']]
fasta <- out[['fasta']]
print(paste(nrow(meta), ' records found'))
write.csv(meta, paste('raw', opt$metadata, sep = '_'), row.names = FALSE)
write.fasta(fasta, names(fasta), paste('raw', opt$fasta, sep = '_'))

####################################
# Filter dataframe and write to file
####################################

# Add sp for empty species values
meta$species_name[which(meta$species_name == "")] <- paste(meta$genus_name[which(meta$species_name == "")], 'sp', sep = "_")

# Add column for fasta IDs
meta$fasta_id <- paste(meta$processid, meta$family_name, meta$subfamily_name, meta$species_name, sep = '_')
meta$fasta_id <- gsub(" ", "_", meta$fasta_id)

# Filter out GenBank sequences
f_meta <- meta %>% drop_na(genbank_accession)
print(paste(nrow(f_meta), ' records remaining after those also on GenBank removed'))

f_meta <- f_meta %>% filter(markercode=="COI-5P"|markercode=="COI-3P")
print(paste(nrow(f_meta), ' records remaining with COI-5P sequences'))

f_meta <- f_meta %>% distinct(bin_uri, .keep_all = TRUE)
print(paste(nrow(f_meta), ' unique bins remaining'))

# Write metadata to CSV
write.csv(f_meta, opt$metadata, row.names = FALSE)
print(paste('Metadata saved to ', opt$meta))

# Get in same format as genbank metadata/mitogenome metadata?


#######################
# Filter and save fasta
#######################

# Save process IDs
ids <- filter$processid

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
print(paste('Sequences saved to ', opt$fasta))


