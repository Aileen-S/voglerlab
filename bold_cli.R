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
  'barcode',    'b', 2, 'character', 'Save only barcodes, delete other genes',
  'metadata',   'm', 1, 'character', 'Read genbank metadata file to get TXIDs'
), byrow = T, ncol = 5)

opt <- getopt(spec)

names <- c('18S-3P' = '18S',
            '28S-D2' = '28S',
            'ArgK'   = 'AK',
            'COI-3P' = 'COX1b',
            'COI-5P' = 'COX1a',
            'COI-5PNMT1' = 'COX1a',
            'COII'   = 'COX2',
            'COIII' = 'COX3',
            'EF1-alpha' = 'EF1A',
            'ND5-0' = 'ND5',
            'Wnt1' = 'Wg')


##############################
# Get fastas and metadata
##############################

# Check either taxon or input file path have been 
#if ( !is.null(opt$taxon) & (!is.null(opt&csv) | (!is.null(opt$sequences) ) ) ) {
#  stop('Choose either taxon or csv&sequences, not both.')
#}

  # BOLD web search
if ( !is.null(opt$taxon) ) {
  out <- bold_seqspec(taxon=opt$taxon, sepfasta = TRUE)
  meta <- out[['data']]
  fasta <- out[['fasta']]
  print(paste(nrow(meta), 'records found for', opt$taxon))
  write.csv(meta, 'raw_metadata.csv', row.names = FALSE)
  print('Saved metadata to raw_metadata.csv')
  write.fasta(fasta, names(fasta), 'raw.fasta')
  print('Saved sequences to raw.fasta')
  
  #out <- bold_seqspec(taxon='Cybister', sepfasta = TRUE)
  
  # Search existing files  
  } else {
  meta <- read.csv(opt$csv)
  fasta <- read.fasta(opt$sequences)
  # Filter csv to remove records not present in fasta
  meta <- subset(meta, processid %in% names(fasta))
  } 

####################################
# Filter dataframe and write to file
####################################

# Add sp for empty species values
meta$species_name[which(meta$species_name == "")] <- paste(meta$genus_name[which(meta$species_name == "")], 'sp', sep = "_")

# Filter out GenBank sequences
if ( !is.null(opt$genbank) ) {
  f_meta <- meta %>% filter(genbank_accession=="")
  print(paste(nrow(f_meta), 'records remaining after those also on GenBank removed'))
} else {
    f_meta <- meta
  }
# Keep only COI sequences
if ( !is.null(opt$barcode)) {
  f_meta <- f_meta %>% filter(markercode=="COI-5P")
  print(paste(nrow(f_meta), 'records with COI-5P sequences'))
}

# Remove records without names genes or BINS
f_meta <- f_meta %>% filter(markercode!="")
f_meta <- f_meta %>% filter(bin_uri!="")


# Keep one record per bin
f_meta <- f_meta %>% distinct(bin_uri, markercode, .keep_all = TRUE)
print(paste(nrow(f_meta), 'unique BINs. Saved one sequence for each BIN.'))

# Get NCBI TXIDs
if ( !is.null(opt$metadata)) {
  ncbi <- read.csv(opt$metadata)
  ncbi <- ncbi %>% distinct(TXID, .keep_all = TRUE)
  ncbi <- ncbi %>% select(TXID, Species)
  names(ncbi)[names(ncbi)=="Species"] <- 'species_name'
  
  # Merge NCBI TXIDs with BOLD data
  f_meta <- merge(f_meta, ncbi, by = 'species_name', all.x = TRUE)
# 
} else {
  f_meta[ , 'TXID'] <- ''
}

# Add column for fasta IDs
f_meta$fasta_id <- ifelse(is.na(f_meta$TXID), 
                          paste(f_meta$processid, '/_', f_meta$family_name, '_', f_meta$subfamily_name, '_', f_meta$species_name, sep = ''), 
                          paste(f_meta$TXID, '/_', f_meta$processid, '/_', f_meta$family_name, '_', f_meta$subfamily_name, '_', f_meta$species_name, sep = ''))

f_meta$fasta_id <- gsub(" ", "_", f_meta$fasta_id)

#Edit dataframe to match genbank output
f_meta <- f_meta %>% select(fasta_id, processid, TXID, bin_uri, markercode, species_name, family_name, subfamily_name, genus_name, lat, lon)
empty <- c('Length', 'Suborder', 'Infraorder', 
           'Superfamily', 'Tribe', 'Description', 'Date Last Modified', 'Date Collected', 
           'Country', 'Region', 'latlon')
f_meta[ , empty] <- ''
f_meta <- f_meta %>% select(fasta_id, processid, TXID, bin_uri, species_name, markercode, Length, Suborder, Superfamily,
                            family_name, subfamily_name, genus_name, Description, 
                            'Date Last Modified', 'Date Collected', Country, Region, latlon, lat, lon)

# Write metadata to CSV
write.csv(f_meta, 'metadata.csv', row.names = FALSE)
print('Metadata saved to metadata.csv')

# Get list of gene names
genes <-c(unique(f_meta$markercode))
print('Genes found:')
print(genes)

# Replace with standard names
f_meta <- f_meta %>%  mutate(markercode = str_replace_all(markercode, names))
genes <-unique(f_meta$markercode)


#######################
# Filter and save fasta
#######################

for (gene in genes) {
  
  # Filter dataframe for gene
  df <- f_meta %>% filter(markercode==gene)

    # Save process IDs
  ids <- df$processid
  
  # Filter fasta
  fasta <- fasta[ids]
  
  # Add taxonomy to fasta
  # Update the identifiers in your list
  for (i in seq_along(fasta)) {
    old_identifier <- names(fasta)[i]
    new_identifier <- df$fasta_id[df$processid == old_identifier]
    # If a corresponding new identifier is found, update the name in the list
    if (!is.na(new_identifier) && length(new_identifier) == 1) {
      names(fasta)[i] <- new_identifier
    }
  }
  
  # Write the named list to a FASTA file
  #file = paste
  #write.fasta(fasta, names(fasta), paste(gene, '.fasta', sep = '' ))
  #print(paste('Sequences saved to', cat(gene, '.fasta', sep = '' )))

  file <- paste(gene, '.fasta', sep = '')
  write.fasta(fasta, names(fasta), file)
  cat(length(ids), "sequences saved to", file, "\n")
  
}
