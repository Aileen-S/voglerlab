# Install packages if necessary
if (!require("bold")) {install.packages("bold")}
if (!require("tidyverse")) {install.packages("tidyverse")}
if (!require("getopt")) {install.packages("getopt")}
# Load packages
library('bold')
library('tidyverse')
library('getopt')
#library('seqinr')

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
  'genbank',    'g', 2, 'logical',   'Remove sequences also on GenBank',
  'barcode',    'b', 2, 'character', 'Save only barcodes, delete other genes',
  'metadata',   'm', 1, 'character', 'Read genbank metadata file to get TXIDs',
  'txid',       'x', 2, 'logical',   'Add NCBI taxon IDs to start of fasta ID where available'
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

#out <- bold_seqspec(taxon='Agabus', sepfasta = TRUE)


  # BOLD web search
if ( !is.null(opt$taxon) ) {
  meta <- bold_seqspec(taxon=opt$taxon)
  print(paste(nrow(meta), 'records found for', opt$taxon))
  write.csv(meta, 'raw_metadata.csv', row.names = FALSE)
  print('Saved metadata to raw_metadata.csv')

  # Search existing files  
  } else {
  meta <- read.csv(opt$csv)
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
f_meta <- f_meta %>% filter(nucleotides!="")

# Get sequence lengths
f_meta$sequence_length <- nchar(f_meta$nucleotides)

# Keep longest sequence for each bin, for each gene
#f_meta <- f_meta %>% distinct(bin_uri, markercode, .keep_all = TRUE)
f_meta <- f_meta %>%
  arrange(bin_uri, markercode, desc(sequence_length)) %>%
  distinct(bin_uri, markercode, .keep_all = TRUE)
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
} #else {
  #f_meta[ , 'TXID'] <- ''
#}

# Add column for fasta IDs
if ( !is.null(opt$txid)) {
  f_meta$fasta_id <- ifelse((is.na(f_meta$TXID) || f_meta$TXID == ''), 
                            paste(f_meta$processid, '/_', f_meta$family_name, '_', f_meta$subfamily_name, '_', f_meta$species_name, sep = ''), 
                            paste(f_meta$TXID, '/_', f_meta$processid, '/_', f_meta$family_name, '_', f_meta$subfamily_name, '_', f_meta$species_name, sep = ''))
} else {
  f_meta$fasta_id <-paste(f_meta$processid, '/_', f_meta$family_name, '_', f_meta$subfamily_name, '_', f_meta$species_name, sep = '')
}
f_meta$fasta_id <- gsub(" ", "_", f_meta$fasta_id)

# Remove gaps from sequences
f_meta$nucleotides <- gsub("-","",as.character(f_meta$nucleotides))

#Edit dataframe to match genbank output
#empty <- c('Suborder', 'Infraorder', 
#           'Superfamily', 'Tribe', 'Description', 'Date Last Modified', 'Date Collected', 
#           'Country', 'Region', 'latlon')
#f_meta[ , empty] <- ''
#csv <- f_meta %>% select(fasta_id, processid, TXID, bin_uri, species_name, markercode, sequence_length, Suborder, Superfamily,
#                            family_name, subfamily_name, genus_name, Description, 
#                            'Date Last Modified', 'Date Collected', Country, Region, latlon, lat, lon)
csv <- f_meta %>% select(fasta_id, processid, bin_uri, species_name, markercode, sequence_length, genbank_accession, phylum_taxID,	phylum_name,	
                         class_taxID,	class_name,	order_taxID,	order_name,	family_taxID,	family_name,	subfamily_taxID,	subfamily_name,	
                         genus_taxID,	genus_name,	species_taxID,	species_name,	subspecies_taxID,	subspecies_name,	habitat,	
                         lat,	lon,	coord_source,	coord_accuracy,	elev,	depth,	elev_accuracy,	depth_accuracy,	country,	province_state,	
                         region,	sector,	exactsite	)

# Write metadata to CSV
write.csv(csv, 'metadata.csv', row.names = FALSE)
print('Metadata saved to metadata.csv')

# Get list of gene names
genes <-c(unique(f_meta$markercode))
print('Genes found:')
print(genes)

# Replace with standard names
f_meta <- f_meta %>%  mutate(markercode = str_replace_all(markercode, names))
genes <-unique(f_meta$markercode)
#print(genes)

#######################
# Filter and save fasta
#######################

for (gene in genes) {
  
  # Filter dataframe for gene
  df <- f_meta %>% filter(markercode==gene)
  
  file_out <- file(paste(gene, ".fasta", sep = ''), "w")
  # Write to fasta
  for (i in 1:nrow(df)) {
    # Write header line
    cat(">", df$fasta_id[i], "\n", file = file_out, sep = '')
    # Write sequence data with line breaks
    cat(df$nucleotides[i], "\n", file = file_out, sep = '')
  }
  
  # Close the file connection
  close(file_out)
  
  #write.fasta(sequences = vec, names = names(vec), file.out = paste(gene, ".fasta", sep = ''))
  cat(length(df$nucleotides), ' sequences writen to ', gene, '.fasta\n', sep = '')
  
}
