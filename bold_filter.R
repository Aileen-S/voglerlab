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
  'tsv',        'i', 2, 'character', 'Input BOLD metadata tsv to be filtered',
  'genbank',    'g', 2, 'logical',   'Remove sequences also on GenBank',
  'barcode',    'b', 2, 'character', 'Save only barcodes, delete other genes',
  'raw',        'r', 2, 'logical',   'Save raw metadata, before any processing.',
  'bins',       'f', 2, 'logical',   'Keep only one sequence per BIN',
  'species',    's', 2, 'logical',   'Keep only one sequence per named species'
), byrow = T, ncol = 5)

opt <- getopt(spec)

# setwd('~/onedrive/treebuilding/Coleoptera/bold/')
# opt <- data.frame(
#   tsv   = 'raw_metadata.tsv',
#   genbank = TRUE,
#   barcode = TRUE,
#   bins = FALSE
# )

gene_names <- c('18S-3P' = '18S',
                '28S-D2' = '28S',
                'ArgK'   = 'AK',
                'COI-3P' = 'COX1b',
                'COI-5P' = 'COX1a',
                'COI-5PNMT1' = 'COX1a',
                'COII'   = 'COX2',
                'COXIII' = 'COX3',
                'EF1-alpha' = 'EF1A',
                'ND5-0' = 'ND5',
                'Wnt1' = 'Wg')


##############################
# Get fastas and metadata
##############################


  # BOLD web search
if ( !is.null(opt$taxon) ) {
  meta <- bold_seqspec(taxon=opt$taxon, cleanData = TRUE, fill=TRUE)
  print(paste(nrow(meta), 'records found for', opt$taxon))
  if ( !is.null(opt$raw) ) {
    write.table(meta, 'raw_metadata.tsv', row.names = FALSE, sep = '\t', quote = FALSE)
    print('Saved metadata to raw_metadata.tsv')
    }

  # Search existing files  
  } else {
  #meta <- read.table(opt$tsv, sep = '\t', header = TRUE)
  meta <- read.csv(opt$tsv, header = TRUE, sep = "\t", quote = "")  
  } 

  meta[meta == 'None'] <- NA

####################################
# Filter dataframe and write to file
###################################

# Filter out GenBank sequences

if (!is.null(opt$genbank)) {
  # Filter for non-empty and non-missing genbank_accession
  # f_meta <- meta %>% filter(is.na(genbank_accession) | genbank_accession == "")
  # print(paste(nrow(f_meta), 'records remaining after removing those with GenBank accession'))
  if ('institution_storing' %in% names(meta)) {
    f_meta <- meta[!grepl("GenBank", meta$institution_storing), ]
  } else if ('sequence_run_site' %in% names(meta)) {
    f_meta <- meta[!grepl("GenBank", meta$sequence_run_site), ]
  } else {
    f_meta <- meta
  }
  print(paste(nrow(f_meta), 'records remaining after removing those from GenBank'))
} else {
  f_meta <- meta
}

# Remove records without sequences or named genes
if ('marker_code' %in% names(f_meta)) {
  f_meta <- f_meta %>% mutate(markercode = marker_code)
}
f_meta <- f_meta %>% filter(markercode!="")
#f_meta <- f_meta %>% filter(bin_uri!="")
if ('nuc' %in% names(f_meta)) {
  f_meta <- f_meta %>% mutate(nucleotides = nuc)
}
f_meta <- f_meta %>% filter(nucleotides!="")

# Keep only COI sequences
if ( !is.null(opt$barcode)) {
  f_meta <- f_meta %>% filter(markercode=="COI-5P")
  print(paste(nrow(f_meta), 'records with COI-5P sequences'))
}

# Get sequence lengths
f_meta$sequence_length <- nchar(f_meta$nucleotides)

# Keep longest sequence for each bin, for each gene
if ( !is.null(opt$bins)) {
  # Filter rows where bin_uri is NA or empty
  na_or_empty <- f_meta %>% filter(is.na(bin_uri) | bin_uri == "") 
  # Filter rows where bin_uri is not NA or empty and keep one row per unique bin_uri 
  unique_bin_uri <- f_meta %>% 
    arrange(bin_uri, markercode, desc(sequence_length)) %>% 
    distinct(bin_uri, .keep_all = TRUE) 
  # Combine the rows 
  f_meta <- bind_rows(na_or_empty, unique_bin_uri)

  cat(paste(nrow(unique_bin_uri), 'unique BINs\nSaved one sequence for each BIN\n'))
  cat(paste(nrow(na_or_empty), 'records without a BIN\n'))
  cat(paste('Saved', nrow(f_meta), 'records total\n'))
}

# Keep longest sequence for each bin, for each gene
if ( !is.null(opt$species)) {
  # Filter rows where bin_uri is NA or empty
  na_or_empty <- f_meta %>% filter(is.na(species_name) | species_name == "") 
  # Filter rows where bin_uri is not NA or empty and keep one row per unique bin_uri 
  unique_species <- f_meta %>% 
    arrange(species_name, markercode, desc(sequence_length)) %>% 
    distinct(species_name, .keep_all = TRUE) 
  # Combine the rows 
  f_meta <- bind_rows(na_or_empty, unique_species)

  cat(paste(nrow(unique_species), 'unique speciess\nSaved one sequence for each species\n'))
  cat(paste(nrow(na_or_empty), 'records without a species name\n'))
  cat(paste('Saved', nrow(f_meta), 'records total\n'))
}


# Get NCBI TXIDs
if ( !is.null(opt$metadata)) {
  ncbi <- read.csv(opt$metadata)
  ncbi <- ncbi %>% distinct(ncbi_taxid, .keep_all = TRUE)
  ncbi <- ncbi %>% select(ncbi_taxid, species)
  names(ncbi)[names(ncbi)=="species"] <- 'species_name'

  # Merge NCBI TXIDs with BOLD data
  f_meta <- merge(f_meta, ncbi, by = 'species_name', all.x = TRUE)
# 
} else {
  f_meta[ , 'ncbi[!is.na(df$coords)]_taxid'] <- ''
}

#f_meta <- mutate(f_meta, TXID = replace(TXID, NA, ""))
#f_meta <- mutate(f_meta, markercode = replace(markercode, "COXIII", "COX3"))

# Remove gaps from sequences
f_meta$nucleotides <- gsub("-","",as.character(f_meta$nucleotides))

#Edit dataframe to match lab metadata
empty <- c('lab_id', 'subgenus', 'subtribe',	'tribe', 'superfamily',	'infraorder',	'suborder',	'genbank_accession')
f_meta[ , empty] <- ''
if ('family_name' %in% names(f_meta)) {
  csv <- f_meta %>% select(genbank_accession,	processid, bin_uri,	lab_id, suborder,	infraorder, superfamily, family_name, 
                         subfamily_name, tribe, subtribe, genus_name, subgenus, species_name, subspecies_name, country,	lat,	lon)
} else {
  f_meta$coord <- gsub("\\[|\\]", "", f_meta$coord)
  coords <- strsplit(f_meta$coord[!is.na(f_meta$coord)], ", ")
  f_meta$lat[!is.na(f_meta$coord)] <- as.numeric(sapply(coords, function(x) x[1]))
  f_meta$lon[!is.na(f_meta$coord)] <- as.numeric(sapply(coords, function(x) x[2]))
  colnames(f_meta)[colnames(f_meta) == "country.ocean"] <- "country"
  csv <- f_meta %>% select(genbank_accession,	processid, bin_uri,	lab_id, suborder,	infraorder, superfamily, family, 
  subfamily, tribe, subtribe, genus, subgenus, species, subspecies, country,	lat, lon)

}

new_names <- c("genbank_accession",	"bold_id",	"bold_bin",	"lab_id",	"suborder", "infraorder",	"superfamily",	
               "family", "subfamily",	"tribe", "subtribe",	"genus",	"subgenus",	"species", "subspecies", "country",	"latitude",	"longitude")
names(csv) <- new_names

csv <- csv %>% unique()


# Write metadata to CSV
write.csv(csv, 'metadata.csv', row.names = FALSE)
print('Metadata saved to metadata.csv')

# Get list of gene names
genes <-c(unique(f_meta$markercode))
print('Genes found:')
print(genes)

# Replace with standard names
f_meta <- f_meta %>%  mutate(markercode = str_replace_all(markercode, gene_names))
genes <-unique(f_meta$markercode)
print(genes)

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
    cat(">", df$processid[i], "\n", file = file_out, sep = '')
    # Write sequence data with line breaks
    cat(df$nucleotides[i], "\n", file = file_out, sep = '')
  }
  
  # Close the file connection
  close(file_out)
  
  #write.fasta(sequences = vec, names = names(vec), file.out = paste(gene, ".fasta", sep = ''))
  cat(length(df$nucleotides), ' sequences writen to ', gene, '.fasta\n', sep = '')
  
}
