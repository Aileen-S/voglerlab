# Filter output from BOLD v5 API

# Load packages
library(getopt)
library(dplyr)

spec <- matrix(c(
  'input',      'i', 1, 'character', 'Input BOLD metadata tsv to be filtered',
  'genbank',    'g', 2, 'logical',   'Remove sequences also on GenBank',
  'barcode',    'c', 2, 'logical',   'Save only barcodes, delete other genes',
  'longest',    'l', 2, 'character', 'Keep only longest sequence per BIN or per species (choices are "bin" or "species")'
), byrow = T, ncol = 5)

opt <- getopt(spec)


# Load data
load_data <- function(input) {
  meta <- read.csv(input, header = TRUE, sep = "\t", quote = "")  
  meta[meta == 'None'] <- NA
  meta[meta == ''] <- NA  
  # Remove records without sequences
  meta <- meta %>% filter(!is.na(marker_code) & !is.na(nuc))
  return(meta)
}


# Filter out GenBank sequences
filter_genbank <- function(meta) {
  meta <- meta %>% filter(inst != 'Mined from GenBank, NCBI')
  cat(nrow(meta), 'records remaining after removing those from GenBank\n')
  return(meta)
}


# Keep only COI 5' sequences
filter_barcodes <- function(meta) {
  meta <- meta %>% filter(markercode=="COI-5P")
  cat('Keeping only COI-5P sequences:', nrow(meta), 'found\n')
  return(meta)
}


filter_longest <- function(meta, longest) {
  # Get sequence lengths
  meta$nucleotides <- gsub("-","",as.character(meta$nucleotides))
  meta$sequence_length <- nchar(meta$nucleotides)

  # Keep longest sequence for each bin, for each gene
  if ( longest == 'bin') {
    no_bin <- meta %>% filter(is.na(bin_uri)) 
    unique_bin <- meta %>% 
      arrange(bin_uri, markercode, desc(sequence_length)) %>% 
      distinct(bin_uri, .keep_all = TRUE) 
    meta <- bind_rows(no_bin, unique_bin)

    cat(paste('Keeping one sequences for each unique BIN:', nrow(unique_bin), 'found\n'))
    cat(paste(nrow(no_bin), 'records without a BIN\n'))
    cat(paste('Keeping', nrow(meta), 'records total\n\n'))
  }  
  # Keep longest sequence for each species, for each gene
  else if ( longest == 'species') {
    no_species <- meta %>% filter(is.na(species)) 
    unique_species <- meta %>% 
      arrange(species, markercode, desc(sequence_length)) %>% 
      distinct(species, .keep_all = TRUE) 
    meta <- bind_rows(no_species, unique_species)

    cat(paste('Keeping one sequences for each unique species name:', nrow(unique_species), 'found\n'))
    cat(paste(nrow(no_species), 'records without a species name\n'))
    cat(paste('Keeping', nrow(meta), 'records total\n\n'))
  }
  return(meta)
}


split_coordinates <- function(meta) {
  meta$coord <- gsub("\\[|\\]", "", meta$coord)
  coords <- strsplit(meta$coord[!is.na(meta$coord)], ", ")
  meta$latitude[!is.na(meta$coord)] <- as.numeric(sapply(coords, function(x) x[1]))
  meta$longitude[!is.na(meta$coord)] <- as.numeric(sapply(coords, function(x) x[2]))
  return(meta)
}


# Write CSV metadata file
write_csv <- function(meta) {
  empty <- c('subgenus', 'subtribe',	'tribe', 'superfamily',	'infraorder',	'suborder')
  meta[ , empty] <- ''
  csv <- meta %>% 
    select(	processid, bin_uri,	suborder,	infraorder, superfamily, family, subfamily, 
          tribe, subtribe, genus, subgenus, species, subspecies, country.ocean,	latitude,	longitude) %>%
    distinct()
  new_names <- c("bold_id",	"bold_bin",	"suborder", "infraorder",	"superfamily", "family", "subfamily",	"tribe", 
                "subtribe",	"genus",	"subgenus",	"species", "subspecies", "country",	"latitude",	"longitude")
  names(csv) <- new_names

  # Write metadata to CSV
  write.csv(csv, 'bold_metadata.csv', row.names = FALSE)
  cat('Metadata saved to bold_metadata.csv\n\n')  
}


# Write fasta file for each gene
write_fasta <- function(meta){
  # Get list of gene names
  genes <-c(unique(meta$marker_code))
  cat('Genes found: ', genes, '\n')

  # Loop through genes and write to file
  for (gene in genes) {
    df <- meta %>% filter(marker_code==gene)
    file_out <- file(paste(gene, ".fasta", sep = ''), "w")
    for (i in 1:nrow(df)) {
      cat(">", df$processid[i], "\n", file = file_out, sep = '')
      cat(df$nuc[i], "\n", file = file_out, sep = '')
    }
    close(file_out)
    cat(length(df$nuc), ' sequences writen to ', gene, '.fasta\n', sep = '')
  }  
}


main <- function(opt) {
  meta <- load_data(opt$input)
  if (!is.null(opt$genbank)) {
    meta <- filter_genbank(meta)
  }
  if (!is.null(opt$barcode)) {
    meta <- filter_barcodes(meta)
  }
  if (!is.null(opt$longest)) {
    if (opt$longest %in% c("bin", "species")) {
      meta <- filter_longest(meta, opt$longest)
    } else {
      stop("Invalid choice for --longest. Must be 'bin' or 'species'.", call. = FALSE)
    }
  }
  meta <- split_coordinates(meta)
  write_csv(meta)
  write_fasta(meta)
}

main(opt)
