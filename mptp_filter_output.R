library(seqinr)
library(dplyr)
library(getopt)

# Filter fasta based on MPTP output - keep longest sequence for each delimited species
# Option to add CSV file with species names/IDs, to keep longest for each even if removed by PTP result
spec <- matrix(c(
  'input',    'i', 1, 'character', 'Input fasta (without taxonomy)',
  'output',   'o', 1, 'character', 'Output fasta',
  'mptp',     'm', 2, 'character', 'mPTP output .txt file',
  'csv',      'c', 2, 'character', 'Taxonomy CSV',
  'out_csv',  'd', 2, 'character', 'Output CSV with selected taxa',
  'list',     'l', 2, 'character', 'Output list of selected taxa',
  'tips',     't', 2, 'character', 'Column name in CSV with fasta IDs',
  'filter',   'f', 2, 'character', 'Column name with species names or ID to filter by',
  'strip',    's', 2, 'logical',   'Strip taxonomy from IDs in mptp file'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)


if (is.null(opt$mptp) & is.null(opt$csv)) {
  stop('Must specify either mPTP output or CSV file with species names/IDs')
}

# opt <- data.frame(input = 'COX1b.fasta',
#                   output = 'COX1b_filtered.fasta',
#                   mptp = 'COX1b.mptp.txt',
#                   strip = TRUE)

# Read fasta
aln <- read.fasta(opt$input)
paste(length(aln), 'sequences', 'in', opt$input)
lengths <- data.frame(rec_id = character(), length = integer())
for (id in names(aln)) {
  seq_nogap <- aln[[id]][aln[[id]] != "-"]
  len <- getLength(seq_nogap)
  lengths <- add_row(lengths, rec_id = id, length = len)
}

# Process mPTP file
if (!is.null(opt$mptp)) {
  mptp_file <- as.list(readLines(opt$mptp))
  mptp_df <- data.frame('ptp_species' = character(), 'rec_id' = character())
  start = FALSE
  for (line in mptp_file) {
    if (grepl('Species', line)) {
      sp_no <- as.numeric(gsub('\\D', '', line))
      spec <- sprintf("%0004d", sp_no)
      spec <- paste0('PTP', spec, collapse = '')
      start = TRUE }
    else if (start == TRUE & line != '') {
      if (!is.null(opt$strip)) {
        id_spl <- (strsplit(line, '_'))
        id_spl <- lapply(id_spl, grep, pattern = "[a-z]|OUTGROUP", invert = TRUE, value = TRUE)
        id_spl <- lapply(id_spl, function(x) x[x!=''])
        rec_id <- paste(id_spl[[1]], collapse = '_')      
      } else
        rec_id <- line
      mptp_df <- add_row(mptp_df, ptp_species = spec, rec_id = rec_id)
    }
  }
  print(paste(nrow(mptp_df), 'records and', as.character(sp_no), 'PTP groups in', opt$mptp))
  add <- merge(mptp_df, lengths, by = 'rec_id')

  # Find longest for each mPTP group
  mptp_select <- add %>%
    arrange(ptp_species, desc(length)) %>% 
    distinct(ptp_species, .keep_all = TRUE) %>% 
    ungroup()
  add <- add %>%
    mutate(selected = ifelse(rec_id %in% mptp_select$rec_id, 'PTP', NA))
} else {
  add <- lengths
  add$selected <- NA
}

# Add taxonomy filter
if (!is.null(opt$csv)) {
  meta <- read.csv(opt$csv)
  meta[meta == ''] <- NA
  meta <- meta %>% 
    select(opt$tips, opt$filter, source, binomial, order, suborder, infraorder, superfamily, family, subfamily, tribe, subtribe, genus, subgenus, species, subspecies, country, region, latitude, longitude) %>%
    rename(rec_id = opt$tips, filter = opt$filter)
  add <- merge(add, meta, by = 'rec_id', all.x = TRUE)

# Find longest for each id
  filter_select <- add %>%
    filter(!is.na(filter)) %>%
    arrange(filter, desc(length)) %>% 
    distinct(filter, .keep_all = TRUE) %>% 
    ungroup()
  print(paste(nrow(filter_select), 'unique values in', opt$filter, 'column of', opt$csv))
  if (!is.null(opt$mptp)) {
    filter_select <- filter_select %>%
      filter(!filter %in% mptp_select$filter)    
  }
  add <- add %>%
    mutate(selected = ifelse(rec_id %in% filter_select$rec_id & is.na(selected), 'taxonomy', selected))
}



add[grep("SRAB00086", add$rec_id), ]
lengths[grep("SRAB00086", lengths$rec_id), ]
mptp_df[grep("SRAB00086", mptp_df$rec_id), ]
meta[grep("SRAB00086", meta$rec_id), ]


# Write CSV
if (!is.null(opt$out_csv)) {
  add[is.na(add)] <- ''
  write.csv(add, opt$out_csv, row.names = FALSE)  
  print(paste('Output CSV written to', opt$out_csv))
}

# Write selected ID list
selected <- add %>% filter(selected != '')

if (!is.null(opt$list)) {
  writeLines(selected$rec_id, opt$list)
  print(paste('Output list of selected IDS written to', opt$list))
}

# Write fasta
output_fasta <- aln[c(which(names(aln) %in% selected$rec_id))]
write.fasta(output_fasta, names = names(output_fasta), file.out = opt$output)
paste('Output fasta with', length(output_fasta) ,'sequences written to', opt$output)
