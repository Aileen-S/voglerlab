library(seqinr)
library(dplyr)
library(getopt)

# Filter fasta based on MPTP output - keep longest sequence for each delimited species
# Option to add CSV file with species names/IDs, to keep longest for each even if removed by PTP result
spec <- matrix(c(
  'input',    'i', 1, 'character', 'Input fasta (without taxonomy)',
  'output',   'o', 1, 'character', 'Output fasta',
  'mptp',     'm', 1, 'character', 'mPTP output .txt file',
  'csv',      'c', 2, 'character', 'Taxonomy CSV',
  'out_csv',  'm', 2, 'character', 'Output CSV with selected taxa',
  'list',     'l', 2, 'character', 'Output list of selected taxa',
  'tips',     't', 2, 'character', 'Column name in CSV with fasta IDs',
  'filter',   'f', 2, 'character', 'Column name with species names or ID to filter by',
  'strip',    's', 2, 'logical',   'Strip taxonomy from IDs in mptp file'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)

opt <- data.frame(
  input = c('~/scratch/adephaga/delimitation/01_raxml/12S.fasta'),
  output = c('~/scratch/adephaga/delimitation/12S_test.fasta'),
  mptp = c('~/scratch/adephaga/delimitation/02_mptp/12S.mptp.txt'),
  csv = c('~/scratch/adephaga/metadata/standardised_combined_metadata.csv'),
  out_csv = c('~/scratch/adephaga/delimitation/12S_ptp.csv'),
  list = c('~/scratch/adephaga/delimitation/12S_ptp.list'),
  tips = c('rec_id'),
  filter = c('binomial'),
  strip = c(TRUE)
)

# Process mPTP file
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
paste(nrow(mptp_df), 'records and', as.character(sp_no), 'PTP groups in', opt$mptp)

# Add sequence length
aln <- read.fasta(opt$input)
lengths <- data.frame(rec_id = character(), length = integer())
for (id in names(aln)) {
  seq_nogap <- aln[[id]][aln[[id]] != "-"]
  len <- getLength(seq_nogap)
  lengths <- add_row(lengths, rec_id = id, length = len)
}
mptp_df <- merge(mptp_df, lengths, by = 'rec_id')

# Find longest for each mPTP group
mptp_select <- mptp_df %>%
  arrange(ptp_species, desc(length)) %>% 
  distinct(ptp_species, .keep_all = TRUE) %>% 
  ungroup()
mptp_df <- mptp_df %>%
  mutate(selected = ifelse(rec_id %in% mptp_select$rec_id, 'PTP', NA))

# Add taxonomy filter
if (!is.null(opt$csv)) {
  meta <- read.csv(opt$csv)
  meta[meta == ''] <- NA
  meta <- meta %>% 
    select(opt$tips, opt$filter, source, binomial, order, suborder, infraorder, superfamily, family, subfamily, tribe, subtribe, genus, subgenus, species, subspecies, country, region, latitude, longitude) %>%
    rename(rec_id = opt$tips, filter = opt$filter)
  mptp_df <- merge(mptp_df, meta, by = 'rec_id', all.x = TRUE)

# Find longest for each id
  filter_select <- mptp_df %>%
    filter(!is.na(filter)) %>%
    arrange(filter, desc(length)) %>% 
    distinct(filter, .keep_all = TRUE) %>% 
    ungroup() %>%
    filter(!filter %in% mptp_select$filter)
  mptp_df <- mptp_df %>%
    mutate(selected = ifelse(rec_id %in% filter_select$rec_id & is.na(selected), 'taxonomy', selected))
}

# Write CSV
if (!is.null(opt$out_csv)) {
  mptp_df[is.na(mptp_df)] <- ''
  write.csv(mptp_df, opt$out_csv, row.names = FALSE)  
}

# Write selected ID list
if (!is.null(opt$list)) {
  selected <- mptp_df %>% filter(selected != '')
  writeLines(selected$rec_id, opt$list)
}

# Write fasta
output_fasta <- aln[c(which(names(aln) %in% all$rec_id))]
write.fasta(output_fasta, names = names(output), file.out = opt$output)
