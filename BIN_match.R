library(bold)
library(seqinr)
library(dplyr)
library(getopt)

spec = matrix(c(
  'input',  'i', 1, "character", "Input fasta",
  'output', 'o', 1, "character", "Output CSV"
), byrow=T, ncol=5)
opt = getopt(spec)


seqs <- read.fasta(opt$input)

matches <- data.frame(matrix(ncol = 12, nrow = 0))
matches <- matches %>% mutate(across(everything(), as.character))
colnames(matches) <- c('db_id', 'BIN_match', 'ID', 'sequencedescription', 'database', 'citation', 
                       'taxonomicidentification', 'similarity', 'specimen_url',
                       'specimen_country', 'specimen_lat', 'specimen_lon')

for (name in names(seqs)) {
  seq <- paste(seqs[[name]], collapse = "")
  print(name)
  bin_info <- bold_identify(seq)
  print('completed bold_identify')
  match <- bin_info[[1]][1,]
  match_id <- ifelse(is.null(match$ID) || is.na(match$ID), "N/A", match$ID)
  cat(sprintf("Top result: %s\n", match_id))
  match['db_id'] = name
  match$similarity <- as.numeric(match$similarity)
  if (is.na(match$similarity)) {
    match['BIN_match'] = 'no'
  } else {
    if (match$similarity > 0.978) {
      match['BIN_match'] = 'yes'
    } else {
      match['BIN_match'] = 'no'
  }}
  match <- match %>% mutate(across(everything(), as.character))
  matches <- bind_rows(matches, match)
}

write.csv(matches, opt$output, row.names = FALSE)

