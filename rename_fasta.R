suppressMessages(library(dplyr))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(library(getopt))

spec <- matrix(c(
  'input',    'i', 1, 'character', 'Input fasta',
  'csv',      'c', 2, 'character', 'Metadata CSV. Specify custom labels with -l flag, or default assumes new names in first column, old names in second',
  'tips',     't', 2, 'character', 'Column name with original IDs (if not second column)',
  'label',    'l', 2, 'character', 'Custom labels: specify columns to use in labels: comma separated column names in order
                                    Default without this option is new names in first column, old names in second column',
  'keep_dups', 'k', 2, 'logical', 'Keep duplicate IDs in output under old IDs (default remove)',
  'output',   'o', 1, 'character', 'Output fasta',
  'renamed',  'r', 2, 'character',   'CSV output with old and new tips names',
  'drop_old', 'd', 2, 'logical',   'Drop original tip names (default keep old name at start of new name)'
  # 'prefix',   'p', 2, 'character', 'Add specified string to the start of all labels'
), byrow = TRUE, ncol = 5)
opt <- getopt(spec)


# Rename fasta
rename_fasta <- function(seqs, df, opt) {
  if (is.null(opt$tips)){
    df <- df %>% 
      filter(df[,2] %in% names(seqs)) %>%
      rename(old_id = names(.)[2]) %>%
      mutate(old_id = old_id)
  } else {
    df <- df %>% 
      filter(df[,opt$tips] %in% names(seqs)) %>%
      rename(old_id = opt$tips) %>%
      mutate(old_id = old_id)
  }

  # Get new IDs
  if (is.null(opt$label)) {
    df <- df %>%
      rename(fasta_id = names(df)[1])
  } else {
    column_names <- strsplit(opt$label, ",")[[1]]    
    if (is.null(opt$drop_old)) {
      df <- df %>%
        rowwise() %>%
        mutate(fasta_id = paste(c(old_id, c_across(all_of(column_names))), collapse = "_")) %>%
        ungroup()
    } else {
      df <- df %>%
        rowwise() %>%
        mutate(fasta_id = paste(c(c_across(all_of(column_names))), collapse = "_")) %>%
        ungroup()
    }
  }
  df <- df %>%
    mutate(fasta_id = gsub("_+$", "", fasta_id)) %>%
      select(old_id, fasta_id) %>%
      distinct()
  
  # Check for duplicates
  if (nrow(df) > length(unique(df$fasta_id))){
    lengths <- data.frame(rec_id = character(), seq_length = integer())
    for (id in names(seqs)) {
      seq_nogap <- seqs[[id]][!seqs[[id]] %in% c("-", "N")]
      len <- getLength(seq_nogap)
      lengths <- add_row(lengths, rec_id = id, seq_length = len)
    }
    df <- merge(df, lengths, by.x='old_id', by.y='rec_id', all=T)
    longest <- df %>%
      arrange(fasta_id, desc(seq_length)) %>% 
      distinct(fasta_id, .keep_all = TRUE) %>% 
      ungroup()

    if (is.null(opt$keep_dups) || opt$keep_dups == FALSE) {
      cat('Removing duplicate IDs')
      df <- df %>%
        mutate(removed = ifelse(old_id %in% longest$old_id, '', TRUE))
    } else {
      cat('Keeping duplicate IDs under old names')
      df <- df %>%
        mutate(fasta_id = ifelse(old_id %in% longest$old_id, fasta_id, old_id))
    }
  }
  df <- df %>% arrange(fasta_id)
  
  # Rename sequences
  lookup <- setNames(df$fasta_id, df$old_id)
  matched <- lookup[names(seqs)]
  renamed <- seqs
  new_labels <- ifelse(is.na(matched), names(renamed), matched)
  names(renamed) <- new_labels
  return(list(renamed_seqs = renamed, renamed_df = df))
}

# opt <- data.frame(
#   input = c('test.fasta'),
#   csv = c('~/scratch/adephaga/metadata/standardised_ptp.csv'),
#   tips = c('rec_id'),
#   label = c('ptp_id'),
#   output = c('test.out'),
#   renamed = c('test.csv'),
#   drop_old = c(T),
#   keep_dups = c(F)
# )

seqs <- read.fasta(opt$input)
tax <- read.csv(opt$csv)
output <- rename_fasta(seqs, tax, opt)
write.fasta(output$renamed_seqs, names(output$renamed_seqs), opt$output)
cat('Output written to', opt$output, '\n')

if (!is.null(opt$renamed)) {
  write.csv(output$renamed_df, opt$renamed, row.names=FALSE)
  cat('Output CSV written to', opt$renamed)
}
