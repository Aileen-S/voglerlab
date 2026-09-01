library(seqinr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(stringr, quietly = TRUE)
library(getopt, quietly = TRUE)

spec <- matrix(c(
  'input',    'i', 1, 'character', 'Input fasta',
  'csv',      'c', 2, 'character', 'Metadata CSV. Specify custom labels with -l flag, or default assumes new names in first column, old names in second',
  'tips',     't', 2, 'character', 'Column name with original IDs (if not second column)',
  'label',    'l', 2, 'character', 'Custom labels: specify columns to use in labels: comma separated column names in order
                                    Default without this option is new names in first column, old names in second column',
  'output',   'o', 1, 'character', 'Output fasta',
  # 'renamed',  'r', 2, 'logical',   'CSV output with old and new tips names',
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
  
  lookup <- setNames(df$fasta_id, df$old_id)
  matched <- lookup[names(seqs)]
  renamed <- seqs
  new_labels <- ifelse(is.na(matched), names(renamed), matched)
  names(renamed) <- new_labels
  return(renamed)
}

# opt <- data.frame(
#   input = c('test.fasta'),
#   csv = c('~/scratch/adephaga/metadata/standardised_all_metadata.csv'),
#   tips = c('rec_id'),
#   label = c('suborder,superfamily,family,subfamily,tribe,genus,species,subspecies'),
#   output = c('12S_renamed.fasta')
# )

seqs <- read.fasta(opt$input)
df <- read.csv(opt$csv)
renamed <- rename_fasta(seqs, df, opt)
write.fasta(renamed, names(renamed), opt$output)
