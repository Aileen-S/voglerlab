# Merge rename_fasta.py CSV files to get combined CSV with old IDs and new IDs for each file
# Eg. merge_renamed_fasta_csv.R df1.csv df2.csv > output.csv

library(dplyr)

# Read in CSVs from rename_fasta.py
process_csv <- function(file) {
  data <- read.csv(file)
  data[is.na(data)] <- ''
  data <- data %>%
    filter(removed == "") %>%
    select(-removed)
  gene <- tools::file_path_sans_ext(basename(file))
  colnames(data)[colnames(data) == 'old_id'] <- gene
  message("Processed ", gene)
  return(data)
}

# Merge CSVs
merge_data_frames <- function(data_list) {
  Reduce(function(x, y) merge(x, y, by = "new_id", all = TRUE), data_list)
}


main <- function(file_paths) {
  cleaned_list <- lapply(file_paths, process_csv)
  combined <- merge_data_frames(cleaned_list)
  combined[is.na(combined)] <- ''
  message(nrow(combined), " taxa written to output file")
  write.csv(combined, stdout(), row.names = FALSE)
}

file_paths <- commandArgs(trailingOnly = TRUE)
main(file_paths)
