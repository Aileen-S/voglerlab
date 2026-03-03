# library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

# Check if file paths and output were provided
if (length(args) < 2) {
  stop("Please provide a list of CSV files and an output file as arguments")
}

# List of CSV files as positional arguments
file_paths <- args[1:length(args)-1]

# Output file
output <- args[length(args)]

combined_df <- data.frame()

for (file in file_paths) {
  data <- read.csv(file)
  print(file)
  Tree <- c(sub(".*\\/", "", file))
  Tips <- c(sum(data$ntips))
  TCI  <- c(mean(data$TCI, na.rm = TRUE))
  TRI  <- c(mean(data$TRI, na.rm = TRUE))
  df <- data.frame(Tree, Tips, TCI, TRI)
  print(df)
  combined_df <- rbind(combined_df, df)
}

write.csv(combined_df, output, row.names = FALSE)
