library(dplyr)
library(tidyr)
library(getopt)

spec = matrix(c(
  'directory', 'd', 1, "character", "Directory with taxonomicindices.R output CSVs",
  'output'   , 'o', 1, "character", "Output CSV"
), byrow=T, ncol=5)
opt = getopt(spec)

# Directory containing CSV files
setwd(opt$directory)
directory <- getwd()
# directory <- sub("/$", "", directory)
getwd()

# List of CSV files in the directory
file_paths <- list.files(directory, pattern = "\\.raxml.csv$", full.names = TRUE)

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

write.csv(combined_df, opt$output, row.names = FALSE)
