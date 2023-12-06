library('ape')
library('tidyverse')
library('getopt')

# Get command line arguments
# column 1 = long flag
# column 2 = short flag
# column 3 = argument mask (0=no argument, 1=required argument, 2=optional argument)
# column 4 = argument type (logical, integer, double, complex, character)
# column 5 = optional, description/help message
spec <- matrix(c(
  'input',  'i', 1, 'character', 'Input tree file',
  'tips',   't', 1, 'character', 'File with tips to remove',
  'output', 'o', 1, 'character', 'Output tree file'
), byrow = T, ncol = 5)
opt <- getopt(spec)

# Open file with tips to be removed
tips <- strsplit(read_file(opt$tips), '\n')

# Save tips as vector
tips <- unlist(tips)

# Open tree
tree <- read.tree(opt$input)

# Remove tips
new <- drop.tip(tree, tips)

# Write new tree
write.tree(new, opt$output)

