#!/usr/bin/Rscript

# Example: Rscript tree_tips.R --nexus < tree.nex > tree_tips.txt
# TO get IDs and full names, use --csv

suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
is_nexus <- "--nexus" %in% args
csv_out <- "--csv" %in% args

input <- file('stdin', 'r')

if (!is_nexus) {
  tree <- read.tree(input)
} else {
  tree <- read.nexus(input)
}

if (csv_out) {
  

  
} else{
  cat(tree$tip.label, sep = "\n")
}
