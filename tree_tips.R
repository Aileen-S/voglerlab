#!/usr/bin/Rscript

# Example: Rscript tree_tips.R --nexus < tree.nex > tree_tips.txt

library(ape, quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
is_nexus <- "--nexus" %in% args

input <- file('stdin', 'r')

if (!is_nexus) {
  tree <- read.tree(input)
} else {
  tree <- read.nexus(input)
}

cat(tree$tip.label, sep = "\n")