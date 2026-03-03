suppressMessages(library(ape))
suppressMessages(library(tidyverse))
suppressMessages(library(getopt))

# Get command line arguments
# column 1 = long flag
# column 2 = short flag
# column 3 = argument mask (0=no argument, 1=required argument, 2=optional argument)
# column 4 = argument type (logical, integer, double, complex, character)
# column 5 = optional, description/help message
spec <- matrix(c(
  'input',  'i', 1, 'character', 'Input tree file',
  'tips',   't', 1, 'character', 'File with tip names to search for: either in list or newick tree. Default is to remove these tips from tree',
  'c_tips', 'c', 1, 'character', 'Comma separated list of tip names. Default is to remove these tips from tree',
  'keep',   'k', 2, 'logical',   'Keep specified tips and remove the rest',
  'output', 'o', 1, 'character', 'Output tree file'
), byrow = T, ncol = 5)
opt <- getopt(spec)


read_tips_file <- function(tips_file, tips_list) {
  # Tips file
  if (!is.null(tips_file)) {
    # Try to read tips file as a tree file
    tips <- tryCatch(
      expr = {
        tree <- read.tree(tips_file)
        tips <- tree$tip.label
      },
        warning = function(e) {
          tips <- trimws(strsplit(readLines(tips_file), '\n'))
          tips <- unlist(tips)
      }
    )
  } else if (!is.null(tips_list)) {
    tips <- strsplit(tips_list, ',')
    tips <- unlist(tips)
  } else {
    stop('Please provide tip name(s) with either -t or -c')
  }
  cat(length(tips), 'IDs found in list of tips to filter\n')
  return(tips)
}


read_tree_file <- function(tree_file) {
  first_line <- readLines(tree_file, n = 1)
  
  # Read nexus
  if (grepl("#NEXUS", first_line, ignore.case = TRUE)) {
    nexus <- TRUE
    tree <- ape::read.nexus(tree_file)
    
    # Multiple trees
    if (inherits(tree, "multiPhylo")) {
      cat("Tree file contains", length(tree), "nexus trees\n")
      # Apply tip label correction to each tree in the list
      tree <- lapply(tree, function(t) {
        t$tip.label <- gsub("(^['\"]|['\"]$)", '', t$tip.label)
        return(t)
      })
      class(tree) <- "multiPhylo"
      
    # Single tree
    } else {
      cat("Tree has", length(tree$tip.label), "tips\n")
      tree$tip.label <- gsub("(^['\"]|['\"]$)", '', tree$tip.label)
    }

  } else {
    # Read newick
    nexus <- FALSE
    tree_lines <- readLines(tree_file)
    tree_lines <- tree_lines[nchar(tree_lines) > 0]
    
    # Single tree
    if (length(tree_lines) == 1) {
      tree <- ape::read.tree(text = tree_lines)
      cat("Tree has", length(tree$tip.label), "tips\n")
      tree$tip.label <- gsub("(^['\"]|['\"]$)", '', tree$tip.label)
      
    # Multiple trees
    } else {
      tree <- lapply(tree_lines, function(x) {
        t <- ape::read.tree(text = x)
        t$tip.label <- gsub("(^['\"]|['\"]$)", '', t$tip.label)
        return(t)
      })
      class(tree) <- "multiPhylo"
      cat("Tree file contains", length(tree), "trees from a Newick file.\n")
    }
  }
  
  return(list(tree = tree, nexus = nexus))
}


filter_tips <- function(tree, tips, keep) {
  common_tips <- intersect(tree$tip.label, tips)
  # Keep tips in tips file if -f flag used
  if ( !is.null(keep) ) {
    cat('Keeping tips present in input list', '\n')
    new <- keep.tip(tree, common_tips)
  # Otherwise remove tips in tips file
  } else {
    cat('Removing tips present in input list', '\n')
    new <- drop.tip(tree, common_tips)
  }  
  cat(length(new$tip.label), 'tips remaining\n')
  return(new)
}


write_tree <- function(tree, output, nexus) {
  # tree <- rename_tips(tree, df)
  if (nexus == FALSE) {
    write.tree(tree, file=output)
    cat("Output written to", output, '\n')
  } else {
    write.nexus(tree, file=output)
    cat("Output written to", output, '\n')
  }
}


main <- function(opt) {
  tips <- read_tips_file(opt$tips, opt$c_tips)
  # Check if input is one or multiple trees
  tree <- read_tree_file(opt$input)
  if (class(tree$tree) == "multiPhylo") {
    output_tree <- lapply(tree$tree, function(t) filter_tips(t, tips, opt$keep))
    class(output_tree) <- "multiPhylo"
  } else {
    output_tree <- filter_tips(tree$tree, tips, opt$keep)
  }  
  write_tree(output_tree, opt$output, tree$nexus)
}


main(opt)