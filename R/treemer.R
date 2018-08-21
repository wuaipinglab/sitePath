#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

#' @title Group Tips
#' @description
#' Tree tips are grouped by their sequence similarity and members in a group
#' are constrained to share a same ancestral node.
#' Similarity between two tips is derived from their multiple sequence alignment.
#' The site doesn't count into total length if both are gap.
#' So similarity is calculated as number of matched divided by revised total length
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @param similarity similarity threshold for tree trimming
#' @param tipnames return tip as names or node number
#' @importFrom ape nodepath
#' @return grouping of tips
#' @export
groupTips <- function(tree, align, similarity, tipnames = TRUE) {
  grouping <- trimTree(
    nodepath(tree), checkMatched(tree, align),
    similarity, TRUE
  )
  if (tipnames) {
    return(lapply(grouping, function(g) {
      tree$tip.label[g]
    }))
  } else {
    return(grouping)
  }
}

#' @title Get sitePath
#' @description
#' Under development
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @param similarity similarity threshold for tree trimming
#' @importFrom ape nodepath
#' @return path as node numbers
#' @export
getSitePath <- function(tree, align, similarity) {
  # Get the nodepath after trimming
  paths <- unique(trimTree(
    nodepath(tree), checkMatched(tree, align),
    similarity, FALSE
  ))
  # Get the bifurcated nodes and their path to the root in a trimmed tree
  # Those paths are the so-called sitePaths
  paths <- lapply(paths, function(p) p[1:(length(p) - 1)])
  paths <- unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
  # Turning nodes are the end nodes of sitePaths plus the divergent nodes
  attr(paths, "turningNodes") <- c(
    sapply(paths, tail, 1),
    unique(divergentNode(paths))
  )
  attr(paths, "tree") <- tree
  attr(paths, "class") <- "sitePath"
  return(paths)
}

#' @export
print.sitePath <- function(sitePath) {
  cat(length(sitePath), "paths\n")
}

checkMatched <- function(tree, align, returnList = TRUE) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be of class phylo")
  } else if (!is(align, "alignment")) {
    stop("align must be of class alignment")
  }
  align$seq <- toupper(align$seq[match(tree$tip.label, align$nam)])
  if (any(is.na(align$seq))) {
    stop("tree tips and alignment names are not matched")
  } else if (length(unique(nchar(align$seq))) > 1) {
    stop("Sequence lengths are not the same in alignment")
  }
  if (returnList) {
    return(align$seq)
  } else {
    return(align)
  }
}
