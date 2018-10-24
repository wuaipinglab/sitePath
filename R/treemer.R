#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

#' @title Trim Tree
#' @description
#' Tree tips are grouped by their sequence similarity and members in a group
#' are constrained to share a same ancestral node.
#' Similarity between two tips is derived from their multiple sequence alignment.
#' The site doesn't count into total length if both are gap.
#' So similarity is calculated as number of matched divided by revised total length
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @param similarity similarity threshold for tree trimming
#' @param tipnames if return as tipnames
#' @importFrom ape nodepath
#' @return grouping of tips
#' @export
groupTips <- function(tree, align, similarity, tipnames = TRUE) {
  grouping <- trimTree(
    nodepath(tree), checkMatched(tree, align),
    similarity, TRUE
  )
  if (tipnames) {
    return(lapply(grouping, function(g) {tree$tip.label[g]}))
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
#' @return path represent by tip names
#' @export
sitePath <- function(tree, align, similarity) {
  align <- checkMatched(tree, align)
  # nodepath after trimming
  trimmedPaths <- unique(trimTree(nodepath(tree), align, similarity, FALSE))
  # get the bifurcated nodes and their path to the root in a trimmed tree
  # those paths are the so-called sitePaths (isolated)
  paths <- lapply(trimmedPaths, function(p) p[1:(length(p) - 1)])
  paths <- unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
  if (length(paths) == 0) {
    stop("Similarity threshold is too low resulting in no sitePath")
  }
  attr(paths, "tree") <- tree
  attr(paths, "align") <- align
  attr(paths, "class") <- "sitePath"
  return(paths)
}

#' @export
print.sitePath <- function(sitePath) {
  cat(length(sitePath), "paths\n")
}

checkMatched <- function(tree, align) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be of class phylo")
  } else if (!is(align, "alignment")) {
    stop("align must be of class alignment")
  }
  align <- toupper(align$seq[match(tree$tip.label, align$nam)])
  if (any(is.na(align))) {
    stop("tree tips and alignment names are not matched")
  } else if (length(unique(nchar(align))) > 1) {
    stop("Sequence lengths are not the same in alignment")
  }
  return(align)
}
