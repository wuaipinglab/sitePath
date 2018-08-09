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
#' @param treeAlignMatch a \code{\link{treeAlignMatch}} object
#' @param similarity similarity threshold for tree trimming
#' @return grouping of tips
#' @export
groupTips <- function(treeAlignMatch, similarity, tipnames = TRUE) {
  grouping <- trimTree(
    ape::nodepath(treeAlignMatch$tree),
    treeAlignMatch$align$seq,
    similarity, TRUE
  )
  if (tipnames) {
    return(lapply(grouping, function(g) {
      treeAlignMatch$tree$tip.label[g]
    }))
  } else {
    return(grouping)
  }
}

#' @title Get sitePath
#' @description 
#' Under development
#' @param treeAlignMatch a \code{\link{treeAlignMatch}} object
#' @param similarity similarity threshold for tree trimming
#' @return path represent by tip names
#' @export
getSitePath <- function(treeAlignMatch, similarity) {
  paths <- trimTree(
    ape::nodepath(treeAlignMatch$tree),
    treeAlignMatch$align$seq,
    similarity, FALSE
  )
  paths <- lapply(unique(paths), function(p) p[1:(length(p) - 1)])
  paths <- unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
  paths <- unique(divergentNode(paths))
  endNodes <- sapply(paths, tail, 1)
  return(paths)
}
