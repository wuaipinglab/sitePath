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
#' @importFrom ape nodepath
#' @return grouping of tips
#' @export
groupTips <- function(treeAlignMatch, similarity, tipnames = TRUE) {
  grouping <- trimTree(
    nodepath(treeAlignMatch$tree),
    treeAlignMatch$align,
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
#' @importFrom ape nodepath
#' @importFrom phangorn Descendants
#' @return path represent by tip names
#' @export
getSitePath <- function(treeAlignMatch, similarity) {
  # trimTree gets the nodepath after trimming
  paths <- trimTree(
    nodepath(treeAlignMatch$tree),
    treeAlignMatch$align,
    similarity, FALSE
  )
  # nodepath after trimming is redundant
  paths <- unique(paths)
  # nodeTips are the tips which were trimmed
  nodeTips <- lapply(terminalNode(paths), function(node) {
    childrenTips <- unlist(Descendants(
      x = treeAlignMatch$tree, node = node, type = "tips"
    ))
    res <- treeAlignMatch$align[childrenTips]
    names(res) <- treeAlignMatch$tree$tip.label[childrenTips]
    return(res)
  })
  # Get the bifurcated nodes and their path to the root in a trimmed tree
  # Those paths are the so-called sitePaths
  paths <- lapply(paths, function(p) p[1:(length(p) - 1)])
  paths <- unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
  # fuse the paths
  fusedPaths <- unique(pathBeforeDivergence(paths))
  paths <- list(isolated = paths, fused = fusedPaths)
  attr(paths, "nodeTips") <- nodeTips
  attr(paths, "class") <- "sitePath"
  return(paths)
}

#' @export
print.sitePath <- function(sitePath) {
  cat(length(sitePath[["isolated"]]), "isolated paths\n")
  cat(length(sitePath[["fused"]]), "fused paths\n")
}
