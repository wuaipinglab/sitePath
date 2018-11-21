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
#' @param simMatrix a diagonal matrix of similarity between sequences
#' @param tipnames if return as tipnames
#' @importFrom ape nodepath
#' @return grouping of tips
#' @export
groupTips <- function(tree, align, similarity, simMatrix = NULL, tipnames = TRUE) {
  simMatrix <- sortSimMatrix(tree, simMatrix)
  grouping <- trimTree(
    nodepath(tree), checkMatched(tree, align),
    simMatrix, similarity, TRUE
  )
  if (tipnames) {
    return(lapply(grouping, function(g) {tree$tip.label[g]}))
  } else {
    return(grouping)
  }
}

#' @title Get sitePath
#' @description
#' Find the sitePath of a phylogenetic tree providing the corresponding sequence alignment.
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @param similarity similarity threshold for tree trimming
#' @param simMatrix a diagonal matrix of similarity between sequences
#' @importFrom ape nodepath
#' @return path represent by tip names
#' @export
sitePath <- function(tree, align, similarity, simMatrix = NULL) {
  simMatrix <- sortSimMatrix(tree, simMatrix)
  align <- checkMatched(tree, align)
  # nodepath after trimming
  trimmedPaths <- unique(trimTree(nodepath(tree), align, simMatrix, similarity, FALSE))
  # get the bifurcated pre-terminal nodes and their path to the root in a trimmed tree
  # those paths are the so-called sitePaths (isolated)
  paths <- lapply(trimmedPaths, function(p) p[1:(length(p) - 1)])
  paths <- unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
  if (length(paths) == 0) {
    warning(paste0(similarity, " is too low of a cutoff resulting in no sitePath"))
  }
  paths <- lapply(paths, function(p) {
    sn <- p[length(p)]
    extended <- lapply(ChildrenTips(tree, sn), function(t) nodepath(tree, sn, t))
    el <- lengths(extended)
    ml <- max(el)
    longest <- extended[which(el == ml)]
    extended <- lapply(1:ml, function(i) unique(sapply(longest, "[[", i)))
    c(p, unlist(extended[which(lengths(extended) == 1)]))
  })
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
  if (!is(align, "alignment")) {
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

sortSimMatrix <- function(tree, simMatrix) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be of class phylo")
  }
  colMatch <- match(tree$tip.label, colnames(simMatrix))
  rowMatch <- match(tree$tip.label, rownames(simMatrix))
  if (is.null(simMatrix)) {
    return(matrix(
      NA,
      ncol = length(tree$tip.label),
      nrow = length(tree$tip.label)
    ))
  } else {
    return(simMatrix[rowMatch, colMatch])
  }
}

ChildrenTips <- function(tree, node) {
  maxTip <- length(tree$tip.label)
  children <- integer(0)
  getChildren <- function(edges, parent) {
    children <<- c(children, parent[which(parent <= maxTip)])
    i <- which(edges[, 1] %in% parent)
    if (length(i) == 0L) {
      return(children)
    } else {
      parent <- edges[i, 2]
      return(getChildren(edges, parent))
    }
  }
  return(getChildren(tree$edge, node))
}
