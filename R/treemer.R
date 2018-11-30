#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

#' @rdname treemer
#' @name treemer
#' @title Topology-dependent tree trimming
#' @description
#' \code{groupTips} uses sequence similarity to group tree tips. Members in a group
#' are always constrained to share a same ancestral node.
#' Similarity between two tips is derived from their multiple sequence alignment.
#' The site will not be counted into total length if both are gap.
#' Similarity is calculated as number of matched divided by the corrected total length.
#' So far the detection of divergence is based on one simple rule: the similarity between
#' two most distant sequence. The two branches are decided to be divergent if the similarity
#' is lower than the threshold. (Other more statistical involved approaches 
#' such as Kolmogorov-Smirnov Tests among pair-wise distance could be introduced in the future)
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

#' @rdname treemer
#' @description
#' \code{sitePath} finds the sitePath of a phylogenetic tree providing the corresponding sequence alignment.
#' This is done by trimming the tree to the ancestor node of tips in each group and then
#' find the the bifurcated terminals of the trimmed tree. The \code{\link{nodepath}} between root node and
#' the bifurcated terminals is the sitePath. In order to extend the search of mutational site. The sitePath
#' will tag some of its trailing nodes. Here nodes up to the ancestor of the tips with
#' the longest \code{\link{nodepath}} are added.
#' @importFrom ape nodepath
#' @return path represent by node number
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
    stop("align is not class alignment")
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
    stop("tree is not class phylo")
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
