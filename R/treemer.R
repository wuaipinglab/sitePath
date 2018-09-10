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
#' @importFrom ape nodepath getMRCA
#' @return path represent by tip names
#' @export
getSitePath <- function(tree, align, similarity) {
  align <- checkMatched(tree, align)
  # nodepath after trimming
  trimmedPaths <- unique(trimTree(nodepath(tree), align, similarity, FALSE))
  # get the bifurcated nodes and their path to the root in a trimmed tree
  # those paths are the so-called sitePaths (isolated)
  paths <- lapply(trimmedPaths, function(p) p[1:(length(p) - 1)])
  paths <- unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
  paths <- lapply(paths, function(p) {class(p) <- "isolated"; return(p)})
  # fuse sitePaths and add them to result
  paths <- c(paths, pathBeforeDivergence(paths))
  # children nodes of nodes on sitePaths
  root <- getMRCA(tree, tree$tip.label)
  attr(paths, "nodeTips") <- lapply(terminalNode(trimmedPaths), function(node) {
    childrenTips <- ChildrenTips(tree, node = node)
    res <- align[childrenTips]
    names(res) <- tree$tip.label[childrenTips]
    return(res)
  })
  attr(paths, "class") <- "sitePath"
  return(paths)
}

#' @export
print.sitePath <- function(sitePath) {
  cat(length(which(sapply(sitePath, class) == "isolated")), "isolated paths\n")
  cat(length(which(sapply(sitePath, class) == "fused")), "fused paths\n")
}

#' @export
print.isolated <- function(isolated) {
  class(isolated) <- NULL
  print(isolated)
}

#' @export
print.fused <- function(fused) {
  class(fused) <- NULL
  print(fused)
}


ChildrenTips <- function(tree, node) {
  maxTip <- length(tree$tip.label)
  if (all(node <= maxTip)) return(node)
  children <- integer(0)
  getChildren <- function(edges, parent) {
    i <- which(edges[, 1] %in% parent)
    if (length(i) == 0L) {
      return(children)
    } else {
      parent <- edges[i, 2]
      children <<- c(children, parent[which(parent <= maxTip)])
      return(getChildren(edges, parent))
    }
  }
  return(getChildren(tree$edge, node))
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

