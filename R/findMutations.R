#' @name findMutations
#' @rdname findMutations
#' @title Find Mutations
#' @description
#' Under development
#' @param paths a \code{sitePath} object
#' @importFrom ape nodepath getMRCA
#' @export
findFixed.sitePath <- function(paths, reference = NULL, gapChar = '-') {
  tree <- attr(paths, "tree")
  align <- attr(paths, "align")
  if (!is.null(reference)) {
    reference <- getReference(align[which(tree$tip.label == reference)], gapChar)
  } else {
    reference <- 1:nchar(align[1])
  }
  findMutation <- function(path, divNodes) {
    afterTips <- ChildrenTips(tree, tail(path, 1))
    pathBefore <- path[1:(length(path) - 1)]
    excludedTips <- sapply(pathBefore[which(pathBefore %in% divNodes)], function(node) {
      children <- tree$edge[which(tree$edge[, 1] == node), 2]
      children <- children[which(children > length(tree$tip.label) & !children %in% path)]
      return(ChildrenTips(tree, children))
    })
    beforeTips <- which(!1:length(tree$tip.label) %in% c(afterTips, excludedTips))
    after <- strsplit(align[afterTips], "")
    before <- strsplit(align[beforeTips], "")
    mutations <- character(0)
    for (i in 1:length(before[[1]])) {
      b <- unique(sapply(before, "[[", i))
      a <- unique(sapply(after, "[[", i))
      if (length(b) == 1 && length(a) == 1 && a != b) {
        mutations <- c(mutations, paste(b, i, a, sep = ""))
      }
    }
    if (length(mutations) == 0) {
      return(character(0))
    } else {
      return(list(
        from  = tree$tip.label[beforeTips],
        to = tree$tip.label[afterTips],
        mutations = mutations
      ))
    }
  }
  res <- list()
  for (n in 1:max(lengths(paths))) {
    ancestralPaths <- unique(ancestralPaths(paths, n))
    ancestralPaths <- ancestralPaths[which(lengths(ancestralPaths) > 1)]
    mutations <- lapply(ancestralPaths, findMutation, attr(paths, "divNodes"))
    res <- c(res, mutations)
  }
  return(res[which(lengths(res) != 0)])
}

#' @export
findFixed <- function(x, ...) UseMethod("findFixed")

ChildrenTips <- function(tree, node) {
  maxTip <- length(tree$tip.label)
  children <- node[which(node <= maxTip)]
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
