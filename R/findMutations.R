#' @name SimilarityPlot
#' @rdname Pre-assessment
#' @title Pre-assessment
#' @description
#' Under development
#' @importFrom ape cophenetic.phylo
#' @export
sneakPeek <- function(tree, align, step = NULL) {
  distances <- cophenetic(tree)
  mostDist <- which(distances == max(distances), arr.ind = TRUE)
  mostDist <- compare(align$seq[[mostDist[1, 1]]], align$seq[[mostDist[1, 2]]])
  if (is.null(step)) {
    step <- round(1 - mostDist, 3) / 50
  }
  similarity <- numeric(0)
  pathNum <- integer(0)
  for (s in seq(mostDist, 1, step)) {
    similarity <- c(similarity, s)
    paths <- sitePath(tree, align, s)
    pathNum <- c(pathNum, length(paths))
  }
  plot(similarity, pathNum)
}

#' @export
sneakPeek2 <- function(tree, align, step = NULL) {
  distances <- cophenetic(tree)
  mostDist <- which(distances == max(distances), arr.ind = TRUE)
  mostDist <- compare(align$seq[[mostDist[1, 1]]], align$seq[[mostDist[1, 2]]])
  if (is.null(step)) {
    step <- round(1 - mostDist, 3) / 50
  }
  similarity <- numeric(0)
  mutationNum <- integer(0)
  for (s in seq(mostDist, 1, step)) {
    similarity <- c(similarity, s)
    paths <- sitePath(tree, align, s)
    mutations <- findFixed(paths)
    mutationNum <- c(mutationNum, sum(lengths(sapply(mutations, tail, 1))))
  }
  plot(similarity, mutationNum)
}

#' @name findMutations
#' @rdname findMutations
#' @title Find Mutations
#' @description
#' Under development
#' @param paths a \code{sitePath} object
#' @export
findFixed.sitePath <- function(paths, tolerance = 0, reference = NULL, gapChar = '-') {
  tree <- attr(paths, "tree")
  align <- attr(paths, "align")
  divNodes <- unique(divergentNode(paths))
  if (!is.null(reference)) {
    reference <- getReference(align[which(tree$tip.label == reference)], gapChar)
  } else {
    reference <- 1:nchar(align[1])
  }
  findMutation <- function(path) {
    afterTips <- ChildrenTips(tree, tail(path, 1))
    pathBefore <- path[1:(length(path) - 1)]
    excludedTips <- sapply(pathBefore[which(pathBefore %in% divNodes)], function(node) {
      children <- tree$edge[which(tree$edge[, 1] == node), 2]
      children <- children[which(children > length(tree$tip.label) & !children %in% path)]
      return(ChildrenTips(tree, children))
    })
    beforeTips <- which(!1:length(tree$tip.label) %in% c(afterTips, unlist(excludedTips)))
    after <- strsplit(align[afterTips], "")
    before <- strsplit(align[beforeTips], "")
    mutations <- character(0)
    for (i in 1:length(reference)) {
      bsum <- table(sapply(before, "[[", reference[i]))
      asum <- table(sapply(after, "[[", reference[i]))
      b <- names(bsum[1])
      a <- names(asum[1])
      if (sum(c(bsum[-1], asum[-1])) <= tolerance && b != a) {
        mutations <- c(mutations, paste(b, i, a, sep = ""))
      }
    }
    if (length(mutations) == 0) {
      return(character(0))
    } else {
      return(list(from  = tree$tip.label[beforeTips], to = tree$tip.label[afterTips], mutations = mutations))
    }
  }
  res <- list()
  for (n in 2:max(lengths(paths))) {
    ancestralPaths <- unique(ancestralPaths(paths, n))
    mutations <- lapply(ancestralPaths, findMutation)
    res <- c(res, mutations)
  }
  return(res[which(lengths(res) != 0)])
}

#' @export
findFixed <- function(x, ...) UseMethod("findFixed")

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
