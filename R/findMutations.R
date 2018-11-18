#' @rdname Pre-assessment
#' @title Similarity Matrix
#' @description 
#' To calculate similarity between aligned sequences
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @return a diagonal matrix of similarity between sequences
#' @export
getSimMatrix <- function(tree, align) {
  if (!inherits(tree, "phylo")) {
    stop("tree must be of class phylo")
  } else if (!is(align, "alignment")) {
    stop("align must be of class alignment")
  }
  sim <- similarityMatrix(checkMatched(tree, align))
  dimnames(sim) <- list(tree$tip.label, tree$tip.label)
  return(sim)
}

#' @name SimilarityPlot
#' @rdname Pre-assessment
#' @title Pre-assessment
#' @description
#' This function is intended to plot similarity as a threshold
#' against number of output sitePath. This plot is intended to give user
#' a feel of how many sitePaths they should expect from the similarity threshold
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @param maxPath maximum number of path to show in the plot
#' @export
sneakPeek <- function(tree, align, step = NULL, maxPath = NULL) {
  simMatrix <- getSimMatrix(tree, align)
  minSim <- min(simMatrix)
  if (is.null(step)) {
    step <- round(minSim - 1, 3) / 50
  }
  if (is.null(maxPath)) {
    maxPath <- length(tree$tip.label) / 20
  } else {
    if (maxPath <= 0) {
      stop("Maximum path number should be set greater than 0")
    }
  }
  similarity <- numeric(0)
  pathNum <- integer(0)
  for (s in seq(1, minSim, step)) {
    paths <- sitePath(tree, align, s, simMatrix)
    if (maxPath < length(paths)) {
      next
    } else if (length(paths) == 0) {
      break
    }
    similarity <- c(similarity, s)
    pathNum <- c(pathNum, length(paths))
  }
  plot(similarity, pathNum)
}

#' @name findMutations
#' @rdname findMutations
#' @title Find Mutations
#' @description
#' After finding the \code{\link{sitePath}} of a phylogenetic tree, we can use the result to find
#' those sites that have fixed on some sitePath but do not show fixation on the tree as a whole.
#' @param paths a \code{sitePath} object
#' @param minSize minimum number of seuqence involved before or after mutation
#' @param reference name of reference for site numbering. The default uses the intrinsic alignment numbering
#' @param gapChar the character to indicate gap.
#' @param tolerance maximum number of allowed sequence with a different amino acid
#' @export
findFixed.sitePath <- function(paths, minSize = 10, reference = NULL, gapChar = '-', tolerance = 0) {
  tree <- attr(paths, "tree")
  align <- attr(paths, "align")
  if (minSize <= 0 || minSize >= length(tree$tip.label)) {
    stop("minSize can't be negative, zero or greater than total number of sequence")
  }
  if (is.null(reference)) {
    reference <- 1:nchar(align[1])
  } else {
    reference <- getReference(align[which(tree$tip.label == reference)], gapChar)
  }
  divNodes <- unique(divergentNode(paths))
  res <- list()
  for (minLen in 2:max(lengths(paths))) { # literate all sitePath at the same time
    mutations <- lapply(unique(ancestralPaths(paths, minLen)), function(path) {
      afterTips <- ChildrenTips(tree, tail(path, 1))
      if (length(afterTips) < minSize) {
        return(NULL)
      }
      pathBefore <- path[1:(length(path) - 1)]
      excludedTips <- sapply(pathBefore[which(pathBefore %in% divNodes)], function(node) {
        children <- tree$edge[which(tree$edge[, 1] == node), 2]
        children <- children[which(children > length(tree$tip.label) & !children %in% path)]
        return(ChildrenTips(tree, children))
      })
      beforeTips <- which(!1:length(tree$tip.label) %in% c(afterTips, unlist(excludedTips)))
      if (length(beforeTips) < minSize) {
        return(NULL)
      }
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
        return(NULL)
      } else {
        return(list(from = tree$tip.label[beforeTips], to = tree$tip.label[afterTips], mutations = mutations))
      }
    })
    mutations <- mutations[which(lengths(mutations) != 0)]
    if (length(mutations) != 0) {res <- c(res, mutations)}
  }
  return(res)
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
