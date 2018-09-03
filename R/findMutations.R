#' @rdname findMutations
#' @name findMutations
#' @title Find mutation
#' @description
#' Find the fixation site according to the evolutionary relationship.
#' @param paths a \code{sitePath} object
#' @param reference the name of the reference sequence for numbering
#' @param gapChar the character to indicate gap
#' @return sites and tip names
#' @export
findFixed.sitePath <- function(paths, reference = NULL, gapChar = "-") {
  tree <- attr(paths, "tree")
  align <- attr(paths, "align")
  paths <- unique(unlist(paths))
  if (!is.null(reference)) {
    reference <- getReference(align[which(tree$tip.label == reference)], gapChar)
  } else {
    reference <- 1:nchar(align[1])
  }
  fixed <- lapply(paths, function(node) {
    children <- ChildrenTips(tree$edge, node)
    ancestors <- which(!1:length(tree$tip.label) %in% children)
    ancestorSeq <- strsplit(align[ancestors], "")
    childrenSeq <- strsplit(align[children], "")
    mutations <- character(0)
    for (i in 1:length(reference)) {
      a <- unique(sapply(ancestorSeq, "[[", reference[i]))
      c <- unique(sapply(childrenSeq, "[[", reference[i]))
      if (length(a) == 1 && length(c) == 1 && a != c) {
        mutations <- c(mutations, paste(a, i, c, sep = ""))
      }
    }
    if (length(mutations) == 0) {
      return(character(0))
    } else {
      return(list(
        from = tree$tip.label[ancestors],
        to = tree$tip.label[children],
        mutations = mutations
      ))
    }
  })
  names(fixed) <- paths
  return(fixed[which(lengths(fixed) != 0)])
}

#' @export
findFixed <- function(x, ...) UseMethod("findFixed")

ChildrenTips <- function(edges, node) {
  maxTip <- length(tree$tip.label)
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
  return(getChildren(edges, node))
}
