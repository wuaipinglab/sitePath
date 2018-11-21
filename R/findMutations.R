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

#' @rdname Pre-assessment
#' @title Site with SNP
#' @description 
#' Test whether the frequency of amino acids in each site is enough to be an SNP
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @param minSNP minimum number of amino acid to be a SNP
#' @return a list of qualified site
#' @export
SNPsite <- function(tree, align, minSNP = NULL, reference = NULL, gapChar = '-') {
  if (is.null(minSNP)) {
    minSNP <- length(tree$tip.label) / 10
  }
  alignedSeq <- toupper(align$seq)
  if (is.null(reference)) {
    reference <- 1:nchar(alignedSeq[1])
  } else {
    reference <- getReference(alignedSeq[which(tree$tip.label == reference)], gapChar)
  }
  seqLen <- unique(nchar(alignedSeq))
  if (length(seqLen) != 1) stop("Sequence length not equal")
  qualified <- integer(0)
  alignedSeq <- strsplit(alignedSeq, "")
  for (i in  1:length(reference)) {
    SNP <- table(sapply(alignedSeq, "[[", reference[i]))
    if (sum(SNP > minSNP) >= 2) {
      qualified <- c(qualified, i)
    }
  }
  return(qualified)
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
#' @param reference name of reference for site numbering. The default uses the intrinsic alignment numbering
#' @param gapChar the character to indicate gap.
#' @param minSizeBefore minimum number of seuqence involved before mutation
#' @param minSizeAfter minimum number of seuqence involved after mutation
#' @param toleranceBefore maximum number of allowed sequence with a different amino acid before mutation
#' @param toleranceAfter maximum number of allowed sequence with a different amino acid after mutation
#' @export
findFixed.sitePath <- function(
  paths, reference = NULL, gapChar = '-',
  minSizeBefore = 10, minSizeAfter = 10,
  toleranceBefore = 2, toleranceAfter = 2
) {
  tree <- attr(paths, "tree")
  align <- attr(paths, "align")
  if (minSizeBefore <= 0 || minSizeBefore >= length(tree$tip.label)) {
    stop("minSizeBefore can't be negative, zero or greater than total number of sequence")
  } else if (minSizeAfter <= 0 || minSizeAfter >= length(tree$tip.label)) {
    stop("minSizeAfter can't be negative, zero or greater than total number of sequence")
  }
  if (is.null(reference)) {reference <- 1:nchar(align[1])} else {
    reference <- getReference(align[which(tree$tip.label == reference)], gapChar)
  }
  divNodes <- unique(divergentNode(paths))
  mutations <- list()
  for (minLen in 2:max(lengths(paths))) { # literate all sitePath at the same time
    for (path in unique(ancestralPaths(paths, minLen))) {
      afterTips <- ChildrenTips(tree, tail(path, 1))
      if (length(afterTips) < minSizeAfter) {
        next
      }
      pathBefore <- path[1:(length(path) - 1)]
      excludedTips <- sapply(pathBefore[which(pathBefore %in% divNodes)], function(node) {
        children <- tree$edge[which(tree$edge[, 1] == node), 2]
        children <- children[which(children > length(tree$tip.label) & !children %in% path)]
        return(ChildrenTips(tree, children))
      })
      beforeTips <- which(!1:length(tree$tip.label) %in% c(afterTips, unlist(excludedTips)))
      if (length(beforeTips) < minSizeBefore && length(excludedTips) == 0) {
        next
      }
      after <- strsplit(align[afterTips], "")
      before <- strsplit(align[beforeTips], "")
      for (i in 1:length(reference)) {
        bsum <- table(sapply(before, "[[", reference[i]))
        asum <- table(sapply(after, "[[", reference[i]))
        b <- names(bsum[1])
        a <- names(asum[1])
        if (sum(bsum[-1]) <= toleranceBefore && sum(asum[-1]) <= toleranceAfter && b != a) {
          mut <- paste(b, i, a, sep = "")
          from <- tree$tip.label[beforeTips]
          to <- tree$tip.label[afterTips]
          if (!mut %in% names(mutations) || length(to) > length(mutations[[mut]]$to)) {
            mutations[[mut]] <- list(from = from, to = to)
          }
        }
      }
    }
  }
  return(mutations)
}

#' @export
findFixed <- function(x, ...) UseMethod("findFixed")
