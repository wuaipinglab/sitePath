#' @rdname Pre-assessment
#' @name Pre-assessment
#' @title Things can be done before the analysis
#' @description \code{similarityMatrix} calculates similarity between aligned sequences
#' The similarity matrix can be used in \code{\link{groupTips}} or \code{\link{sitePath}}
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @return \code{similarityMatrix} returns a diagonal matrix of similarity between sequences
#' @export
similarityMatrix <- function(tree, align) {
  if (!inherits(tree, "phylo")) {
    stop("tree is not class phylo")
  } else if (!is(align, "alignment")) {
    stop("align is not class alignment")
  }
  sim <- getSimilarityMatrix(checkMatched(tree, align))
  dimnames(sim) <- list(tree$tip.label, tree$tip.label)
  return(sim)
}

#' @rdname Pre-assessment
#' @description
#' \code{sneakPeek} is intended to plot similarity as a threshold
#' against number of output sitePath. This plot is intended to give user
#' a feel about how many sitePaths they should expect from the similarity threshold.
#' The number of sitePath should not be too many or too few. The result excludes
#' where the number of sitePath is greater than number of tips divided by 20 or
#' self-defined maxPath. The zero sitePath result will also be excluded
#' @param step the similarity window
#' @param maxPath maximum number of path to show in the plot
#' @param makePlot whether make a dot plot when return
#' @return \code{sneakPeek} return the similarity threhold against number of sitePath
#' @export
sneakPeek <- function(tree, align, step = NULL, maxPath = NULL, makePlot = TRUE) {
  simMatrix <- similarityMatrix(tree, align)
  minSim <- min(simMatrix)
  if (is.null(step)) {
    step <- round(minSim - 1, 3) / 50
  }
  if (is.null(maxPath)) {
    maxPath <- length(tree$tip.label) / 20
  } else {
    if (maxPath <= 0) {
      stop("Maximum path number is not greater than 0")
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
  if (makePlot) {plot(similarity, pathNum)}
  return(data.frame(similarity, pathNum))
}

#' @rdname findSites
#' @name findSites
#' @title Finding sites with variation
#' @description 
#' Single nucleotide polymorphism (SNP) in the whole package refers to variation of amino acid.
#' \code{findSNPsite} will try to find SNP in the multiple sequence alignment. A reference sequence
#' and gap character may be specified to number the site. This us irrelevant to the intended analysis
#' but might be helpful to evaluate the performance of \code{fixationSites}
#' @param tree a \code{phylo} object
#' @param align an \code{alignment} object
#' @param reference name of reference for site numbering. 
#' The name has to be one of the sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar the character to indicate gap.
#' The numbering will skip the gapChar if reference sequence if specified.
#' @param minSNP minimum number of amino acid variation to be a SNP
#' @return \code{findSNPsite} returns a list of qualified SNP site
#' @export
SNPsites <- function(tree, align, reference = NULL, gapChar = '-', minSNP = NULL) {
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

#' @rdname findSites
#' @name fixationSites
#' @description
#' After finding the \code{\link{sitePath}} of a phylogenetic tree, we use the result to find
#' those sites that show fixation on some if not all sitePath. Parallel evolution is relatively
#' common in RNA virus. There is chance that some site be fixed in one lineage but does not show
#' fixation because of different sequence context.
#' @param paths a \code{sitePath} object
#' @param minSizeBefore minimum tree tips involved before mutation
#' @param minSizeAfter minimum tree tips involved after mutation
#' @param toleranceBefore maximum amino acid variation before mutation
#' @param toleranceAfter maximum amino acid variation after mutation
#' @return 
#' \code{fixationSites} returns a list of mutations with names of the tips involved.
#' The name of each list element is the discovered mutation. A mutation has two vectors of
#' tip names: 'from' before the fixation and 'to' after the fixation.
#' @export
fixationSites.sitePath <- function(
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
  refSeqName <- reference
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
  attr(mutations, "tree") <- tree
  attr(mutations, "align") <- align
  attr(mutations, "reference") <- reference
  attr(mutations, "refSeqName") <- refSeqName
  class(mutations) <- "fixationSites"
  return(mutations)
}

#' @export
fixationSites <- function(x, ...) UseMethod("fixationSites")

#' @export
print.fixationSites <- function(fixationSites) {
  # for (m in names(fixationSites)) {
  #   cat(m)
  #   cat("\n")
  #   from <- fixationSites[[m]]$from
  #   nfrom <- length(from)
  #   cat(" from ", nfrom, " tips", sep = "")
  #   if (nfrom >= 6) {
  #     cat(": ", paste(from[1:5], collapse = ", "), " ...\n", sep = "")
  #   } else {
  #     cat(": ", paste(from, collapse = ", "), "\n", sep = "")
  #   }
  #   to <- fixationSites[[m]]$to
  #   nto <- length(to)
  #   cat(" to ", nto, " tips", sep = "")
  #   if (nto >= 6) {
  #     cat(":", paste(to[1:5], collapse = ", "), " ...\n", sep = "")
  #   } else {
  #     cat(":", paste(to, collapse = ", "), "\n", sep = "")
  #   }
  #   cat("\n")
  # }
  cat(paste(names(fixationSites), collapse = " "))
  refSeqName <- attr(fixationSites, "refSeqName")
  if (is.null(refSeqName)) {
    cat("\nNo reference sequence specified. Using alignment numbering\n")
  } else {
    cat("\nReference sequence: ", refSeqName, "\n", sep = "")
  }
}
