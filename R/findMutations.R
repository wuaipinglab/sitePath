#' @rdname Pre-assessment
#' @name Pre-assessment
#' @title Things can be done before the analysis
#' @description \code{similarityMatrix} calculates similarity between aligned sequences
#' The similarity matrix can be used in \code{\link{groupTips}} or \code{\link{sitePath}}
#' @param tree 
#' a \code{phylo} object. This commonly can be from tree paring function 
#' in \code{ape} or \code{ggtree}
#' @param align 
#' an \code{alignment} object. This commonly can be from sequence parsing function
#' in \code{ape} or \code{seqinr} and many others
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
#' @param step the similarity window for calculating and ploting
#' @param maxPath 
#' maximum number of path to show in the plot. The number of path 
#' in the raw tree can be far greater than trimmed tree. To better 
#' see the impact of chaning threshold on path number. This should be
#' specified. The default is one 20th of tree tip number.
#' @param makePlot whether make a dot plot when return
#' @return 
#' \code{sneakPeek} return the similarity threhold against number of sitePath.
#' There will be a simple dot plot between threshold and path number if
#' \code{makePlot} is TRUE.
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
      stop("Invalid maxPath")
    }
  }
  similarity <- numeric(0)
  pathNum <- integer(0)
  for (s in seq(1, minSim, step)) {
    paths <- sitePath(tree, align, s, simMatrix, FALSE)
    if (maxPath < length(paths)) {
      next
    } else if (length(paths) <= 1) {
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
#' @param tree 
#' a \code{phylo} object. This commonly can be from tree paring function 
#' in \code{ape} or \code{ggtree}
#' @param align 
#' an \code{alignment} object. This commonly can be from sequence parsing function
#' in \code{ape} or \code{seqinr} and many others
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
#' @param paths a \code{sitePath} object returned from \code{\link{sitePath}} function
#' @param tolerance
#' A vector of two integers to specify maximum amino acid variation before/after mutation. 
#' Otherwise the mutation will not be counted into the return.
#' If more than one number is given, the ancestral takes the first and
#' descendant takes the second as the maximum.
#' If only given one number, it' <- s the maximum for both ancestral and descendant.
#' @param minEffectiveSize
#' A vector of two integers to specifiy minimum tree tips involved before/after mutation. 
#' Otherwise the mutation will not be counted into the return.
#' If more than one number is given, the ancestral takes the first and
#' descendant takes the second as the minimum.
#' If only given one number, it's the minimum for both ancestral and descendant. 
#' @param extendedSearch 
#' Whether to extend the search. The terminal of each \code{sitePath} is
#' a cluster of tips. To look for the fixation mutation in the cluster, 
#' the common ancestral node of farthest tips (at least two) will be
#' the new terminal search point.
#' @return 
#' \code{fixationSites} returns a list of mutations with names of the tips involved.
#' The name of each list element is the discovered mutation. A mutation has two vectors of
#' tip names: 'from' before the fixation and 'to' after the fixation.
#' @export
fixationSites.sitePath <- function(
  paths, reference = NULL, gapChar = '-',
  tolerance = c(0, 0), minEffectiveSize = NULL,
  extendedSearch = TRUE
) {
  tree <- attr(paths, "tree")
  align <- attr(paths, "align")
  refSeqName <- reference
  if (is.null(reference)) {reference <- 1:nchar(align[1])} else {
    if (!is.character(gapChar) || nchar(gapChar) != 1 || length(gapChar) != 1) {
      stop("gapChar only accepts one single character")
    }
    reference <- getReference(align[which(tree$tip.label == reference)], gapChar)
  }
  if (!is.numeric(tolerance)) {
    stop("tolerance only accepts numeric")
  } else {
    toleranceAnc <- tolerance[1]
    toleranceDesc <- if (length(tolerance) == 1) toleranceAnc else tolerance[2]
  }
  if (is.null(minEffectiveSize)) {
    minAnc <- length(tree$tip.label) / 10
    minDesc <- length(tree$tip.label) / 10
  } else if (!is.numeric(minEffectiveSize)) {
    stop("minEffectiveSize only accepts numeric")
  } else {
    minAnc <- minEffectiveSize
    minDesc <- if (length(minEffectiveSize) == 1) minAnc else minEffectiveSize[2]
  }
  divNodes <- unique(divergentNode(paths))
  if (extendedSearch) {
    paths <- lapply(paths, function(p) {
      sn <- p[length(p)]
      extended <- lapply(ChildrenTips(tree, sn), function(t) nodepath(tree, sn, t)[-1])
      el <- lengths(extended)
      ml <- max(el)
      longest <- extended[which(el == ml)]
      extended <- lapply(1:ml, function(i) unique(sapply(longest, "[[", i)))
      c(p, unlist(extended[which(lengths(extended) == 1)]))
    })
  }
  mutations <- list()
  for (minLen in 2:max(lengths(paths))) { # literate all sitePath at the same time
    for (path in unique(ancestralPaths(paths, minLen))) {
      afterTips <- ChildrenTips(tree, tail(path, 1))
      if (length(afterTips) < minDesc) {
        next
      }
      pathBefore <- path[1:(length(path) - 1)]
      excludedTips <- sapply(pathBefore[which(pathBefore %in% divNodes)], function(node) {
        children <- tree$edge[which(tree$edge[, 1] == node), 2]
        children <- children[which(children > length(tree$tip.label) & !children %in% path)]
        return(ChildrenTips(tree, children))
      })
      beforeTips <- which(!1:length(tree$tip.label) %in% c(afterTips, unlist(excludedTips)))
      if (length(beforeTips) < minAnc && length(excludedTips) == 0) {
        next
      }
      after <- align[afterTips]
      before <- align[beforeTips]
      for (i in 1:length(reference)) {
        s <- reference[i] - 1
        b <- summarizeAA(before, s, toleranceAnc)
        a <- summarizeAA(after, s, toleranceDesc)
        if (!(is.na(b) || is.na(a)) && a != b) {
          mut <- paste(b, i, a, sep = "")
          from <- tree$tip.label[beforeTips]
          to <- tree$tip.label[afterTips]
          if (!mut %in% names(mutations) || length(to) > length(mutations[[mut]]$to)) {
            mutations[[mut]] <- list(ancestral = from, descendant = to)
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
as.data.frame.fixationSites <- function(mutations) {
  tree <- attr(mutations, "tree")
  reference <- attr(mutations, "reference")
  res <- as.data.frame(matrix(
    nrow = length(tree$tip.label),
    ncol = length(reference),
    dimnames = list(tree$tip.label, 1:length(reference))
  ))
  for (m in names(mutations)) {
    site <- as.integer(substr(m, 2, nchar(m) - 1))
    fixedAA <- substr(m, nchar(m), nchar(m))
    for (desc in mutations[[m]][[2]]) {
      res[desc, site] <- AA_FULL_NAMES[tolower(fixedAA)]
    }
  }
  whichNA <- is.na(res)
  res <- res[rowSums(whichNA) < ncol(res), colSums(whichNA) < nrow(res)]
  return(res)
}

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
  if (length(fixationSites) == 0) {
    cat("No mutation found\n")
  } else {
    cat(paste(names(fixationSites), collapse = " "), "\n")
    refSeqName <- attr(fixationSites, "refSeqName")
    if (is.null(refSeqName)) {
      cat("No reference sequence specified. Using alignment numbering\n")
    } else {
      cat("\nReference sequence: ", refSeqName, "\n", sep = "")
    }
  }
}
