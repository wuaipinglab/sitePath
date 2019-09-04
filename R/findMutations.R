#' @rdname findSites
#' @name findSites
#' @title Finding sites with variation
#' @description
#' Single nucleotide polymorphism (SNP) in the whole package refers to
#' variation of amino acid. \code{findSNPsite} will try to find SNP in
#' the multiple sequence alignment. A reference sequence
#' and gap character may be specified to number the site. This is
#' irrelevant to the intended analysis but might be helpful to evaluate
#' the performance of \code{fixationSites}.
#' @param tree The return from \code{\link{addMSA}} function
#' @param minSNP Minimum number of amino acid variation to be a SNP
#' @examples
#' data('zikv_tree_reduced')
#' data('zikv_align_reduced')
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' SNPsites(tree)
#' @return \code{SNPsite} returns a list of qualified SNP site
#' @export
SNPsites <- function(tree, minSNP = NULL) {
    if (is.null(minSNP)) {
        minSNP <- length(tree$tip.label)/10
    }
    alignedSeq <- attr(tree, "align")
    if (is.null(alignedSeq)) {
        stop("No alignment found in \"tree\"")
    }
    seqLen <- unique(nchar(alignedSeq))
    if (length(seqLen) != 1) 
        stop("Sequence length not equal")
    reference <- attr(tree, "reference")
    qualified <- integer(0)
    alignedSeq <- strsplit(alignedSeq, "")
    for (i in seq_along(reference)) {
        SNP <- table(vapply(alignedSeq, FUN = "[[", FUN.VALUE = character(1), reference[i]))
        if (sum(SNP > minSNP) >= 2) {
            qualified <- c(qualified, i)
        }
    }
    return(qualified)
}

#' @export
print.sitePath <- function(x, ...) {
    cat("Site", attr(x, "site"), "may experience fixation on", length(x), "path(s):\n\n")
    # A 'sitePath' composes of all the fixation paths for a single site. So each 'm'
    # represent a single fixation path
    for (m in x) {
        if (length(m) == 2) {
            mutName <- paste0(attr(m[[1]], "AA"), attr(x, "site"), attr(m[[2]], "AA"))
            cat(mutName, paste0("(", length(m[[1]]), "->", length(m[[2]]), ")"), 
                "\n")
        } else {
            mutName <- character(0)
            for (tips in m) {
                aa <- attr(tips, "AA")
                mutName <- c(mutName, paste0(aa, "(", length(tips), ")"))
            }
            cat(paste0(mutName, collapse = " -> "), "\n")
        }
    }
    cat("\nIn the bracket are the number of tips", "involved before and after the fixation\n")
}

.compareMutPathAA <- function(e, s) {
    len <- length(e)
    if (len != length(s)) {
        return(FALSE)
    } else {
        for (i in seq_len(len)) {
            if (attr(e[[i]], "AA") != attr(s[[i]], "AA")) {
                return(FALSE)
            }
        }
    }
    return(TRUE)
}

#' @rdname findSites
#' @name fixationSites
#' @description
#' After finding the \code{\link{lineagePath}} of a phylogenetic tree,
#' \code{fixationSites} uses the result to find those sites that show
#' fixation on some, if not all, of the lineages. Parallel evolution is
#' relatively common in RNA virus. There is chance that some site be fixed
#' in one lineage but does not show fixation because of different
#' sequence context.
#' @param paths
#' a \code{lineagePath} object returned from \code{\link{lineagePath}} function
#' @param minEffectiveSize
#' A vector of two integers to specifiy minimum tree tips involved
#' before/after mutation. Otherwise the mutation will not be counted into
#' the return. If more than one number is given, the ancestral takes the first
#' and descendant takes the second as the minimum. If only given one number,
#' it's the minimum for both ancestral and descendant.
#' @param searchDepth
#' The function uses heuristic search but the termination of the search
#' cannot be intrinsically decided. \code{searchDepth} is needed to tell
#' the search when to stop.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' fixationSites(lineagePath(tree))
#' @return
#' \code{fixationSites} returns a list of mutations
#' with names of the tips involved. The name of each list element
#' is the discovered mutation. A mutation has two vectors of tip names:
#' 'from' before the fixation and 'to' after the fixation.
#' @importFrom utils tail
#' @export
fixationSites.lineagePath <- function(paths, minEffectiveSize = NULL, searchDepth = 1, 
    ...) {
    tree <- attr(paths, "tree")
    nTips <- length(tree$tip.label)
    align <- attr(paths, "align")
    # Generate the site mapping from reference
    reference <- attr(paths, "reference")
    # Exclude the invariant sites
    loci <- which(vapply(X = seq_along(reference), FUN = function(s) {
        length(unique(substr(align, s, s))) > 1
    }, FUN.VALUE = logical(1)))
    # Get the 'minEffectiveSize' for each fixation
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- nTips/length(unique(unlist(paths)))
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    minEffectiveSize <- ceiling(minEffectiveSize)
    # Get the 'searchDepth' for heuristic search
    if (searchDepth < 1) {
        stop("\"searchDepth\" should be at least 1")
    } else {
        searchDepth <- ceiling(searchDepth)
    }
    divNodes <- unique(divergentNode(paths))
    paths <- .extendPaths(paths, tree)
    # Get all the nodes that are not at divergent point
    nodes <- setdiff(unlist(paths), divNodes)
    # Get the sequence of the children tips that are descendant of the nodes.  Assign
    # the tip index to the sequences for retrieving the tip name
    nodeAlign <- lapply(nodes, function(n) {
        childrenNode <- tree$edge[which(tree$edge[, 1] == n), 2]
        childrenNode <- setdiff(childrenNode, c(divNodes, nodes))
        tips <- .childrenTips(tree, childrenNode)
        res <- align[tips]
        names(res) <- tips
        return(res)
    })
    # Assign the node names to the 'nodeAlign' list
    names(nodeAlign) <- nodes
    # Store all the tips by node and their fixed AA to avoid repeating calculation.
    # Each entry in the list stores the info for a single site. Under each site are
    # groups of 'nodeTips' where tip names and the summary of their AAs at the site.
    nodeAAsum <- list()
    # 'res' is going to be the return of this function. Each entry in the list is the
    # 'sitePath' for a site. Each site ('sitePath') consists of 'mutPath' that is
    # named by the starting node name.  The fixed AA and number of non-dominant AA is
    # also stored.
    res <- list()
    # Iterate each path in the 'lineagePath' object
    for (path in paths) {
        # Try every node as the terminal node
        for (maxLen in seq_along(path)[-1]) {
            # For the sequences after the terminal node, examine them as a whole
            afterTips <- as.integer(.childrenTips(tree, path[maxLen]))
            # Group the sequences by nodes in the path before
            pathBefore <- setdiff(path[seq_len(maxLen - 1)], divNodes)
            # Iterate every site
            for (i in loci) {
                s <- reference[i] - 1
                # 'tableAA' is similar to R function 'table' Here the AA at site 's' for all tips
                # is summarized
                afterSummary <- tableAA(align[afterTips], s)
                # TODO: 'tolerance' is used to track total number of non-dominant AA for the site
                # initialized with the number of non-dominant AA in the 'afterTips'.  May jump to
                # next site if exceeding tolerance
                tolerance <- sum(afterSummary) - max(afterSummary)
                if (tolerance > length(afterTips) * 0.01) {
                  next
                }
                attr(afterTips, "aaSummary") <- afterSummary
                # If AA is purely fixed for the node
                if (length(afterSummary) == 1) {
                  afterAA <- names(afterSummary)
                } else {
                  afterAA <- NULL
                }
                site <- as.character(i)
                # 'nodeTips' is vector of tips with an attribute of 'AA' to store fixed AA, and
                # an attribute of 'nonDominant' to store non-dominant AA.
                nodeTips <- integer(0)
                # We need a 'previousAA' and a 'currentAA' to track the AA along the
                # 'lineagePath'. In the summary stage, we only focus on the purely fixed AA.
                previousAA <- NULL
                currentAA <- NULL
                previousNode <- NULL
                # 'nodeSummaries' groups tips by node in the summary stage. In the case that the
                # adjacent ndoes have the same AA purely fixed, they will be grouped into one.
                nodeSummaries <- list()
                # Iterate each single node in the 'pathBefore'. This is the summary stage
                for (node in as.character(pathBefore)) {
                  # If the node has a record in nodeAAsum
                  nodeTips <- nodeAAsum[[site]][[node]]
                  # Summarize the node if not existed
                  if (is.null(nodeTips)) {
                    # Get the related descendant tips from 'nodeAlign'
                    nodeTips <- as.integer(names(nodeAlign[[node]]))
                    aaSummary <- tableAA(nodeAlign[[node]], s)
                    # Assign the 'aaSummary' to the tip names
                    attr(nodeTips, "aaSummary") <- aaSummary
                    # Store the result to avoid repeating calculation TODO: Bug fix: 'nodeTips' will
                    # be allocated as a named vector if its length is of 1. An S4 class should be
                    # created for 'nodeAAsum' to solve this typing problem. Here a vector of c(1,2)
                    # is used to first guarantee a type of list.
                    nodeAAsum[[site]][[node]] <- c(1, 2)
                    nodeAAsum[[site]][[node]] <- nodeTips
                  }
                  # Extract 'aaSummary' if existed
                  aaSummary <- attr(nodeTips, "aaSummary")
                  # If AA is purely fixed for the node
                  if (length(aaSummary) == 1) {
                    currentAA <- names(aaSummary)
                  } else {
                    currentAA <- NULL
                  }
                  # Attach the node to the 'preivous' node if they're both purely fixed and the
                  # fixed AA is the same
                  if (!is.null(previousAA) && !is.null(currentAA) && previousAA == 
                    currentAA) {
                    node <- previousNode
                    nodeTips <- c(nodeSummaries[[node]], nodeTips)
                    attr(nodeTips, "aaSummary") <- attr(nodeSummaries[[node]], "aaSummary") + 
                      aaSummary
                    # Oddly, R uses the name of the first element in the numeric vector when adding
                    # two named number. So there is no name (AA) assignment
                  }
                  # Assign or re-assign the nodeTips with 'aaSummary' to the 'nodeSummaries'
                  nodeSummaries[[node]] <- nodeTips
                  previousAA <- currentAA
                  previousNode <- node
                }
                # Attach the 'afterTips' to 'nodeSummaries'
                if (!is.null(previousAA) && !is.null(afterAA) && previousAA == afterAA) {
                  nodeTips <- c(nodeSummaries[[previousNode]], afterTips)
                  attr(nodeTips, "aaSummary") <- attr(nodeSummaries[[previousNode]], 
                    "aaSummary") + afterSummary
                  nodeSummaries[[previousNode]] <- nodeTips
                } else {
                  nodeSummaries[[as.character(path[maxLen])]] <- afterTips
                }
                # Skip to the next 'site' if AA of 'pathBefore' is purely fixed or the terminal
                # AA is the same as fixed AA of 'afterTips'. This avoids repetition.
                if (length(nodeSummaries) <= 1) {
                  next
                } else if (site %in% names(res)) {
                  aTips <- nodeSummaries[[length(nodeSummaries)]]
                  exist <- vapply(res[[site]], FUN = function(ep) {
                    all(aTips %in% ep[[length(ep)]])
                  }, FUN.VALUE = logical(1))
                  if (any(exist)) {
                    next
                  }
                }
                seg <- minimizeEntropy(nodeSummaries, minEffectiveSize, searchDepth)
                targetIndex <- NULL
                if (length(seg) < 2) {
                  next
                } else if (!site %in% names(res)) {
                  targetIndex <- 1
                } else {
                  # Some site may have multiple fixation on multiple lineages. The following is for
                  # deciding at which index should it be assigned in the 'res[[site]]' Retrieve the
                  # existing mutation path of the site
                  existPath <- res[[site]]
                  # Add new fixation for the site. Assume none of the existing mutation path has
                  # the same mutations as 'seg'
                  targetIndex <- length(existPath) + 1
                  # Which mutaiton path has the same mutations as 'seg'
                  existIndex <- which(vapply(existPath, FUN = .compareMutPathAA, 
                    FUN.VALUE = logical(1), seg))
                  if (length(existIndex) > 0) {
                    # afterTips <- unlist(tail(seg, 1)) 'Adding state' for each existing mutation
                    # path that has the same mutations as 'seg' does.  0L: 'seg' is a subset of
                    # existing path and won't be added to the 'res' 1L: The path is different from
                    # 'seg'. They just happen to share the same mutations. A new mutation path will
                    # be created 2L: The existing path is a subset of 'seg' and needs to be replaced
                    # by 'seg'
                    qualified <- vapply(existPath[existIndex], FUN = function(ep) {
                      existTips <- unlist(tail(ep, 1))
                      if (all(existTips %in% afterTips)) {
                        return(2L)
                      }
                      if (!all(afterTips %in% existTips)) {
                        return(1L)
                      }
                      return(0L)
                    }, FUN.VALUE = integer(1))
                    r <- qualified == 2L
                    if (any(r)) {
                      # Remove existing paths with an 'adding state' of 2L and add the 'seg'
                      res[[site]] <- res[[site]][-which(r)]
                      targetIndex <- length(res[[site]]) + 1
                    } else if (all(qualified == 0)) {
                      targetIndex <- NULL
                    }
                  }
                }
                if (is.null(targetIndex)) {
                  next
                }
                # Assign the result to the 'res[[site]]'
                res[[site]][[targetIndex]] <- seg
                # Assign or re-assign class and 'site' to the 'res[[site]]'
                class(res[[site]]) <- "sitePath"
                attr(res[[site]], "site") <- i
            }
        }
    }
    attr(res, "paths") <- paths
    attr(res, "reference") <- reference
    class(res) <- "multiFixationSites"
    return(res)
}

#' @export
fixationSites <- function(paths, ...) UseMethod("fixationSites")

#' @export
print.fixationSites <- function(x, ...) {
    cat("Result for", length(attr(x, "paths")), "paths:\n\n")
    if (length(x) == 0) {
        cat("No multi-fixation found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(x, "refSeqName")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified. Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}

#' @rdname findSites
#' @name multiFixationSites
#' @description
#' After finding the \code{\link{lineagePath}} of a phylogenetic tree,
#' \code{multiFixationSites} uses the result to find those sites that show
#' multiple fixations on some, if not all, of the lineages.
#' @examples
#' data(h3n2_tree_reduced)
#' data(h3n2_align_reduced)
#' tree <- addMSA(h3n2_tree_reduced, alignment = h3n2_align_reduced)
#' multiFixationSites(lineagePath(tree))
#' @return
#' \code{multiFixationSites} returns sites with multiple fixations.
#' @export
multiFixationSites.lineagePath <- function(paths, minEffectiveSize = NULL, searchDepth = 1, 
    ...) {
    tree <- attr(paths, "tree")
    nTips <- length(tree$tip.label)
    align <- attr(paths, "align")
    # Generate the site mapping from reference
    reference <- attr(paths, "reference")
    # Exclude the invariant sites
    loci <- which(vapply(X = seq_along(reference), FUN = function(s) {
        length(unique(substr(align, s, s))) > 1
    }, FUN.VALUE = logical(1)))
    # Get the 'minEffectiveSize' for each fixation
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- nTips/length(unique(unlist(paths)))
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    minEffectiveSize <- ceiling(minEffectiveSize)
    # Get the 'searchDepth' for heuristic search
    if (searchDepth < 1) {
        stop("\"searchDepth\" should be at least 1")
    } else {
        searchDepth <- ceiling(searchDepth)
    }
    divNodes <- unique(divergentNode(paths))
    paths <- .extendPaths(paths, tree)
    # Get all the nodes that are not at divergent point
    nodes <- setdiff(unlist(paths), divNodes)
    # Get the sequence of the children tips that are descendant of the nodes.  Assign
    # the tip index to the sequences for retrieving the tip name
    nodeAlign <- lapply(nodes, function(n) {
        childrenNode <- tree$edge[which(tree$edge[, 1] == n), 2]
        childrenNode <- setdiff(childrenNode, c(divNodes, nodes))
        tips <- .childrenTips(tree, childrenNode)
        res <- align[tips]
        names(res) <- tips
        return(res)
    })
    # Assign the node names to the 'nodeAlign' list
    names(nodeAlign) <- nodes
    # Store all the tips by node and their fixed AA to avoid repeating calculation.
    # Each entry in the list stores the info for a single site. Under each site are
    # groups of 'nodeTips' where tip names and the summary of their AAs at the site.
    nodeAAsum <- list()
    # 'res' is going to be the return of this function. Each entry in the list is the
    # 'sitePath' for a site. Each site ('sitePath') consists of 'mutPath' that is
    # named by the starting node name.  The fixed AA and number of non-dominant AA is
    # also stored.
    res <- list()
    # Iterate each path in the 'lineagePath' object
    for (path in paths) {
        # Try every node as the terminal node
        for (maxLen in seq_along(path)[-1]) {
            # For the sequences after the terminal node, examine them as a whole
            afterTips <- as.integer(.childrenTips(tree, path[maxLen]))
            # Group the sequences by nodes in the path before
            pathBefore <- setdiff(path[seq_len(maxLen - 1)], divNodes)
            # Iterate every site
            for (i in loci) {
                s <- reference[i] - 1
                # 'tableAA' is similar to R function 'table' Here the AA at site 's' for all tips
                # is summarized
                afterSummary <- tableAA(align[afterTips], s)
                # TODO: 'tolerance' is used to track total number of non-dominant AA for the site
                # initialized with the number of non-dominant AA in the 'afterTips'.  May jump to
                # next site if exceeding tolerance
                tolerance <- sum(afterSummary) - max(afterSummary)
                if (tolerance > length(afterTips) * 0.01) {
                  next
                }
                attr(afterTips, "aaSummary") <- afterSummary
                # If AA is purely fixed for the node
                if (length(afterSummary) == 1) {
                  afterAA <- names(afterSummary)
                } else {
                  afterAA <- NULL
                }
                site <- as.character(i)
                # 'nodeTips' is vector of tips with an attribute of 'AA' to store fixed AA, and
                # an attribute of 'nonDominant' to store non-dominant AA.
                nodeTips <- integer(0)
                # We need a 'previousAA' and a 'currentAA' to track the AA along the
                # 'lineagePath'. In the summary stage, we only focus on the purely fixed AA.
                previousAA <- NULL
                currentAA <- NULL
                previousNode <- NULL
                # 'nodeSummaries' groups tips by node in the summary stage. In the case that the
                # adjacent ndoes have the same AA purely fixed, they will be grouped into one.
                nodeSummaries <- list()
                # Iterate each single node in the 'pathBefore'. This is the summary stage
                for (node in as.character(pathBefore)) {
                  # If the node has a record in nodeAAsum
                  nodeTips <- nodeAAsum[[site]][[node]]
                  # Summarize the node if not existed
                  if (is.null(nodeTips)) {
                    # Get the related descendant tips from 'nodeAlign'
                    nodeTips <- as.integer(names(nodeAlign[[node]]))
                    aaSummary <- tableAA(nodeAlign[[node]], s)
                    # Assign the 'aaSummary' to the tip names
                    attr(nodeTips, "aaSummary") <- aaSummary
                    # Store the result to avoid repeating calculation TODO: Bug fix: 'nodeTips' will
                    # be allocated as a named vector if its length is of 1. An S4 class should be
                    # created for 'nodeAAsum' to solve this typing problem. Here a vector of c(1,2)
                    # is used to first guarantee a type of list.
                    nodeAAsum[[site]][[node]] <- c(1, 2)
                    nodeAAsum[[site]][[node]] <- nodeTips
                  }
                  # Extract 'aaSummary' if existed
                  aaSummary <- attr(nodeTips, "aaSummary")
                  # If AA is purely fixed for the node
                  if (length(aaSummary) == 1) {
                    currentAA <- names(aaSummary)
                  } else {
                    currentAA <- NULL
                  }
                  # Attach the node to the 'preivous' node if they're both purely fixed and the
                  # fixed AA is the same
                  if (!is.null(previousAA) && !is.null(currentAA) && previousAA == 
                    currentAA) {
                    node <- previousNode
                    nodeTips <- c(nodeSummaries[[node]], nodeTips)
                    attr(nodeTips, "aaSummary") <- attr(nodeSummaries[[node]], "aaSummary") + 
                      aaSummary
                    # Oddly, R uses the name of the first element in the numeric vector when adding
                    # two named number. So there is no name (AA) assignment
                  }
                  # Assign or re-assign the nodeTips with 'aaSummary' to the 'nodeSummaries'
                  nodeSummaries[[node]] <- nodeTips
                  previousAA <- currentAA
                  previousNode <- node
                }
                # Attach the 'afterTips' to 'nodeSummaries'
                if (!is.null(previousAA) && !is.null(afterAA) && previousAA == afterAA) {
                  nodeTips <- c(nodeSummaries[[previousNode]], afterTips)
                  attr(nodeTips, "aaSummary") <- attr(nodeSummaries[[previousNode]], 
                    "aaSummary") + afterSummary
                  nodeSummaries[[previousNode]] <- nodeTips
                } else {
                  nodeSummaries[[as.character(path[maxLen])]] <- afterTips
                }
                # Skip to the next 'site' if AA of 'pathBefore' is purely fixed or the terminal
                # AA is the same as fixed AA of 'afterTips'. This avoids repetition.
                if (length(nodeSummaries) <= 1) {
                  next
                } else if (site %in% names(res)) {
                  aTips <- nodeSummaries[[length(nodeSummaries)]]
                  exist <- vapply(res[[site]], FUN = function(ep) {
                    all(aTips %in% ep[[length(ep)]])
                  }, FUN.VALUE = logical(1))
                  if (any(exist)) {
                    next
                  }
                }
                seg <- minimizeEntropy(nodeSummaries, minEffectiveSize, searchDepth)
                targetIndex <- NULL
                if (length(seg) < 2) {
                  next
                } else if (!site %in% names(res)) {
                  targetIndex <- 1
                } else {
                  # Some site may have multiple fixation on multiple lineages. The following is for
                  # deciding at which index should it be assigned in the 'res[[site]]' Retrieve the
                  # existing mutation path of the site
                  existPath <- res[[site]]
                  # Add new fixation for the site. Assume none of the existing mutation path has
                  # the same mutations as 'seg'
                  targetIndex <- length(existPath) + 1
                  # Which mutaiton path has the same mutations as 'seg'
                  existIndex <- which(vapply(existPath, FUN = .compareMutPathAA, 
                    FUN.VALUE = logical(1), seg))
                  if (length(existIndex) > 0) {
                    # afterTips <- unlist(tail(seg, 1)) 'Adding state' for each existing mutation
                    # path that has the same mutations as 'seg' does.  0L: 'seg' is a subset of
                    # existing path and won't be added to the 'res' 1L: The path is different from
                    # 'seg'. They just happen to share the same mutations. A new mutation path will
                    # be created 2L: The existing path is a subset of 'seg' and needs to be replaced
                    # by 'seg'
                    qualified <- vapply(existPath[existIndex], FUN = function(ep) {
                      existTips <- unlist(tail(ep, 1))
                      if (all(existTips %in% afterTips)) {
                        return(2L)
                      }
                      if (!all(afterTips %in% existTips)) {
                        return(1L)
                      }
                      return(0L)
                    }, FUN.VALUE = integer(1))
                    r <- qualified == 2L
                    if (any(r)) {
                      # Remove existing paths with an 'adding state' of 2L and add the 'seg'
                      res[[site]] <- res[[site]][-which(r)]
                      targetIndex <- length(res[[site]]) + 1
                    } else if (all(qualified == 0)) {
                      targetIndex <- NULL
                    }
                  }
                }
                if (is.null(targetIndex)) {
                  next
                }
                # Assign the result to the 'res[[site]]'
                res[[site]][[targetIndex]] <- seg
                # Assign or re-assign class and 'site' to the 'res[[site]]'
                class(res[[site]]) <- "sitePath"
                attr(res[[site]], "site") <- i
            }
        }
    }
    attr(res, "paths") <- paths
    attr(res, "reference") <- reference
    class(res) <- "multiFixationSites"
    return(res)
}

#' @export
multiFixationSites <- function(paths, ...) UseMethod("multiFixationSites")

#' @export
print.multiFixationSites <- function(x, ...) {
    cat("Result for", length(attr(x, "paths")), "paths:\n\n")
    if (length(x) == 0) {
        cat("No multi-fixation found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(x, "refSeqName")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified. Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}
