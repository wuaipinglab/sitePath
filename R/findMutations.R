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
#' @param reference
#' Name of reference for site numbering. The name has to be one of the
#' sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar
#' The character to indicate gap. The numbering will skip the gapChar
#' for the reference sequence.
#' @param minSNP Minimum number of amino acid variation to be a SNP
#' @examples
#' data("zikv_tree_reduced")
#' data("zikv_align_reduced")
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' SNPsites(tree)
#' @return \code{SNPsite} returns a list of qualified SNP site
#' @export
SNPsites <- function(tree,
                     reference = NULL,
                     gapChar = '-',
                     minSNP = NULL) {
    if (is.null(minSNP)) {
        minSNP <- length(tree$tip.label) / 10
    }
    alignedSeq <- attr(tree, "alignment")
    if (is.null(alignedSeq)) {
        stop("No alignment found in \"tree\"")
    }
    seqLen <- unique(nchar(alignedSeq))
    if (length(seqLen) != 1)
        stop("Sequence length not equal")
    reference <-
        checkReference(tree, alignedSeq, reference, gapChar)
    qualified <- integer(0)
    alignedSeq <- strsplit(alignedSeq, "")
    for (i in  seq_along(reference)) {
        SNP <- table(vapply(
            alignedSeq,
            FUN = "[[",
            FUN.VALUE = character(1),
            reference[i]
        ))
        if (sum(SNP > minSNP) >= 2) {
            qualified <- c(qualified, i)
        }
    }
    return(qualified)
}

#' @export
print.sitePath <- function(x, ...) {
    cat("Site",
        attr(x, "site"),
        "may experience fixation on",
        length(x),
        "path(s):\n\n")
    # A "sitePath" composes of all the fixation paths for
    # a single site. So each "m" represent a single fixation
    # path
    for (m in x) {
        if (length(m) == 2) {
            mutName <- paste0(attr(m[[1]], "AA"),
                              attr(x, "site"),
                              attr(m[[2]], "AA"))
            cat(mutName,
                paste0("(", length(m[[1]]), "->", length(m[[2]]), ")"),
                "\n")
        } else {
            mutName <- character(0)
            for (tips in m) {
                aa <- attr(tips, "AA")
                mutName <-
                    c(mutName, paste0(aa, "(", length(tips), ")"))
            }
            cat(paste0(mutName, collapse = " -> "), "\n")
        }
    }
    cat("\nIn the bracket are the number of tips",
        "involved before and after the fixation\n")
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
#' @param tolerance
#' A vector of two integers to specify maximum amino acid variation
#' before/after mutation. Otherwise the mutation will not be counted into
#' the return. If more than one number is given, the ancestral takes the first
#' and descendant takes the second as the maximum. If only given one number,
#' it's the maximum for both ancestral and descendant. The default is 0.01
#' @param minEffectiveSize
#' A vector of two integers to specifiy minimum tree tips involved
#' before/after mutation. Otherwise the mutation will not be counted into
#' the return. If more than one number is given, the ancestral takes the first
#' and descendant takes the second as the minimum. If only given one number,
#' it's the minimum for both ancestral and descendant.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' fixationSites(
#'     lineagePath(tree),
#'     tolerance = c(1, 1),
#'     minEffectiveSize = c(50, 50)
#' )
#' @return
#' \code{fixationSites} returns a list of mutations
#' with names of the tips involved. The name of each list element
#' is the discovered mutation. A mutation has two vectors of tip names:
#' 'from' before the fixation and 'to' after the fixation.
#' @importFrom utils tail
#' @export
fixationSites.lineagePath <- function(paths,
                                      reference = NULL,
                                      gapChar = '-',
                                      tolerance = 0.01,
                                      minEffectiveSize = NULL,
                                      ...) {
    tree <- attr(paths, "tree")
    align <- attr(paths, "align")
    if (!is.numeric(tolerance)) {
        stop("\"tolerance\" only accepts numeric")
    } else if (any(tolerance < 0)) {
        stop("\"tolerance\" can only be positive number")
    } else {
        toleranceAnc <- tolerance[1]
        toleranceDesc <-
            if (length(tolerance) == 1) {
                toleranceAnc
            } else {
                tolerance[2]
            }
    }
    if (is.null(minEffectiveSize)) {
        minAnc <- length(tree$tip.label) / 10
        minDesc <- length(tree$tip.label) / 10
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    } else if (any(minEffectiveSize <= 0)) {
        stop("\"minEffectiveSize\" can only be positive number")
    } else {
        minAnc <- minEffectiveSize[1]
        minDesc <-
            if (length(minEffectiveSize) == 1) {
                minAnc
            } else {
                minEffectiveSize[2]
            }
    }
    refSeqName <- reference
    reference <- checkReference(tree, align, reference, gapChar)
    divNodes <- unique(divergentNode(paths))
    paths <- extendPaths(paths, tree)
    # "mutations" is the final return of this function.
    # Each item ("sitePath") is named by the site predicted to
    # experience fixation and an attribute names "site" also
    # stores the site index. The info for a fixation ("mutPath")
    # consists of a list of grouped tips with their AA as the
    # name and "AA" attribute.
    mutations <- list()
    for (minLen in 2:max(lengths(paths))) {
        # Iterate all lineagePath at the same time
        for (path in unique(ancestralPaths(paths, minLen))) {
            # First, the number of tips should meet the minimum size
            # otherwise we'll go for the next node
            afterTips <- ChildrenTips(tree, tail(path, 1))
            if (length(afterTips) < minDesc) {
                next
            }
            pathBefore <- path[seq_len(length(path) - 1)]
            # The "excludedTips" are the descendants of the "divNodes"
            # which are not belong to the "lineagePath"
            excludedTips <- lapply(
                intersect(pathBefore, divNodes),
                FUN = function(node) {
                    children <- tree$edge[which(tree$edge[, 1] == node), 2]
                    children <-
                        children[which(children > length(tree$tip.label) &
                                           !children %in% path)]
                    return(ChildrenTips(tree, children))
                }
            )
            beforeTips <-
                setdiff(x = seq_along(tree$tip.label),
                        y = c(afterTips, unlist(excludedTips)))
            if (length(beforeTips) < minAnc &&
                length(excludedTips) == 0) {
                next
            }
            # Get the sequences for ancestral and descendant group
            after <- align[afterTips]
            before <- align[beforeTips]
            # Iterate through each site and compare the two groups
            for (i in seq_along(reference)) {
                s <- reference[i] - 1
                # Get the dominant AA of the site for the two groups
                b <- summarizeAA(before, s, toleranceAnc)
                a <- summarizeAA(after, s, toleranceDesc)
                if (is.na(b) || is.na(a) || a == b) {
                    # The AA should be fixed for both ancestral and
                    # descendant group. Plus the fixed AA can't be
                    # same. Otherwise, jump to the next site.
                    next
                }
                site <- as.character(i)
                # Use the site to retrieve possible existing "sitePath"
                existPath <- mutations[[site]]
                # Get the tips names for two groups
                # from <- tree$tip.label[beforeTips]
                # to <- tree$tip.label[afterTips]
                # Use "targetIndex" to deciden where to insert
                # the tip/AA names in the "sitePath" entry.
                if (is.null(existPath)) {
                    # Create an entry for the site if it's not existed
                    targetIndex <- 1
                } else {
                    # Use "mut" to retrieve possible existing "sitePath".
                    # The reason for not using "mut" as the name of each
                    # mutation is that same mutaiton fixation could
                    # theoretically happen on more than one lineages.
                    # Using "mut" as the name will cause the conflict
                    # and override some fixation events.
                    mut <- paste(c(b, a), collapse = "")
                    existIndex <- which(vapply(
                        existPath,
                        FUN = function(i) {
                            paste(names(i), collapse = "") == mut
                        },
                        FUN.VALUE = logical(1)
                    ))
                    if (length(existIndex) == 0) {
                        # Add new fixation for the site if it's not existed
                        targetIndex <- length(existPath) + 1
                    } else {
                        qualified <- vapply(
                            existIndex,
                            FUN = function(ei) {
                                existMut <- existPath[[ei]]
                                existTips <-
                                    unlist(tail(existMut, 1))
                                if (all(existTips %in% afterTips)) {
                                    return(2L)
                                }
                                if (!all(afterTips %in% existTips)) {
                                    return(1L)
                                }
                                return(0L)
                            },
                            FUN.VALUE = integer(1)
                        )
                        # Extend the "sitePath" when all existing paths
                        # have an "adding state" of 1L, being different
                        # from the "s"
                        if (all(qualified == 1L)) {
                            targetIndex <- length(existPath) + 1
                        } else if (any(r <- qualified == 2L)) {
                            # Remove existing paths with an "adding state"
                            # of 2L and add the "s"
                            mutations[[site]] <-
                                mutations[[site]][-which(r)]
                            targetIndex <-
                                length(mutations[[site]]) + 1
                        }
                    }
                }
                # Construct the fixation path
                attr(beforeTips, "AA") <- b
                attr(afterTips, "AA") <- a
                mutPath <- list(beforeTips, afterTips)
                names(mutPath) <- c(b, a)
                # Store the entry by the appropiate index of "sitePath"
                mutations[[site]][[targetIndex]] <- mutPath
                # (Re)assign the class and site attribute for the site
                class(mutations[[site]]) <- "sitePath"
                attr(mutations[[site]], "site") <- as.integer(site)
            }
        }
    }
    attr(mutations, "paths") <- paths
    attr(mutations, "reference") <- reference
    attr(mutations, "refSeqName") <- refSeqName
    class(mutations) <- "fixationSites"
    return(mutations)
}

#' @export
fixationSites <- function(paths, reference, gapChar, ...)
    UseMethod("fixationSites")

#' @export
print.fixationSites <- function(x, ...) {
    cat("Result for", length(attr(x, "paths")), "paths:\n\n")
    if (length(x) == 0) {
        cat("No fixation found\n")
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

compareMutPathAA <- function(e, s) {
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
#' @name multiFixationSites
#' @description
#' After finding the \code{\link{lineagePath}} of a phylogenetic tree,
#' \code{multiFixationSites} uses the result to find those sites that show
#' multiple fixations on some, if not all, of the lineages.
#' @param searchDepth
#' The function uses heuristic search but the termination of the search
#' cannot be intrinsically decided. \code{searchDepth} is needed to tell
#' the search when to stop.
#' @examples
#' data(h3n2_tree_reduced)
#' data(h3n2_align_reduced)
#' tree <- addMSA(h3n2_tree_reduced, alignment = h3n2_align_reduced)
#' multiFixationSites(lineagePath(tree))
#' @return
#' \code{multiFixationSites} returns sites with multiple fixations.
#' @export
multiFixationSites.lineagePath <- function(paths,
                                           reference = NULL,
                                           gapChar = '-',
                                           minEffectiveSize = NULL,
                                           searchDepth = 1,
                                           ...) {
    tree <- attr(paths, "tree")
    nTips <- length(tree$tip.label)
    align <- attr(paths, "align")
    # Generate the site mapping from reference
    refSeqName <- reference
    reference <- checkReference(tree, align, reference, gapChar)
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    # Get the "minEffectiveSize" for each fixation
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- nTips / length(unique(unlist(paths)))
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    minEffectiveSize <- ceiling(minEffectiveSize)
    # Get the "searchDepth" for heuristic search
    if (searchDepth < 1) {
        stop("\"searchDepth\" should be at least 1")
    } else {
        searchDepth <- ceiling(searchDepth)
    }
    divNodes <- unique(divergentNode(paths))
    paths <- extendPaths(paths, tree)
    # Get all the nodes that are not at divergent point
    nodes <- setdiff(unlist(paths), divNodes)
    # Get the sequence of the children tips that are descendant of the nodes.
    # Assign the tip index to the sequences for retrieving the tip name
    nodeAlign <- lapply(nodes, function(n) {
        childrenNode <- tree$edge[which(tree$edge[, 1] == n), 2]
        childrenNode <- setdiff(childrenNode, c(divNodes, nodes))
        tips <- ChildrenTips(tree, childrenNode)
        res <- align[tips]
        names(res) <- tips
        return(res)
    })
    # Assign the node names to the "nodeAlign" list
    names(nodeAlign) <- nodes
    # Store all the tips by node and their fixed AA to avoid repeating
    # calculation. Each entry in the list stores the info for a single
    # site. Under each site are groups of "nodeTips" where tip names
    # and the summary of their AAs at the site.
    nodeAAsum <- list()
    # "res" is going to be the return of this function. Each entry in
    # the list is the "sitePath" for a site. Each site ("sitePath")
    # consists of "mutPath" that is named by the starting node name.
    # The fixed AA and number of non-dominant AA is also stored.
    res <- list()
    # Iterate each path in the "lineagePath" object
    for (path in paths) {
        # Try every node as the terminal node
        for (maxLen in seq_along(path)[-1]) {
            # For the sequences after the terminal node,
            # examine them as a whole
            afterTips <-
                as.integer(ChildrenTips(tree, path[maxLen]))
            # Group the sequences by nodes in the path before
            pathBefore <-
                setdiff(path[seq_len(maxLen - 1)], divNodes)
            # Iterate every site
            for (i in loci) {
                s <- reference[i] - 1
                # "tableAA" is similar to R function "table"
                # Here the AA at site "s" for all tips is summarized
                afterSummary <- tableAA(align[afterTips], s)
                # TODO: "tolerance" is used to track total number of
                # non-dominant AA for the site initialized with the
                # number of non-dominant AA in the "afterTips".
                # May jump to next site if exceeding tolerance
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
                # "nodeTips" is vector of tips with an attribute
                # of "AA" to store fixed AA, and an attribute of
                # "nonDominant" to store non-dominant AA.
                nodeTips <- integer(0)
                # We need a "previousAA" and a "currentAA" to track
                # the AA along the "lineagePath". In the summary
                # stage, we only focus on the purely fixed AA.
                previousAA <- NULL
                currentAA <- NULL
                previousNode <- NULL
                # "nodeSummaries" groups tips by node in the summary
                # stage. In the case that the adjacent ndoes have
                # the same AA purely fixed, they will be grouped
                # into one.
                nodeSummaries <- list()
                # Iterate each single node in the "pathBefore". This is
                # the summary stage
                for (node in as.character(pathBefore)) {
                    # If the node has a record in nodeAAsum
                    nodeTips <- nodeAAsum[[site]][[node]]
                    # Summarize the node if not existed
                    if (is.null(nodeTips)) {
                        # Get the related descendant tips from "nodeAlign"
                        nodeTips <-
                            as.integer(names(nodeAlign[[node]]))
                        aaSummary <- tableAA(nodeAlign[[node]], s)
                        # Assign the "aaSummary" to the tip names
                        attr(nodeTips, "aaSummary") <- aaSummary
                        # Store the result to avoid repeating calculation
                        # TODO: Bug fix: 'nodeTips' will be allocated as
                        # a named vector if its length is of 1. An S4 class
                        # should be created for 'nodeAAsum' to solve this
                        # typing problem. Here a vector of c(1,2) is used
                        # to first guarantee a type of list.
                        nodeAAsum[[site]][[node]] <- c(1, 2)
                        nodeAAsum[[site]][[node]] <- nodeTips
                    }
                    # Extract "aaSummary" if existed
                    aaSummary <- attr(nodeTips, "aaSummary")
                    # If AA is purely fixed for the node
                    if (length(aaSummary) == 1) {
                        currentAA <- names(aaSummary)
                    } else {
                        currentAA <- NULL
                    }
                    # Attach the node to the "preivous" node if they're
                    # both purely fixed and the fixed AA is the same
                    if (!is.null(previousAA) &&
                        !is.null(currentAA) &&
                        previousAA == currentAA) {
                        node <- previousNode
                        nodeTips <-
                            c(nodeSummaries[[node]], nodeTips)
                        attr(nodeTips, "aaSummary") <-
                            attr(nodeSummaries[[node]], "aaSummary") +
                            aaSummary
                        # Oddly, R uses the name of the first element
                        # in the numeric vector when adding two named
                        # number. So there is no name (AA) assignment
                    }
                    # Assign or re-assign the nodeTips with "aaSummary"
                    # to the "nodeSummaries"
                    nodeSummaries[[node]] <- nodeTips
                    previousAA <- currentAA
                    previousNode <- node
                }
                # Attach the "afterTips" to "nodeSummaries"
                if (!is.null(previousAA) &&
                    !is.null(afterAA) &&
                    previousAA == afterAA) {
                    nodeTips <- c(nodeSummaries[[previousNode]], afterTips)
                    attr(nodeTips, "aaSummary") <-
                        attr(nodeSummaries[[previousNode]], "aaSummary") +
                        afterSummary
                    nodeSummaries[[previousNode]] <- nodeTips
                } else {
                    nodeSummaries[[as.character(path[maxLen])]] <- afterTips
                }
                # Skip to the next "site" if AA of "pathBefore" is
                # purely fixed or the terminal AA is the same as
                # fixed AA of "afterTips". This avoids repetition.
                if (length(nodeSummaries) <= 1) {
                    next
                } else if (site %in% names(res)) {
                    aTips <- nodeSummaries[[length(nodeSummaries)]]
                    exist <- vapply(
                        res[[site]],
                        FUN = function(ep) {
                            all(aTips %in% ep[[length(ep)]])
                        },
                        FUN.VALUE = logical(1)
                    )
                    if (any(exist)) {
                        next
                    }
                }
                seg <- minimizeEntropy(nodeSummaries,
                                       minEffectiveSize,
                                       searchDepth)
                targetIndex <- NULL
                if (length(seg) < 2) {
                    next
                } else if (!site %in% names(res)) {
                    targetIndex <- 1
                } else {
                    # Some site may have multiple fixation on multiple
                    # lineages. The following is for deciding at which
                    # index should it be assigned in the "res[[site]]"
                    # Retrieve the existing mutation path of the site
                    existPath <- res[[site]]
                    # Add new fixation for the site. Assume none of
                    # the existing mutation path has the same
                    # mutations as "seg"
                    targetIndex <- length(existPath) + 1
                    # Which mutaiton path has the same mutations as "seg"
                    existIndex <- which(
                        vapply(
                            existPath,
                            FUN = compareMutPathAA,
                            FUN.VALUE = logical(1),
                            seg
                        )
                    )
                    if (length(existIndex) > 0) {
                        # afterTips <- unlist(tail(seg, 1))
                        # "Adding state" for each existing mutation
                        # path that has the same mutations as "seg" does.
                        #
                        # 0L: "seg" is a subset of existing path and won't
                        # be added to the "res"
                        #
                        # 1L: The path is different from "seg". They just
                        # happen to share the same mutations. A new
                        # mutation path will be created
                        #
                        # 2L: The existing path is a subset of "seg" and
                        # needs to be replaced by "seg"
                        qualified <- vapply(
                            existPath[existIndex],
                            FUN = function(ep) {
                                existTips <- unlist(tail(ep, 1))
                                if (all(existTips %in% afterTips)) {
                                    return(2L)
                                }
                                if (!all(afterTips %in% existTips)) {
                                    return(1L)
                                }
                                return(0L)
                            },
                            FUN.VALUE = integer(1)
                        )
                        r <- qualified == 2L
                        if (any(r)) {
                            # Remove existing paths with an "adding state"
                            # of 2L and add the "seg"
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
                # Assign the result to the "res[[site]]"
                res[[site]][[targetIndex]] <- seg
                # Assign or re-assign class and "site" to
                # the "res[[site]]"
                class(res[[site]]) <- "sitePath"
                attr(res[[site]], "site") <- i
            }
        }
    }
    attr(res, "paths") <- paths
    attr(res, "reference") <- reference
    attr(res, "refSeqName") <- refSeqName
    class(res) <- "multiFixationSites"
    return(res)
}

#' @export
multiFixationSites <- function(paths, reference, gapChar, ...)
    UseMethod("multiFixationSites")

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
