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
#' data("zikv_tree")
#' data("zikv_align")
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
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
    for (m in x) {
        mutName <- character(0)
        for (tips in m) {
            aa <- attr(tips, "AA")
            mutName <-
                c(mutName, paste0(aa, "(", length(tips), ")"))
        }
        cat(paste0(mutName, collapse = " -> "), "\n")
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
#' it's the maximum for both ancestral and descendant.
#' @param minEffectiveSize
#' A vector of two integers to specifiy minimum tree tips involved
#' before/after mutation. Otherwise the mutation will not be counted into
#' the return. If more than one number is given, the ancestral takes the first
#' and descendant takes the second as the minimum. If only given one number,
#' it's the minimum for both ancestral and descendant.
#' @param extendedSearch
#' Whether to extend the search. The terminal of each \code{lineagePath} is
#' a cluster of tips. To look for the fixation mutation in the cluster,
#' the common ancestral node of farthest tips (at least two) will be
#' the new terminal search point.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' fixationSites(
#'     lineagePath(tree, 0.996),
#'     tolerance = c(1, 1),
#'     minEffectiveSize = c(10, 10)
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
                                      tolerance = 0,
                                      minEffectiveSize = NULL,
                                      extendedSearch = TRUE,
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
    } else if (any(minEffectiveSize) <= 0) {
        stop("\"minEffectiveSize\" can only be positive number")
    } else {
        minAnc <- minEffectiveSize
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
    if (extendedSearch) {
        paths <- extendPaths(paths, tree)
    }
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
                from <- tree$tip.label[beforeTips]
                to <- tree$tip.label[afterTips]
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
                    targetIndex <- which(vapply(
                        existPath,
                        FUN = function(i) {
                            paste(names(i), collapse = "") == mut
                        },
                        FUN.VALUE = logical(1)
                    ))
                    if (length(targetIndex) == 0) {
                        # Add new fixation for the site if it's not existed
                        targetIndex <- 1
                    } else {
                        # Retrieve existing "sitePath" entry to compare
                        # with the current "sitePath" enry
                        existMut <- existPath[[targetIndex]]
                        # A new "sitePath" will be added to the site
                        # or an old "sitePath" will be replaced if
                        # it's actually a subset of the new one.
                        if (all(!to %in% existMut[[a]])) {
                            targetIndex <- length(existPath) + 1
                        } else if (length(to) <= length(existMut[[a]])) {
                            next
                        }
                    }
                }
                # Construct a "mutPath" for each fixation
                attr(from, "AA") <- b
                attr(to, "AA") <- a
                mutPath <- list(from, to)
                names(mutPath) <- c(b, a)
                # Store the "mutPath" to the appropiate index
                mutations[[site]][[targetIndex]] <- mutPath
                # (Re)assign the class and site attribute for
                # the site
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
fixationSites <- function(paths,
                          reference,
                          gapChar,
                          tolerance,
                          minEffectiveSize,
                          extendedSearch,
                          ...)
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

#' @rdname findSites
#' @name multiFixationSites
#' @description
#' After finding the \code{\link{lineagePath}} of a phylogenetic tree,
#' \code{multiFixationSites} uses the result to find those sites that show
#' multiple fixations on some, if not all, of the lineages.
#' @return
#' \code{multiFixationSites} returns sites with multiple fixations.
#' @export
multiFixationSites.lineagePath <- function(paths,
                                           reference = NULL,
                                           gapChar = '-',
                                           tolerance = 0,
                                           minEffectiveSize = NULL,
                                           extendedSearch = TRUE,
                                           ...) {
    nTips <- length(tree$tip.label)
    tree <- attr(paths, "tree")
    align <- attr(paths, "align")
    refSeqName <- reference
    reference <- checkReference(tree, align, reference, gapChar)
    if (!is.numeric(tolerance)) {
        stop("\"tolerance\" only accepts numeric")
    } else {
        tolerance <- if (tolerance < 0) {
            stop("\"tolerance\" can only be positive number")
        } else if (tolerance < 0.5) {
            tolerance * nTips
        } else if (tolerance < 1) {
            (1 - tolerance) * nTips
        } else {
            as.integer(tolerance)
        }
    }
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- nTips / (length(paths) * 10)
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    divNodes <- unique(divergentNode(paths))
    if (extendedSearch) {
        paths <- extendPaths(paths, tree)
    }
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
    # Assign the node index to the "nodeAlign" list
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
        for (maxLen in seq_along(path)[-1]) {
            # For the sequences after the node, examine them as a whole
            afterTips <- ChildrenTips(tree, path[maxLen])
            if (length(afterTips) < minEffectiveSize) {
                next
            }
            after <- align[afterTips]
            # Group the sequences by nodes in the path before
            pathBefore <-
                setdiff(path[seq_len(maxLen - 1)], divNodes)
            # Iterate every site
            for (i in seq_along(reference)) {
                s <- reference[i] - 1
                a <- summarizeAA(after, s, tolerance)
                if (is.na(a)) {
                    # The AA of the tips after are not fixed
                    # so the site is disqualifed and go for
                    # the next site
                    next
                }
                site <- as.character(i)
                # "toleranceSum" is used to track total number of
                # non-dominant AA for the site initialized with the
                # number of non-dominant AA in the "afterTips"
                toleranceSum <- attr(a, "n")
                # "nodeTips" is vector of tips with an attribute
                # of "AA" to store fixed AA, and an attribute of
                # "nonDominant" to store non-dominant AA.
                nodeTips <- integer(0)
                # We need a "previousAA" and a "currentAA" to track
                # the AA along the "lineagePath". In the summary
                # stage, we only focus on the purely fixed AA.
                previousAA <- NULL
                currentAA <- NULL
                # "nodeSummaries" groups tips by node in the summary
                # stage. In the case that the adjacent ndoes have
                # the same AA purely fixed, they will be grouped
                # into one.
                nodeSummaries <- list()
                # Iterate all nodes in the "pathBefore". This is
                # the summary stage
                for (node in pathBefore) {
                    # If the node has a record in nodeAAsum
                    nodeTips <- nodeAAsum[[site]][[node]]
                    if (is.null(nodeTips)) {
                        # Get the related descendant tips
                        nodeTips <- names(nodeAlign[[node]])
                        # "tableAA" is similar to R function "table"
                        aaSummary <- tableAA(nodeAlign[[node]], s)
                        attr(nodeTips, "aaSummary") <- aaSummary
                        # Store the result to avoid repeating calculation
                        nodeAAsum[[site]][[node]] <- nodeTips
                    }
                    aaSummary <- attr(nodeTips, "aaSummary")
                    if (length(aaSummary) == 1) {
                        currentAA <- names(aaSummary)
                    } else {
                        currentAA <- NULL
                    }
                    if (!is.null(previousAA) &&
                        previousAA == currentAA) {
                        nodeTips <- c(nodeSummaries[[node]], nodeTips)
                        attr(nodeTips, "aaSummary") <-
                            attr(nodeSummaries[[node]], "aaSummary") +
                            aaSummary
                    }
                    nodeSummaries[[node]] <- nodeTips
                    previousAA <- currentAA
                }
                if (length(nodeSummaries) == 1 || currentAA == a) {
                    next
                }
                # TODO: integrate the "nodeSummaries"
                # The function should output tips in groups
                # where AA is fixed (mostly pure)
                if (is.null(res[[site]])) {
                    targetIndex <- 1
                } else {
                    targetIndex <- length(res[[site]]) + 1
                }
                class(res[[site]]) <- "sitePath"
                attr(res[[site]], "site") <- as.integer(site)
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
multiFixationSites <- function(paths,
                               reference = NULL,
                               gapChar = '-',
                               tolerance = 0,
                               minEffectiveSize = NULL,
                               extendedSearch = TRUE,
                               ...)
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
