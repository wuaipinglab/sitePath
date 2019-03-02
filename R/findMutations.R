#' @rdname pre-assessment
#' @name pre-assessment
#' @title Things can be done before the analysis
#' @description
#' \code{similarityMatrix} calculates similarity between aligned sequences
#' The similarity matrix can be used in \code{\link{groupTips}}
#' or \code{\link{sitePath}}
#' @param tree The return from \code{\link{addMSA}} function
#' @examples
#' data("zikv_tree")
#' data("zikv_align")
#' tree <- addMSA(zikv_tree, seqs = zikv_align)
#' simMatrix <- similarityMatrix(tree)
#' @return
#' \code{similarityMatrix} returns a diagonal matrix of
#' similarity between sequences
#' @importFrom methods is
#' @export
similarityMatrix <- function(tree) {
    if (!inherits(tree, "phylo")) {
        stop("\"tree\" is not class phylo")
    } else if (is.null(attr(tree, "alignment"))) {
        stop("No alignment found in \"tree\"")
    }
    sim <- getSimilarityMatrix(attr(tree, "alignment"))
    dimnames(sim) <- list(tree$tip.label, tree$tip.label)
    return(sim)
}

#' @rdname pre-assessment
#' @description
#' \code{sneakPeek} is intended to plot similarity as a threshold
#' against number of output sitePath. This plot is intended to give user
#' a feel about how many sitePaths they should expect from
#' the similarity threshold. The number of sitePath should not
#' be too many or too few. The result excludes where the number of sitePath
#' is greater than number of tips divided by 20 or self-defined maxPath.
#' The zero sitePath result will also be excluded
#' @param step
#' the similarity window for calculating and ploting. To better
#' see the impact of threshold on path number. This is preferably
#' specified. The default is one 50th of the difference between 1
#' and minimal pairwise sequence similarity.
#' @param maxPath
#' maximum number of path to return show in the plot. The number of path
#' in the raw tree can be far greater than trimmed tree. To better
#' see the impact of threshold on path number. This is preferably
#' specified. The default is one 20th of tree tip number.
#' @param minPath
#' minimum number of path to return show in the plot. To better
#' see the impact of threshold on path number. This is preferably
#' specified. The default is 1.
#' @param makePlot whether make a dot plot when return
#' @examples
#' sneakPeek(tree)
#' @return
#' \code{sneakPeek} return the similarity threhold against number of sitePath.
#' There will be a simple dot plot between threshold and path number if
#' \code{makePlot} is TRUE.
#' @importFrom methods is
#' @importFrom graphics plot
#' @export
sneakPeek <- function(tree,
                      step = NULL,
                      maxPath = NULL,
                      minPath = 1,
                      makePlot = TRUE) {
    simMatrix <- similarityMatrix(tree)
    minSim <- min(simMatrix)
    if (is.null(step)) {
        step <- round(minSim - 1, 3) / 50
    }
    if (is.null(maxPath)) {
        maxPath <- length(tree$tip.label) / 20
    } else if (maxPath <= 0) {
        stop("Invalid \"maxPath\": less than or equal to zero")
    }
    if (minPath >= maxPath) {
        stop("Invalid \"minPath\": greater than \"maxPath\"")
    } else if (minPath < 0) {
        stop("Invalid \"minPath\": less than zero")
    }
    similarity <- numeric(0)
    pathNum <- integer(0)
    for (s in seq(1, minSim, step)) {
        paths <- sitePath(
            tree,
            similarity = s,
            simMatrix = simMatrix,
            forbidTrivial = FALSE
        )
        if (maxPath < length(paths)) {
            next
        } else if (length(paths) <= minPath) {
            break
        }
        similarity <- c(similarity, s)
        pathNum <- c(pathNum, length(paths))
    }
    if (makePlot) {
        plot(similarity, pathNum)
    }
    return(data.frame(similarity, pathNum))
}

#' @rdname findSites
#' @name SNPsites
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
#' name of reference for site numbering. The name has to be one of the
#' sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar the character to indicate gap.
#' The numbering will skip the gapChar if reference sequence if specified.
#' @param minSNP minimum number of amino acid variation to be a SNP
#' @examples
#' data("zikv_tree")
#' data("zikv_align")
#' tree <- addMSA(zikv_tree, seqs = zikv_align)
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
            FUN.VALUE = "",
            reference[i]
        ))
        if (sum(SNP > minSNP) >= 2) {
            qualified <- c(qualified, i)
        }
    }
    return(qualified)
}

#' @rdname findSites
#' @name fixationSites
#' @description
#' After finding the \code{\link{sitePath}} of a phylogenetic tree, 
#' \code{fixationSites} uses the result to find those sites that show 
#' fixation on some, if not all, of the lineages. Parallel evolution is
#' relatively common in RNA virus. There is chance that some site be fixed
#' in one lineage but does not show fixation because of different 
#' sequence context.
#' @param paths
#' a \code{sitePath} object returned from \code{\link{sitePath}} function
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
#' Whether to extend the search. The terminal of each \code{sitePath} is
#' a cluster of tips. To look for the fixation mutation in the cluster,
#' the common ancestral node of farthest tips (at least two) will be
#' the new terminal search point.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' fixationSites(
#'     sitePath(tree, 0.996),
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
fixationSites.sitePath <- function(paths,
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
    mutations <- list()
    for (minLen in 2:max(lengths(paths))) {
        # literate all sitePath at the same time
        for (path in unique(ancestralPaths(paths, minLen))) {
            # first, the number of tips should meet the minimum size
            # otherwise we'll go for the next node
            afterTips <- ChildrenTips(tree, tail(path, 1))
            if (length(afterTips) < minDesc) {
                next
            }
            pathBefore <- path[seq_len(length(path) - 1)]
            excludedTips <- lapply(
                pathBefore[which(pathBefore %in% divNodes)],
                FUN = function(node) {
                    children <- tree$edge[which(tree$edge[, 1] == node), 2]
                    children <-
                        children[which(children > length(tree$tip.label) &
                                           !children %in% path)]
                    return(ChildrenTips(tree, children))
                }
            )
            beforeTips <-
                which(!seq_along(tree$tip.label) %in% c(afterTips, unlist(excludedTips)))
            if (length(beforeTips) < minAnc &&
                length(excludedTips) == 0) {
                next
            }
            # get the sequences for ancestral and descendant group
            after <- align[afterTips]
            before <- align[beforeTips]
            # iterate through each site and compare the two groups
            for (i in seq_along(reference)) {
                s <- reference[i] - 1
                # get the dominant AA of the site for the two groups
                b <- summarizeAA(before, s, toleranceAnc)
                a <- summarizeAA(after, s, toleranceDesc)
                if (!(is.na(b) || is.na(a)) && a != b) {
                    mut <- paste(b, i, a, sep = "")
                    from <- tree$tip.label[beforeTips]
                    to <- tree$tip.label[afterTips]
                    # add the mutation or update the tips for a mutation
                    # if it's in the mutaiton list.
                    # NOTE: Sometimes, the involved tips could be incompletely
                    # retrived. The replacement assumes the mutation involves
                    # a same set of tips but this might not always be true.
                    if (!mut %in% names(mutations) ||
                        length(to) > length(mutations[[mut]]$to)) {
                        mutations[[mut]] <- list(ancestral = from,
                                                 descendant = to)
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
fixationSites <- function(paths,
                          reference,
                          gapChar,
                          tolerance,
                          minEffectiveSize,
                          extendedSearch,
                          ...)
    UseMethod("fixationSites")

checkReference <- function(tree, align, reference, gapChar) {
    if (is.null(reference)) {
        reference <- seq_len(nchar(align[1]))
    } else {
        if (!is.character(gapChar) ||
            nchar(gapChar) != 1 || length(gapChar) != 1) {
            stop("\"gapChar\" only accepts one single character")
        }
        reference <-
            getReference(align[which(tree$tip.label == reference)], gapChar)
    }
    return(reference)
}

#' @export
as.data.frame.fixationSites <- function(x, ...) {
    tree <- attr(x, "tree")
    reference <- attr(x, "reference")
    res <- as.data.frame(matrix(
        nrow = length(tree$tip.label),
        ncol = length(reference),
        dimnames = list(tree$tip.label, seq_along(reference))
    ))
    for (m in names(x)) {
        site <- as.integer(substr(m, 2, nchar(m) - 1))
        fixedAA <- substr(m, nchar(m), nchar(m))
        for (desc in x[[m]][[2]]) {
            res[desc, site] <- AA_FULL_NAMES[tolower(fixedAA)]
        }
    }
    whichNA <- is.na(res)
    res <-
        res[rowSums(whichNA) < ncol(res), colSums(whichNA) < nrow(res)]
    return(res)
}

#' @export
print.fixationSites <- function(x, ...) {
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
#' @name fixationSites
#' @description
#' After finding the \code{\link{sitePath}} of a phylogenetic tree, 
#' \code{multiFixationSites} uses the result to find those sites that show 
#' multiple fixations on some, if not all, of the lineages.
#' @return 
#' \code{multiFixationSites} returns sites with multiple fixations.
#' @export
multiFixationSites.sitePath <- function(paths,
                                        reference = NULL,
                                        gapChar = '-',
                                        tolerance = 0,
                                        minEffectiveSize = NULL,
                                        extendedSearch = TRUE,
                                        ...) {
    tree <- attr(paths, "tree")
    align <- attr(paths, "align")
    # refSeqName <- reference
    reference <-
        sitePath:::checkReference(tree, align, reference, gapChar)
    if (!is.numeric(tolerance)) {
        stop("\"tolerance\" only accepts numeric")
    } else {
        tolerance <- if (tolerance < 0) {
            stop("\"tolerance\" can only be positive number")
        } else if (tolerance < 0.5) {
            tolerance * length(tree$tip.label)
        } else if (tolerance < 1) {
            (1 - tolerance) * length(tree$tip.label)
        } else {
            as.integer(tolerance)
        }
    }
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- 10
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    divNodes <- unique(sitePath:::divergentNode(paths))
    if (extendedSearch) {
        paths <- sitePath:::extendPaths(paths, tree)
    }
    # get all the nodes that are not at divergent point
    nodes <- setdiff(unlist(paths), divNodes)
    # get the sequence of the children tips that are descendant of the nodes.
    # assign the tip index to the sequences for retrieving the tip name
    nodeAlign <- lapply(nodes, function(n) {
        childrenNode <- tree$edge[which(tree$edge[, 1] == n), 2]
        childrenNode <- setdiff(childrenNode, c(divNodes, nodes))
        tips <- sitePath:::ChildrenTips(tree, childrenNode)
        res <- align[tips]
        names(res) <- tips
        return(res)
    })
    # assign the node index to the 'nodeAlign' list
    names(nodeAlign) <- nodes
    # store all the tips by node and the fixed AA
    # to avoid repeating calculation
    nodeAAsum <- list()
    res <- list()
    for (minLen in 2:max(lengths(paths))) {
        # literate all sitePath at the same time
        for (path in unique(sitePath:::ancestralPaths(paths, minLen))) {
            # for the sequences after the node, examine them as a whole
            afterTips <-
                sitePath:::ChildrenTips(tree, tail(path, 1))
            if (length(afterTips) < minEffectiveSize) {
                next
            }
            after <- align[afterTips]
            # group the sequences before by nodes
            pathBefore <- path[seq_len(length(path) - 1)]
            nodeSeqsBefore <-
                nodeAlign[as.character(setdiff(pathBefore, divNodes))]
            # iterate every site
            for (i in seq_along(reference)) {
                s <- reference[i] - 1
                a <- sitePath:::summarizeAA(after, s, tolerance)
                if (is.na(a)) {
                    # the AA of the tips after are not fixed
                    # so the site is disqualifed and go for
                    # the next site
                    next
                }
                # total number of non-dominant AA for the site
                # initialized with the number of non-dominant AA
                # in the 'afterTips'
                toleranceSum <- attr(a, "n")
                # 'fixationNodes' groups tips by node. It's a growing
                # list or NULL if disqualified. If qualified it's
                # going to be included in the 'res'
                fixationNodes <- list()
                # 'fixationNodes' is vector of tips with
                # an attribute of 'n' to store fixed AA
                nodeTips <- integer(0)
                # 'aaSum' collects fixed AA for the site
                aaSum <- character(0)
                # iterate all nodes in the 'pathBefore' and check
                # if their tips have fixed AA for the site
                for (node in names(nodeSeqsBefore)) {
                    # if the node has a record in nodeAAsum
                    nodeTips <- nodeAAsum[[as.character(i)]][[node]]
                    if (is.null(nodeTips)) {
                        # summarize the AA for the tips in the node.
                        # The output is either an AA (with the number
                        # of non-dominant AA as an attribute) or NA
                        nodeAA <- sitePath:::summarizeAA(
                            seqs = nodeAlign[[node]],
                            siteIndex = s,
                            tolerance = tolerance
                        )
                        if (is.na(nodeAA)) {
                            # mark the node with non-fixed site as NA
                            nodeAAsum[[as.character(i)]][[node]] <-
                                integer(0)
                            # the rest of the node in the path will be
                            # ignored once a node appears to be
                            # 'un-fixed'. And the site is disqualified
                            fixationNodes <- NULL
                            break
                        } else {
                            toleranceSum <- toleranceSum + attr(nodeAA, "n")
                            if (toleranceSum > tolerance) {
                                # mark the node with as NA because the
                                # the accumlated number of non-dominant
                                # AA excceeds total tolerance
                                nodeAAsum[[as.character(i)]][[node]] <-
                                    integer(0)
                                # the rest of the node in the path will be
                                # ignored. And the site
                                # will be disqualified
                                fixationNodes <- NULL
                                break
                            } else {
                                # get the descendant tips for the node
                                nodeTips <-
                                    as.integer(names(nodeAlign[[node]]))
                                # AA for the site, for the node
                                attr(nodeTips, "nonDominant") <-
                                    attr(nodeAA, "n")
                                attr(nodeAA, "n") <- NULL
                                attr(nodeTips, "AA") <- nodeAA
                                # for summarizing AA for the site
                                aaSum <- c(aaSum, nodeAA)
                                nodeTips <- list(nodeTips)
                                names(nodeTips) <- node
                                # store the node with fixed AA
                                nodeAAsum[[as.character(i)]] <-
                                    c(nodeAAsum[[as.character(i)]], nodeTips)
                                # collect tips with fixed AA by nodes
                                # and keep looking for the next node
                                fixationNodes <-
                                    c(fixationNodes, nodeTips)
                            }
                        }
                    } else if (length(nodeTips) == 0) {
                        fixationNodes <- NULL
                        # the site will be disqualified if the node
                        # is marked AA for some reason
                        break
                    } else {
                        toleranceSum <-
                            toleranceSum + attr(nodeTips, "nonDominant")
                        if (toleranceSum > tolerance) {
                            nodeAAsum[[as.character(i)]][[node]] <-
                                integer(0)
                            fixationNodes <- NULL
                            break
                        } else {
                            nodeTips <- list(nodeTips)
                            names(nodeTips) <- node
                            fixationNodes <-
                                c(fixationNodes, nodeTips)
                        }
                    }
                }
                # store the nodes with fixed AA for the site
                if (!is.null(fixationNodes) &&
                    attr(nodeTips[[node]], "AA") != a &&
                    length(unique(aaSum)) > 1) {
                    qualified <- TRUE
                    previousAA <- attr(fixationNodes[[1]], "AA")
                    currentAAnum <- length(fixationNodes[[1]])
                    for (node in fixationNodes[-1]) {
                        currentAA <- attr(node, "AA")
                        if (currentAA == previousAA) {
                            currentAAnum <- currentAAnum + length(node)
                        } else {
                            if (currentAAnum < minEffectiveSize) {
                                qualified <- FALSE
                                break
                            }
                            previousAA <- currentAA
                            currentAAnum <- length(node)
                        }
                    }
                    if (qualified && currentAAnum >= minEffectiveSize) {
                        attr(afterTips, "nonDominant") <- attr(a, "n")
                        attr(a, "n") <- NULL
                        attr(afterTips, "AA") <- a
                        terminalTips <- list(afterTips)
                        names(terminalTips) <- tail(path, 1)
                        res[[as.character(i)]] <-
                            c(fixationNodes, terminalTips)
                    }
                }
            }
        }
    }
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
