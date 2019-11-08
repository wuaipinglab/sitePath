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
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' SNPsites(tree)
#' @return \code{SNPsite} returns a list of qualified SNP site
#' @export
SNPsites <- function(tree, minSNP = NULL) {
    if (is.null(minSNP)) {
        minSNP <- length(tree[["tip.label"]]) / 10
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
        SNP <- table(vapply(
            X = alignedSeq,
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
    # A 'sitePath' composes of all the fixation paths for a single site.
    #  So each 'm' represent a single fixation path
    for (m in x) {
        if (length(m) == 2) {
            mutName <-
                paste0(attr(m[[1]], "AA"), attr(x, "site"), attr(m[[2]], "AA"))
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

.tipSeqsAlongPathNodes <- function(paths, divNodes, tree, align) {
    allNodes <- unlist(paths)
    terminalNodes <- vapply(
        X = paths,
        FUN = function(p) {
            p[length(p)]
        },
        FUN.VALUE = integer(1)
    )
    # Get all the nodes that are not at divergent point
    nodes <- setdiff(allNodes, divNodes)
    # Get the sequence of the children tips that are descendant of
    # the nodes. Assign the tip index to the sequences for
    # retrieving the tip name
    nodeAlign <- lapply(nodes, function(n) {
        isTerminal <- FALSE
        if (n %in% terminalNodes) {
            childrenNode <- n
            isTerminal <- TRUE
        } else {
            childrenNode <- tree$edge[which(tree$edge[, 1] == n), 2]
            # Keep the node that is not on the path.
            childrenNode <- setdiff(childrenNode, allNodes)
        }
        tips <- .childrenTips(tree, childrenNode)
        res <- align[tips]
        attr(res, "isTerminal") <- isTerminal
        names(res) <- tips
        return(res)
    })
    # Assign the node names to the 'nodeAlign' list
    names(nodeAlign) <- nodes
    return(nodeAlign)
}

.findFixationSite <- function(paths,
                              tree,
                              align,
                              nodeAlign,
                              divNodes,
                              reference,
                              minimizeEntropy,
                              minEffectiveSize,
                              searchDepth) {
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    # The variable to store the result from entropy minimization for
    # each path with those purely fixed excluded.
    res <- list()
    # Iterate each path
    rootNode <- attr(paths, "rootNode")
    if (!rootNode %in% names(nodeAlign)) {
        paths <- lapply(paths, function(p)
            p[-1])
    }
    for (path in paths) {
        # Iterate every loci (variant sites)
        for (i in loci) {
            site <- as.character(i)
            # The index to use for cpp code
            s <- reference[i] - 1
            # Assign a variable to store the tip names and their info on
            # amino acids. They are the potential fixation segment
            nodeTips <- integer()
            previousAA <- NULL
            currentAA <- NULL
            previousNode <- NULL
            # The input for entropy minimization calculation
            nodeSummaries <- list()
            # Divergent nodes are not included anywhere in the result
            for (node in as.character(setdiff(path, divNodes))) {
                # Get the related descendant tips and related sequences
                nodeTips <- as.integer(names(nodeAlign[[node]]))
                # Frequency of the amino acids at the site
                aaSummary <- tableAA(nodeAlign[[node]], s)
                # Assoicate the amino acid frequence with the tip names
                attr(nodeTips, "aaSummary") <- aaSummary
                # Decide the current fixed amino acid
                if (length(aaSummary) == 1) {
                    currentAA <- names(aaSummary)
                } else {
                    currentAA <- NULL
                }
                # Attach the node to the preivous node if they're
                # both purely fixed and have the same AA fixed.
                if (!is.null(previousAA) &&
                    !is.null(currentAA) &&
                    previousAA == currentAA) {
                    node <- previousNode
                    # Combine the tips in the previous node
                    nodeTips <- c(nodeSummaries[[node]], nodeTips)
                    # Add up the amino acid frequency
                    attr(nodeTips, "aaSummary") <-
                        attr(nodeSummaries[[node]], "aaSummary") +
                        aaSummary
                    # Oddly, R uses the name of the first variable when
                    # adding two numeric vectors. So there is no need for
                    # names (AA) assignment
                }
                # Assign or re-assign the nodeTips with 'aaSummary'
                # to the 'nodeSummaries'
                nodeSummaries[[node]] <- nodeTips
                previousAA <- currentAA
                previousNode <- node
            }
            # Skip to the next locus if AA is fixed along the whole path
            if (length(nodeSummaries) >= 2) {
                seg <- minimizeEntropy(nodeSummaries,
                                       minEffectiveSize,
                                       searchDepth)
                if (length(seg) >= 2) {
                    targetIndex <- length(res[[site]]) + 1
                    attr(seg, "path") <- path
                    res[[site]][[targetIndex]] <- seg
                    attr(res[[site]], "site") <- i
                }
            }
        }
    }
    return(res)
}

.combineFixations <- function(fixations, tree, align) {
    # 'res' is going to be the return of this function. Each entry in
    # the list is the 'sitePath' for a site. Each site ('sitePath')
    # consists of 'mutPath' that is named by the starting node name.
    # The fixed AA and number of non-dominant AA is also stored.
    res <- list()
    for (segs in fixations) {
        i <- attr(segs, "site")
        site <- as.character(i)
        res[[site]][[1]] <- lapply(segs[[1]], function(tips) {
            attr(tips, "tipsAA") <- substr(x = align[tips],
                                           start = i,
                                           stop = i)
            return(tips)
        })
        attr(res[[site]], "site") <- i
        attr(res[[site]], "tree") <- tree
        class(res[[site]]) <- "sitePath"
        for (seg in segs[-1]) {
            # Assume a new fixation path is to add.
            targetIndex <- length(res[[site]]) + 1
            # The index to extract the terminal tips of the fixation.
            endIndex <- length(seg)
            finalAA <- attr(seg[[endIndex]], "AA")
            # The following is to decide if any fixation path can be
            # combined.
            existPath <- res[[site]]
            # The fixation before the temrinal tips should be identical
            # and the final fixed amino acid should be the same.
            toCombine <- vapply(
                X = existPath,
                FUN = function(ep) {
                    identical(lapply(seg, c)[-endIndex],
                              lapply(ep, c)[-endIndex]) &&
                        finalAA == attr(ep[[endIndex]], "AA")
                },
                FUN.VALUE = logical(1)
            )
            if (any(toCombine)) {
                existIndex <- which(toCombine)
                # These are the candidates to combine. The additional
                # condition be all the descendant tips are included.
                toCombine <- unique(unlist(lapply(
                    X = c(res[[site]][existIndex], list(seg)),
                    FUN = "[[",
                    ... = endIndex
                )))
                allTips <-
                    .childrenTips(tree, getMRCA(tree, toCombine))
                if (all(allTips %in% toCombine)) {
                    seg[[endIndex]] <- toCombine
                    attr(seg[[endIndex]], "AA") <- finalAA
                    res[[site]] <- res[[site]][-existIndex]
                    targetIndex <- length(res[[site]]) + 1
                } else {
                    # Add new fixation for the site if no existing
                    # mutation path can be combined with
                    targetIndex <- length(existPath) + 1
                }
            }
            seg <- lapply(seg, function(tips) {
                attr(tips, "tipsAA") <- substr(x = align[tips],
                                               start = i,
                                               stop = i)
                return(tips)
            })
            res[[site]][[targetIndex]] <- seg
            attr(res[[site]], "site") <- i
            attr(res[[site]], "tree") <- tree
            class(res[[site]]) <- "sitePath"
        }
    }
    return(res)
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
#' @param method
#' The strategy for predicting the fixation. The basic approach is entropy
#' minimization and can be achieved by adding or removing fixation point,
#' or by comparing the two.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' fixationSites(lineagePath(tree))
#' @return
#' \code{fixationSites} returns a list of fixation mutations
#' with names of the tips involved.
#' @importFrom utils tail
#' @export
fixationSites.lineagePath <- function(paths,
                                      minEffectiveSize = NULL,
                                      searchDepth = 1,
                                      method = c("compare", "insert", "delete"),
                                      ...) {
    tree <- attr(paths, "tree")
    nTips <- length(tree[["tip.label"]])
    align <- attr(paths, "align")
    # Generate the site mapping from reference
    reference <- attr(paths, "reference")
    # Decide which miniminzing strategy
    minimizeEntropy <- switch(
        match.arg(method),
        "compare" = minEntropyByComparing,
        "insert" = minEntropyByInserting,
        "delete" = minEntropyByDeleting
    )
    # Get the 'minEffectiveSize' for each fixation
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- nTips / length(unique(unlist(paths)))
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
    divNodes <- divergentNode(paths)
    paths <- .extendPaths(paths, tree)
    nodeAlign <- .tipSeqsAlongPathNodes(
        paths = paths,
        divNodes = divNodes,
        tree = tree,
        align = align
    )
    fixations <- .findFixationSite(
        paths = paths,
        tree = tree,
        align = align,
        nodeAlign = nodeAlign,
        divNodes = divNodes,
        reference = reference,
        minimizeEntropy = minimizeEntropy,
        minEffectiveSize = minEffectiveSize,
        searchDepth = searchDepth
    )
    res <- .combineFixations(fixations, tree, align)
    attr(res, "paths") <- paths
    attr(res, "reference") <- reference
    class(res) <- "fixationSites"
    return(res)
}

#' @export
fixationSites <- function(paths, ...)
    UseMethod("fixationSites")

#' @export
print.fixationSites <- function(x, ...) {
    cat("Result for", length(attr(x, "paths")), "paths:\n\n")
    if (length(x) == 0) {
        cat("No multi-fixation found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(attr(x, "reference"), "refSeqName")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}

.originalPath2sampled <- function(nodeAlign,
                                  original2sampled,
                                  sampledTree) {
    rootNode <- getMRCA(sampledTree, sampledTree[["tip.label"]])
    sampledEdge <- sampledTree[["edge"]]
    # To translate the path of original tree, first we need to find the
    # ancestral node on the sampled tree of the remaining tips in each
    # node of the original path.
    sampledNodeAlign <- list()
    for (nd in nodeAlign) {
        isTerminal <- attr(nd, "isTerminal")
        # Keep the sampled tips/sequences for each node
        nd <-
            nd[which(as.integer(names(nd)) %in% original2sampled)]
        # Remove the node if all the tips are gone during the sampling
        if (length(nd) == 0) {
            next
        }
        # Translate the tips in the original tree to the sampled tree
        tips <- match(as.integer(names(nd)), original2sampled)
        names(nd) <- as.character(tips)
        if (length(tips) == 1) {
            # A lineage path can disappear due to the sampling
            sampledNode <-
                sampledEdge[which(sampledEdge[, 2] == tips), 1]
        } else {
            sampledNode <- getMRCA(sampledTree, tips)
            if (!isTerminal && sampledNode != rootNode) {
                sampledNode <-
                    sampledEdge[which(sampledEdge[, 2] == sampledNode), 1]
            }
        }
        sampledNode <- as.character(sampledNode)
        # Combine the tips with the same ancestral node
        sampledNodeAlign[[sampledNode]] <-
            c(sampledNodeAlign[[sampledNode]], nd)
    }
    sampledNodes <- as.integer(names(sampledNodeAlign))
    # Node path on the sampled tree for each translated node.
    # Tranlate the paths on to the sampled tree
    sampledPaths <- list(nodepath(
        phy = sampledTree,
        from = rootNode,
        to = sampledNodes[1]
    ))
    for (n in sampledNodes[-1]) {
        nodePath <- nodepath(sampledTree, from = rootNode, to = n)
        # Descide if the current node path is a sub-path of existing paths
        # or any existing path is a sub-path to it.
        #
        # 0L: it's a new node path to add.
        #
        # 1L: the current node path is a sub-path and won't be added.
        #
        # 2L: one existing path is a sub-path and needs replacing
        qualified <- vapply(sampledPaths, function(p) {
            if (all(nodePath %in% p)) {
                return(1L)
            } else if (all(p %in% nodePath)) {
                return(2L)
            }
            return(0L)
        }, FUN.VALUE = integer(1))
        r <- qualified == 2L
        if (any(r)) {
            sampledPaths <- sampledPaths[-which(r)]
            sampledPaths <- c(sampledPaths, list(nodePath))
        } else if (all(qualified == 0L)) {
            sampledPaths <- c(sampledPaths, list(nodePath))
        }
    }
    attr(sampledPaths, "nodeAlign") <- sampledNodeAlign
    attr(sampledPaths, "rootNode") <- rootNode
    return(sampledPaths)
}

.sampleSummarize <- function(allMutations, nodeAlign, tree) {
    allSampledTips <- lapply(allMutations, function(m) {
        sampledPaths <- attr(m, "paths")
        sampledTree <- attr(sampledPaths, "tree")
        sampledTree[["tip.label"]]
    })
    cat("Summarizing...\n")
    flush.console()
    pb <- txtProgressBar(min = 0,
                         max = length(allMutations),
                         style = 3)
    # The tip names grouped by ancestral node of the original tree
    originalNodeTips <- lapply(nodeAlign, function(nd) {
        tree[["tip.label"]][as.integer(names(nd))]
    })
    # Assign amino acid for each ancestral node from the sampling result
    assembled <- summarizeAA(
        allMutations = allMutations,
        allSampledTips = allSampledTips,
        originalNodeTips = originalNodeTips,
        setTxtProgressBar = setTxtProgressBar,
        pb = pb
    )
    close(pb)
    # Split the result by node and summarize fixed amino acid of all samples
    res <- lapply(names(assembled), function(site) {
        nodeAAdist <- assembled[[site]]
        site <- as.integer(site)
        siteSummary <- lapply(names(nodeAAdist), function(node) {
            nodeTips <- as.integer(names(nodeAlign[[node]]))
            attr(nodeTips, "samplingSummary") <- nodeAAdist[[node]]
            attr(nodeTips, "aaSummary") <-
                tableAA(nodeAlign[[node]], site - 1)
            return(nodeTips)
        })
        names(siteSummary) <- names(nodeAAdist)
        return(siteSummary)
    })
    names(res) <- names(assembled)
    return(res)
}

.summarizeSamplingAA <- function(nodeTips) {
    tipsAA <- attr(nodeTips, "aaSummary")
    samplingAA <- attr(nodeTips, "samplingSummary")
    if (length(tipsAA) == 1 &&
        names(tipsAA) %in% names(samplingAA)) {
        return(names(tipsAA))
    }
    overlappedAA <- intersect(names(tipsAA), names(samplingAA))
    if (length(overlappedAA) != 0) {
        tipsAA <- tipsAA[overlappedAA]
        return(names(tipsAA)[which.max(tipsAA)])
    } else {
        return(names(samplingAA)[which.max(samplingAA)])
    }
}

.assembleFixation <- function(x, tree, align, paths, divNodes) {
    # Divergent nodes are not included anywhere in the result
    noDivNodesPaths <- lapply(paths, function(p) {
        as.character(setdiff(p, divNodes))
    })
    res <- list()
    for (site in names(x)) {
        summarized <- x[[site]]
        for (path in noDivNodesPaths) {
            # Only when all nodes in a path is covered will it be considered
            if (all(path %in% names(summarized))) {
                seg <- list()
                # Initiate aggregating the nodes with the same amino acids
                # fixed
                previousNode <- path[1]
                nodeTips <- summarized[[previousNode]]
                previousAA <- .summarizeSamplingAA(nodeTips)
                attr(nodeTips, "AA") <- previousAA
                attr(nodeTips, "node") <- previousNode
                seg[[previousNode]] <- nodeTips
                # Iterate through the remaining nodes
                for (node in path[-1]) {
                    nodeTips <- summarized[[node]]
                    currentAA <- .summarizeSamplingAA(nodeTips)
                    # Aggregate the tips if the current fixed amino acid
                    # is the same as the previous
                    if (currentAA == previousAA) {
                        node <- previousNode
                        nodeTips <- c(seg[[node]], nodeTips)
                        attr(nodeTips, "AA") <- currentAA
                    } else {
                        attr(nodeTips, "AA") <- currentAA
                    }
                    seg[[node]] <- nodeTips
                    attr(seg[[node]], "node") <- node
                    previousAA <- currentAA
                    previousNode <- node
                }
                targetIndex <- NULL
                if (length(seg) < 2) {
                    next
                } else if (!site %in% names(res)) {
                    targetIndex <- 1
                } else {
                    # Assume a new fixation path is to add.
                    targetIndex <- length(res[[site]]) + 1
                    # The index to extract the terminal tips of the fixation.
                    endIndex <- length(seg)
                    # Some site may have multiple fixation on multiple
                    # lineages. The following is for deciding at which
                    # index should it be assigned in the 'res[[site]]'
                    # Retrieve the existing mutation path of the site
                    existPath <- res[[site]]
                    # Which mutaiton path has the same mutations as 'seg'
                    toCombine <- vapply(
                        X = existPath,
                        FUN = function(ep) {
                            identical(lapply(seg, c)[-endIndex],
                                      lapply(ep, c)[-endIndex]) &&
                                currentAA == attr(ep[[endIndex]], "AA")
                        },
                        FUN.VALUE = logical(1)
                    )
                    if (any(toCombine)) {
                        existIndex <- which(toCombine)
                        # These are the candidates to combine. The additional
                        # condition be all the descendant tips are included.
                        toCombine <- unlist(lapply(
                            X = c(res[[site]][existIndex], list(seg)),
                            FUN = "[[",
                            ... = endIndex
                        ))
                        allTips <-
                            .childrenTips(tree, getMRCA(tree, toCombine))
                        if (all(allTips %in% toCombine)) {
                            seg[[endIndex]] <- toCombine
                            attr(seg[[endIndex]], "AA") <- currentAA
                            res[[site]] <- res[[site]][-existIndex]
                            targetIndex <- length(res[[site]]) + 1
                        } else {
                            # Add new fixation for the site if no existing
                            # mutation path can be combined with
                            targetIndex <- length(existPath) + 1
                        }
                    }
                }
                if (is.null(targetIndex)) {
                    next
                }
                i <- as.integer(site)
                seg <- lapply(seg, function(tips) {
                    attr(tips, "tipsAA") <- substr(x = align[tips],
                                                   start = i,
                                                   stop = i)
                    return(tips)
                })
                res[[site]][[targetIndex]] <- seg
                attr(res[[site]], "site") <- i
                attr(res[[site]], "tree") <- tree
                class(res[[site]]) <- "sitePath"
            }
        }
    }
    return(res)
}

#' @rdname findSites
#' @name multiFixationSites
#' @description
#' After finding the \code{\link{lineagePath}} of a phylogenetic tree,
#' \code{multiFixationSites} uses random sampling on the original tree
#' and applies the method used in \code{fixationSites} to each sampled
#' tree and summarize the results from all the samples.
#' @param samplingSize
#' The number of tips sampled for each round of resampling. It shoud be
#' at least 10th and at most nine 10ths of the tree tips.
#' @param samplingTimes
#' The total times of random sampling to do. It should be greater than 100.
#' @return
#' \code{multiFixationSites} returns sites with multiple fixations.
#' @importFrom ape drop.tip
#' @importFrom utils flush.console
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
multiFixationSites.lineagePath <- function(paths,
                                           samplingSize = NULL,
                                           samplingTimes = 100,
                                           minEffectiveSize = 0,
                                           searchDepth = 1,
                                           method = c("compare", "insert", "delete"),
                                           ...) {
    # Get the tree and aligned sequences
    tree <- attr(paths, "tree")
    nTips <- length(tree[["tip.label"]])
    align <- attr(paths, "align")
    # Check the parameters for resampling
    if (is.null(samplingSize)) {
        samplingSize <- nTips / 2
    } else if (!is.numeric(samplingSize)) {
        stop("\"samplingSize\" only accept numeric")
    } else if (samplingSize > 9 * nTips / 10 ||
               samplingSize < nTips / 10) {
        stop("\"samplingSize\" should be within",
             " one 10th and nine 10ths of tree size")
    } else {
        samplingSize <- as.integer(samplingSize)
    }
    if (samplingTimes < 100) {
        warning("\"samplingTimes\" is preferably over 100.")
    }
    if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    } else {
        minEffectiveSize <- ceiling(minEffectiveSize)
    }
    # Get the 'searchDepth' for heuristic search
    if (searchDepth < 1) {
        stop("\"searchDepth\" should be at least 1")
    } else {
        searchDepth <- ceiling(searchDepth)
    }
    # Decide which miniminzing strategy
    minimizeEntropy <- switch(
        match.arg(method),
        "compare" = minEntropyByComparing,
        "insert" = minEntropyByInserting,
        "delete" = minEntropyByDeleting
    )
    # Generate the site mapping from reference
    reference <- attr(paths, "reference")
    # Extend the path
    paths <- .extendPaths(paths, tree)
    divNodes <- divergentNode(paths)
    # Get the name and sequence of the children tips that
    # are descendant of the nodes.
    nodeAlign <- .tipSeqsAlongPathNodes(
        paths = paths,
        divNodes = divNodes,
        tree = tree,
        align = align
    )
    cat("\nResampling...\n")
    flush.console()
    # The resampling process
    pb <- txtProgressBar(min = 0,
                         max = samplingTimes,
                         style = 3)
    allMutations <- lapply(seq_len(samplingTimes), function(iter) {
        # Sampled tree
        sampledTree <-
            drop.tip(tree, sample(seq_len(nTips), nTips - samplingSize))
        # The matching between original and sampled tips
        original2sampled <-
            match(sampledTree[["tip.label"]], tree[["tip.label"]])
        # Sampled aligned sequence
        sampledAlign <- align[original2sampled]
        # The corresponding path on the sampled tree subject to
        # path on original tree
        sampledPaths <- .originalPath2sampled(
            nodeAlign = nodeAlign,
            original2sampled = original2sampled,
            sampledTree = sampledTree
        )
        # Find the fixation site for the sampled tree
        sampledDivNodes <- divergentNode(sampledPaths)
        sampledNodeAlign <- attr(sampledPaths, "nodeAlign")
        res <- .findFixationSite(
            paths = sampledPaths,
            tree = sampledTree,
            align = sampledAlign,
            nodeAlign = sampledNodeAlign,
            divNodes = sampledDivNodes,
            reference = reference,
            minimizeEntropy = minimizeEntropy,
            minEffectiveSize = minEffectiveSize,
            searchDepth = searchDepth
        )
        attr(sampledPaths, "tree") <- sampledTree
        attr(res, "paths") <- sampledPaths
        setTxtProgressBar(pb, iter)
        return(res)
    })
    close(pb)
    # Summarize the fixed amino acid for tips from resampling
    sampleSummary <- .sampleSummarize(allMutations, nodeAlign, tree)
    # Rebuild the fixation path
    res <- .assembleFixation(
        x = sampleSummary,
        tree = tree,
        align = align,
        paths = paths,
        divNodes = divNodes
    )
    attr(res, "paths") <- paths
    attr(res, "reference") <- reference
    class(res) <- "multiFixationSites"
    return(res)
}

#' @export
multiFixationSites <- function(paths, ...)
    UseMethod("multiFixationSites")

#' @export
print.multiFixationSites <- function(x, ...) {
    cat("Result for", length(attr(x, "paths")), "paths:\n\n")
    if (length(x) == 0) {
        cat("No multi-fixation found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(attr(x, "reference"), "refSeqName")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}
