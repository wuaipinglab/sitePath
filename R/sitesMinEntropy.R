#' @importFrom utils tail
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster stopCluster
#' @importFrom parallel parLapply
#' @importFrom ape getMRCA

#' @rdname sitesMinEntropy
#' @title Fixation sites prediction
#' @description After finding the \code{\link{lineagePath}} of a phylogenetic
#'   tree, \code{sitesMinEntropy} perform entropy minimization on every site of
#'   the sequence to group the tips according to amino acid/nucleotide.
#' @param x A \code{lineagePath} object returned from \code{\link{lineagePath}}
#'   function.
#' @param minEffectiveSize The minimum number of tips in a group.
#' @param searchDepth The function uses heuristic search but the termination of
#'   the search cannot be intrinsically decided. \code{searchDepth} is needed to
#'   tell the search when to stop.
#' @param method The strategy for predicting the fixation. The basic approach is
#'   entropy minimization and can be achieved by adding or removing fixation
#'   point, or by comparing the two.
#' @param ... further arguments passed to or from other methods.
#' @return A \code{sitesMinEntropy} object.
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' sitesMinEntropy(lineagePath(tree))
sitesMinEntropy <- function(x, ...) {
    UseMethod("sitesMinEntropy")
}

#' @rdname sitesMinEntropy
#' @export
sitesMinEntropy.lineagePath <- function(x,
                                        minEffectiveSize = NULL,
                                        searchDepth = 1,
                                        method = c("compare", "insert", "delete"),
                                        ...) {
    paths <- .phyMSAmatch(x)
    # Set the minimal size of the group during the search
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- attr(x, "minSize")
    } else if (!is.numeric(minEffectiveSize) ||
               minEffectiveSize < 0) {
        stop("'minEffectiveSize' (",
             minEffectiveSize,
             ") is not positive numeric")
    } else {
        minEffectiveSize <- ceiling(minEffectiveSize)
    }
    # Set the search depth for heuristic search
    if (searchDepth < 1) {
        stop("'searchDepth' (", searchDepth, ")  is less than 1")
    } else {
        searchDepth <- ceiling(searchDepth)
    }
    # Decide which minimizing strategy
    minimizeEntropy <- switch(
        match.arg(method),
        "compare" = minEntropyByComparing,
        "insert" = minEntropyByInserting,
        "delete" = minEntropyByDeleting
    )
    # Get the divergent nodes
    divNodes <- divergentNode(paths)
    # The tips and matching
    pathNodeTips <- .tipSeqsAlongPathNodes(paths, divNodes)
    # In case root node does not have any tips
    excludedNodes <- divNodes
    rootNode <- attr(paths, "rootNode")
    if (!rootNode %in% names(pathNodeTips)) {
        excludedNodes <- c(rootNode, excludedNodes)
    }
    # Exclude the invariant sites
    loci <- attr(x, "loci")
    # Turn the site number into index for C++ code
    siteIndices <- attr(paths, "msaNumbering")[loci] - 1
    names(siteIndices) <- as.character(loci)
    # Group the result by path for all loci
    pathsWithSeqs <- lapply(paths, function(path) {
        path <- as.character(setdiff(path, excludedNodes))
        names(path) <- path
        pathNodeAlign <- pathNodeTips[path]
        # The node changed to the previous node if the previous node is the
        # divergent node. This is for the mutation labeling
        for (currIndex in seq_along(pathNodeAlign)[-1]) {
            prevIndex <- currIndex - 1
            prevNode <- names(pathNodeAlign)[prevIndex]
            if (prevNode %in% divNodes) {
                names(pathNodeAlign)[currIndex] <- prevNode
            }
        }
        attr(pathNodeAlign, "pathTipNum") <-
            sum(lengths(pathNodeAlign))
        return(pathNodeAlign)
    })
    pathTipNums <-
        vapply(pathsWithSeqs, attr, integer(1), "pathTipNum")
    # Use the number of tips to order the paths. The path with fewer number of
    # tips come first so its result will be prioritize when merging
    tipNumRank <- order(pathTipNums)
    # The max number of tips a path has
    maxPathTipNum <- tail(pathTipNums[tipNumRank], 1)
    # Test if multiprocessing is turned on
    mc <- getOption("cl.cores")
    if (is.null(mc)) {
        res <- lapply(
            X = pathsWithSeqs[tipNumRank],
            FUN = function(pathNodeAlign) {
                scaledSize <- attr(pathNodeAlign, "pathTipNum") / maxPathTipNum
                # The path with fewer tips will use smaller threshold
                scaledSize <- ceiling(minEffectiveSize * scaledSize)
                attr(pathNodeAlign, "scaledSize") <- scaledSize
                segs <- lapply(
                    X = siteIndices,
                    FUN = .runEntropyMinimization,
                    pathNodeAlign = pathNodeAlign,
                    minimizeEntropy = minimizeEntropy,
                    searchDepth = searchDepth
                )
                attr(segs, "pathNodeTips") <- pathNodeAlign
                return(segs)
            }
        )
        # Calibrate the result from all paths
        res <- .unifyEntropyGrouping(res, paths, NULL)
        # Cluster tips according to fixation sites
        clustersByPath <- lapply(res, function(segs) {
            # Set the site and amino acid/nucleotide info for each group of tips
            group <- lapply(names(segs), function(site) {
                lapply(segs[[site]], function(tips) {
                    siteChar <- attr(tips, "AA")
                    names(siteChar) <- site
                    node <- attr(tips, "node")
                    # Purge the attributes and keep only the node and amino
                    # acid/nucleotide info
                    attributes(tips) <- NULL
                    attr(tips, "AA") <- siteChar
                    attr(tips, "node") <- node
                    return(tips)
                })
            })
            # This will group the tips of the lineage path and the adjacent
            # groups will have at least one fixed site different
            group <- .clusterByFixation(group)
            attr(group, "pathNodeTips") <-
                attr(segs, "pathNodeTips")
            return(group)
        })
    } else {
        cl <- .createCluster(mc, method = FALSE)
        res <- lapply(
            X = pathsWithSeqs[tipNumRank],
            FUN = function(pathNodeAlign) {
                scaledSize <- attr(pathNodeAlign, "pathTipNum") / maxPathTipNum
                # The path with fewer tips will use smaller threshold
                scaledSize <- ceiling(minEffectiveSize * scaledSize)
                attr(pathNodeAlign, "scaledSize") <- scaledSize
                # Entropy minimization result for every locus
                segs <- parLapply(
                    cl = cl,
                    X = siteIndices,
                    fun = .runEntropyMinimization,
                    pathNodeAlign = pathNodeAlign,
                    minimizeEntropy = minimizeEntropy,
                    searchDepth = searchDepth
                )
                attr(segs, "pathNodeTips") <- pathNodeAlign
                return(segs)
            }
        )
        # Calibrate the result from all paths
        res <- .unifyEntropyGrouping(res, paths, cl)
        # Cluster tips according to fixation sites
        clustersByPath <- parLapply(
            cl = cl,
            X = res,
            fun = function(segs) {
                # Set the site and amino acid/nucleotide info for each group of
                # tips
                group <- lapply(names(segs), function(site) {
                    lapply(segs[[site]], function(tips) {
                        siteChar <- attr(tips, "AA")
                        names(siteChar) <- site
                        node <- attr(tips, "node")
                        # Purge the attributes and keep only the node and amino
                        # acid/nucleotide info
                        attributes(tips) <- NULL
                        attr(tips, "AA") <- siteChar
                        attr(tips, "node") <- node
                        return(tips)
                    })
                })
                # This will group the tips of the lineage path and the adjacent
                # groups will have at least one fixed site different
                group <- .clusterByFixation(group)
                attr(group, "pathNodeTips") <-
                    attr(segs, "pathNodeTips")
                return(group)
            }
        )
        stopCluster(cl)
        cat(paste("Multiprocessing ended.\n"))
    }
    clustersByPath <- .mergeClusters(clustersByPath)
    attr(res, "clustersByPath") <-
        .assignClusterNames(clustersByPath)
    attr(res, "paths") <- paths
    class(res) <- "sitesMinEntropy"
    return(res)
}

.runEntropyMinimization <- function(siteIndex,
                                    pathNodeAlign,
                                    minimizeEntropy,
                                    searchDepth) {
    minEffectiveSize <- attr(pathNodeAlign, "scaledSize")
    # Assign a variable to store the tip names and their info on amino acids.
    # They are the potential fixation segment
    nodeTips <- integer()
    previousAA <- NULL
    currentAA <- NULL
    previousNode <- NULL
    # The input for entropy minimization calculation
    nodeSummaries <- list()
    # Divergent nodes are not included anywhere in the result
    for (node in names(pathNodeAlign)) {
        # Get the related descendant tips and related sequences
        nodeTips <- pathNodeAlign[[node]]
        # Frequency of the amino acids at the locus
        aaSummary <- tableAA(attr(nodeTips, "align"), siteIndex)
        # Associate the amino acid frequency with the tip names
        attr(nodeTips, "aaSummary") <- aaSummary
        # Decide the current fixed amino acid
        if (length(aaSummary) == 1) {
            currentAA <- names(aaSummary)
        } else {
            currentAA <- NULL
        }
        # Attach the node to the previous node if they're both purely fixed and
        # have the same AA fixed.
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
            # R uses the name of the first vector variable when adding two
            # numeric vectors. So there is no need for names (AA) assignment
        }
        # Assign or re-assign the 'nodeTips' with 'aaSummary' to the
        # 'nodeSummaries'
        nodeSummaries[[node]] <- nodeTips
        previousAA <- currentAA
        previousNode <- node
    }
    # Return empty value if the site is purely fixed on the lineage
    if (length(nodeSummaries) >= 2) {
        seg <- minimizeEntropy(nodeSummaries,
                               minEffectiveSize,
                               searchDepth)
        names(seg) <- vapply(seg, attr, character(1), "node")
    } else {
        seg <- nodeSummaries
        attr(seg[[1]], "AA") <- names(attr(seg[[1]], "aaSummary"))
        attr(seg[[1]], "node") <- names(seg)
    }
    return(seg)
}

.unifyEntropyGrouping <- function(res, paths, cl) {
    align <- attr(paths, "align")
    if (is.null(cl)) {
        allMergedGroupings <- lapply(
            X = names(res[[1]]),
            FUN = .calibrateFixedAA,
            fixations = res,
            align = align
        )
    } else {
        allMergedGroupings <- parLapply(
            cl = cl,
            X = names(res[[1]]),
            fun = .calibrateFixedAA,
            fixations = res,
            align = align
        )
    }
    names(allMergedGroupings) <- names(res[[1]])
    # Iterate each locus
    for (locus in names(res[[1]])) {
        # The index for C++
        locusIndex <- as.integer(locus) - 1
        # Get the merged groupings of each path for the locus
        mergedGroupings <- allMergedGroupings[[locus]]
        # Link the merged groupings for each lineage path to calibrate and
        # recover the full groupings. The linked merged groupings are stored for
        # the remaining paths.
        linkedMerged <- list()
        for (pathIndex in seq_along(mergedGroupings)) {
            currLinked <- .findGroupingLinkage(
                pathIndex = pathIndex,
                mergedGroupings = mergedGroupings,
                linkedMerged = linkedMerged
            )
            # Store the linked groupings back tracing to the root for the
            # remaining groupings
            linkedMerged <- c(linkedMerged, list(currLinked))
            # The merged result for the current grouping needs to be included
            # for reconstructing but not needed for other groupings
            currLinked[[pathIndex]] <- mergedGroupings[[pathIndex]]
            # Make the tips clusters into one list
            currLinked <- currLinked[which(lengths(currLinked) > 0)]
            currLinked <- unlist(currLinked, recursive = FALSE)
            # Reconstruct the grouping of the path
            reconstructed <- list()
            clusterNodes <- character()
            prevGP <- currLinked[[1]]
            for (gpIndex in seq_along(currLinked)[-1]) {
                currGP <- currLinked[[gpIndex]]
                if (attr(prevGP, "AA") == attr(currGP, "AA")) {
                    currGP <- c(prevGP, currGP)
                    attr(currGP, "aaSummary") <-
                        tableAA(align[currGP], locusIndex)
                    attr(currGP, "AA") <- attr(prevGP, "AA")
                    attr(currGP, "node") <- attr(prevGP, "node")
                } else {
                    # Add a new tip cluster only when there is a change of fixed
                    # amino acid/nucleotide
                    reconstructed <- c(reconstructed, list(prevGP))
                    clusterNodes <-
                        c(clusterNodes, attr(prevGP, "node"))
                }
                prevGP <- currGP
            }
            reconstructed <- c(reconstructed, list(prevGP))
            clusterNodes <- c(clusterNodes, attr(prevGP, "node"))
            # Re-assign the merged result to the original
            names(reconstructed) <- clusterNodes
            res[[pathIndex]][[locus]] <- reconstructed
        }
    }
    return(res)
}

.calibrateFixedAA <- function(locus, fixations, align) {
    # The original entropy minimization result for the locus
    unMergedGroupings <- lapply(fixations, function(segs) {
        res <- segs[[locus]]
        attr(res, "pathNodeTips") <- attr(segs, "pathNodeTips")
        return(res)
    })
    # Calibrate the tip grouping result from all result
    res <- .mergeClusters(unMergedGroupings)
    locusIndex <- as.integer(locus) - 1
    toMergePrevAA <- character()
    # Special case when the root group is the divergent point
    tips <- res[[1]][[1]]
    # The dominant 'AA' of all direct descendant tip groups
    if (!is.null(attr(tips, "toMerge"))) {
        refAA <- unique(vapply(unMergedGroupings, function(seg) {
            return(attr(seg[[1]], "AA"))
        }, character(1)))
        aaSummary <- tableAA(align[tips], locusIndex)
        toMergePrevAA[1] <- refAA[1]
        if (any(names(aaSummary) %in% refAA)) {
            toMergePrevAA[1] <- names(which.max(aaSummary[refAA]))
        }
        attr(res[[1]][[1]], "AA") <- toMergePrevAA[1]
    }
    for (pathIndex in seq_along(res)) {
        mGrouping <- res[[pathIndex]]
        # Iterate to find the tip groups at divergent point and find the most
        # fitting amino acid/nucleotide
        for (gIndex in seq_along(mGrouping)) {
            tips <- mGrouping[[gIndex]]
            toMerge <- attr(tips, "toMerge")
            if (!is.null(toMerge)) {
                otherIndex <- as.integer(names(toMerge))
                toMergePrevAA[otherIndex] <- attr(tips, "AA")
                # Calibrate the fixed amino acid/nucleotide at the divergent
                # point after merging
                originalAA <- vapply(
                    X = unMergedGroupings[otherIndex],
                    FUN = function(seg) {
                        majorFixedAA <- character()
                        maxOverlapNum <- 0
                        # To find the fixed amino acid/nucleotide from the
                        # original entropy minimum result for each path
                        for (original in seg) {
                            overlapNum <- intersect(tips, original)
                            if (length(overlapNum) > maxOverlapNum) {
                                majorFixedAA <- attr(original, "AA")
                                maxOverlapNum <- length(overlapNum)
                            }
                        }
                        return(majorFixedAA)
                    },
                    FUN.VALUE = character(1)
                )
                # All originally assigned 'AA' for the current tip group
                refAA <- unique(c(attr(tips, "AA"), originalAA))
                # In case there is a disagreement between the paths
                if (length(refAA) > 1) {
                    # The remaining tips of the group split at the divergent
                    # point on the current and other paths
                    nextTips <- mGrouping[[gIndex + 1]]
                    nextFixedAA <- attr(nextTips, "AA")
                    for (toMergeIndex in otherIndex) {
                        otherTips <- res[[toMergeIndex]][[1]]
                        nextFixedAA <- c(nextFixedAA,
                                         attr(otherTips, "AA"))
                    }
                    # The fixed amino acid/nucleotide in the divergent point
                    # that does not come from the splitting of the following tip
                    # groups but is rather related to the previous tip group
                    extraFixedAA <- setdiff(refAA, nextFixedAA)
                    if (length(extraFixedAA) > 1) {
                        # This might be quite impossible but just in case
                        refAA <- extraFixedAA[1]
                    } else if (length(extraFixedAA) == 0) {
                        if (gIndex == 1) {
                            # The previous 'AA' of the first group (after
                            # merging) is a little tricky to find
                            refAA <- toMergePrevAA[pathIndex]
                        } else {
                            prevTips <- mGrouping[[gIndex - 1]]
                            aaSummary <-
                                tableAA(align[c(tips, prevTips)], locusIndex)
                            consistentAA <-
                                intersect(refAA, names(aaSummary))
                            if (length(consistentAA) == 0) {
                                # The previous fixed amino acid/nucleotide
                                refAA <- attr(prevTips, "AA")
                            } else if (length(consistentAA) > 1) {
                                # The dominant amino acid/nucleotide
                                refAA <-
                                    names(which.max(aaSummary[refAA]))
                            } else {
                                refAA <- consistentAA
                            }
                        }
                    } else {
                        refAA <- extraFixedAA
                    }
                }
                attr(res[[pathIndex]][[gIndex]], "AA") <- refAA
            }
        }
    }
    return(res)
}

.findGroupingLinkage <- function(pathIndex,
                                 mergedGroupings,
                                 linkedMerged) {
    # The linked groupings for the current path The index of an irrelevant path
    # will hold a NULL while the index of a relevant path will hold the grouping
    # up till the tip cluster with the 'toMerge' index indicating the divergent
    # point.
    currLinked <- rep(list(), len = pathIndex)
    # Iterate the merged grouping of each path to find the 'toMerge' index same
    # as the current path index
    for (mergedIndex in seq_along(linkedMerged)) {
        # The merged grouping that has no overlap
        mGrouping <- mergedGroupings[[mergedIndex]]
        # Iterate each tip cluster in the merged grouping
        for (gIndex in seq_along(mGrouping)) {
            tips <- mGrouping[[gIndex]]
            toMergeIndex <- as.integer(names(attr(tips, "toMerge")))
            # Trace back to the root when the path index of current grouping
            # is found in 'toMergeIndex'
            if (pathIndex %in% toMergeIndex) {
                # The linked merged groupings trace back to the root
                currLinked <- linkedMerged[[mergedIndex]]
                # Attach the index that directly links to the current grouping
                currLinked[[mergedIndex]] <-
                    mGrouping[seq_len(gIndex)]
                # The result of the groups and grouping will be ignored as there
                # could only be one relevant divergent point
                return(currLinked)
            }
        }
    }
    return(currLinked)
}

.clusterByFixation <- function(group) {
    # Group tips according to fixation points
    res <- group[[1]]
    for (seg in group[-1]) {
        # the node name for each group of tips
        nodeNames <- names(seg)
        # Iterate each group of tips to contribute to grouping
        for (n in seq_along(seg)) {
            tips <- seg[[n]]
            # The site number and its amino acid/nucleotide
            site <- attr(tips, "AA")
            # Update grouping for each tips by growing a new list
            newGrouping <- list()
            # Compare with each group of the current grouping
            for (i in seq_along(res)) {
                gp <- res[[i]]
                common <- sort(intersect(tips, gp))
                # No new cluster when the coming tips have no overlap or are
                # identical to tips in an existing cluster
                if (length(common) == 0) {
                    # Keep the current grouping if the coming group has no
                    # overlap yet
                    newGrouping <- res[seq_len(i)]
                } else if (identical(sort(gp), sort(tips))) {
                    # The only effect here is to add the new 'AA' info to
                    # the group
                    attr(gp, "AA") <- c(attr(gp, "AA"), site)
                    # The groups after the current group to be added
                    newGrouping <- c(newGrouping,
                                     list(gp),
                                     tail(res, length(res) - i))
                    break
                } else {
                    # A new cluster formed when there is overlapped between new
                    # coming tips and existing tips in a cluster
                    if (identical(sort(gp), common)) {
                        # The new coming tips includes the current group. The
                        # extra tips stay for the next loop just in case it has
                        # impact on the grouping
                        tips <- setdiff(tips, gp)
                        # Update the SNP site info for the current group
                        attr(gp, "AA") <- c(attr(gp, "AA"), site)
                        newGrouping <- c(newGrouping, list(gp))
                    } else if (identical(sort(tips), common)) {
                        # The new coming tips are included in the group (they
                        # are used up at this point) and they split the original
                        # grouping
                        separate <- setdiff(gp, tips)
                        attributes(separate) <- attributes(gp)
                        attr(separate, "node") <- nodeNames[n + 1]
                        # 'tips' is the common part and inherit the attributes
                        # of the to-be-split original group
                        attr(tips, "AA") <- c(attr(gp, "AA"), site)
                        attr(tips, "node") <- attr(gp, "node")
                        newGrouping <- c(
                            newGrouping,
                            list(tips),
                            list(separate),
                            tail(res, length(res) - i)
                        )
                        # Go for the next new coming tips
                        break
                    } else {
                        stop("Something's not right")
                    }
                }
            }
            # The new coming tips are used up and update the grouping
            res <- newGrouping
        }
    }
    return(res)
}

.mergeClusters <- function(clustersByPath) {
    # 'res' stores all the non-overlapped parts which means all the clusters are
    # unique but it still splits into each path
    res <- list(clustersByPath[[1]])
    # Find the divergent point and remove the overlapped part
    for (gpIndex in seq_along(clustersByPath)[-1]) {
        # 'gp' is the complete path with overlapped parts
        gp <- clustersByPath[[gpIndex]]
        # Assume there is no overlapped tips
        t <- integer()
        # The index of 'res' which to merge with 'gp'
        toMergeIndex <- NULL
        # The index of 'gp' where the divergent point is. Each truncated 'gp' in
        # 'res' will have one but only the deepest will be used
        divergedIndex <- 0L
        # The number of shared tips at divergent point will be used to decide if
        # the two clusters are completely diverged or not
        sharedAtDiv <- integer()
        # Loop through 'res' to find the most related group ('res' is changed
        # after each iteration)
        for (i in seq_along(res)) {
            # All existing tips in another 'gp' to see if overlapped with tips
            # in the 'gp' to be merged
            allTips <- unlist(clustersByPath[[i]])
            # Because all the tip groups are unique in 'res', the first cluster
            # in 'gp' containing tips that cannot be found is the divergent
            # point
            for (j in seq_along(gp)) {
                # Once a potential divergent point having being found (the tip
                # cluster in 'gp' containing cannot-be-found tips), safeguard
                # the current 'gp' have actual overlap with all tips in 'res'
                # with index 'j'
                if (any(!gp[[j]] %in% allTips) &&
                    any(unlist(gp) %in% unlist(res[[i]]))) {
                    t <- intersect(gp[[j]], allTips)
                    # The deepest and most divergent point, which is decided by
                    # the index. When the index is the same as the previous one,
                    # chose the one with more shared tips
                    if (j > divergedIndex ||
                        (j == divergedIndex &&
                         length(t) >= length(sharedAtDiv))) {
                        toMergeIndex <- i
                        divergedIndex <- j
                        sharedAtDiv <- t
                    }
                    break
                }
            }
        }
        # Find the tips when diverged
        divergedTips <- gp[[divergedIndex]]
        tempAttrs <- attributes(divergedTips)
        # The non-shared part of the 'divergedTips'. This part will not be empty
        divergedTips <- setdiff(divergedTips,
                                unlist(clustersByPath[[toMergeIndex]]))
        attributes(divergedTips) <- tempAttrs
        # Re-assign the ancestral node since the tip group is split
        pathNodeTips <- attr(gp, "pathNodeTips")
        attr(divergedTips, "node") <- .calibrateNode(divergedTips,
                                                     pathNodeTips)
        refSites <- attr(divergedTips, "AA")
        # Add the truncated 'gp' (no overlap) to 'res'
        if (divergedIndex == length(gp)) {
            # No more trailing tips besides the non-shared part
            res[[gpIndex]] <- list(divergedTips)
        } else {
            # Non-shared part plus the trailing part
            res[[gpIndex]] <- c(list(divergedTips),
                                gp[(divergedIndex + 1):length(gp)])
        }
        # Find the most related group of 'gp' in 'res'
        toMerge <- res[[toMergeIndex]]
        # To determine where to add the new group (truncated 'gp'). This wasn't
        # done above just in case the merged part might not be the same for the
        # two paths
        gpTips <- unlist(gp)
        for (i in seq_along(toMerge)) {
            # The divergent point of the most related group in 'res', which
            # might already be different from the original 'gp'
            if (any(!toMerge[[i]] %in% gpTips)) {
                # The non-shared part
                divergedTips <- setdiff(toMerge[[i]], gpTips)
                # Give back the attributes, including fixation site and possible
                # merging info
                attributes(divergedTips) <- attributes(toMerge[[i]])
                # Re-assign the ancestral node since the tip group is split
                pathNodeTips <- attr(clustersByPath[[toMergeIndex]],
                                     "pathNodeTips")
                attr(divergedTips, "node") <-
                    .calibrateNode(divergedTips, pathNodeTips)
                # The shared part
                sharedTips <- setdiff(toMerge[[i]], divergedTips)
                toMergeRefSites <- list()
                toMergeRefSites[[as.character(gpIndex)]] <- refSites
                if (length(sharedTips) == 0) {
                    # There is at least one group of tips before divergence
                    attr(toMerge[[i - 1]], "toMerge") <-
                        c(toMergeRefSites,
                          attr(toMerge[[i - 1]], "toMerge"))
                    sharedTips <- list()
                } else {
                    # When 'sharedTips' is not empty, the site and node should
                    # be the only info to give back to
                    attributes(sharedTips) <-
                        attributes(toMerge[[i]])
                    # The original 'toMerge' info should be taken by the
                    # 'divergedTips'
                    attr(sharedTips, "toMerge") <- toMergeRefSites
                    sharedTips <- list(sharedTips)
                }
                # The divergent part
                if (i == length(toMerge)) {
                    # No more trailing tips besides the non-shared part
                    divergedTips <- list(divergedTips)
                } else {
                    # Non-shared part plus the trailing part
                    divergedTips <- c(list(divergedTips),
                                      toMerge[(i + 1):length(toMerge)])
                }
                # Reform the most related group because the divergent tips might
                # be split
                res[[toMergeIndex]] <- c(toMerge[seq_len(i - 1)],
                                         sharedTips,
                                         divergedTips)
                break
            }
        }
    }
    return(res)
}

.calibrateNode <- function(divergedTips, pathNodeTips) {
    # Just in case the ancestral node is not found
    notFound <- TRUE
    for (node in names(pathNodeTips)) {
        if (all(pathNodeTips[[node]] %in% divergedTips)) {
            return(node)
        }
    }
    if (notFound) {
        stop("Something is wrong finding the ancestral node")
    }
}

.assignClusterNames <- function(grouping) {
    # The starting major numbers of all 'gp' in 'grouping'
    startingMajors <- rep(NA_integer_, length(grouping))
    startingMajors[1] <- 1L
    # The maximum minor number for each major number (so can be continued on a
    # new 'gp')
    maxMinors <- 0L
    # Iterate 'grouping' to assign cluster number
    for (i in seq_along(grouping)) {
        # Get the starting major number
        currMajor <- startingMajors[i]
        # Increase the maximum minor number for 'currMajor'
        maxMinors[currMajor] <- maxMinors[currMajor] + 1L
        # Initiate mini number
        currMini <- 1L
        for (j in seq_along(grouping[[i]])) {
            attr(grouping[[i]][[j]], "clsName") <- paste(currMajor,
                                                         maxMinors[currMajor],
                                                         currMini,
                                                         sep = ".")
            currMini <- currMini + 1
            toMerge <- attr(grouping[[i]][[j]], "toMerge")
            if (!is.null(toMerge)) {
                # Create a new major number when encounter a divergent point
                currMajor <- currMajor + 1L
                # Reset mini number
                currMini <- 1L
                nextMajor <- currMajor
                # Assign the starting major number for the 'gp' to be merged
                toMergeIndex <- as.integer(names(toMerge))
                startingMajors[toMergeIndex] <-
                    rep(nextMajor, length(toMergeIndex))
                # Initiate minor number for the new major number
                maxMinors[nextMajor] <- 1L
            }
        }
    }
    return(grouping)
}
