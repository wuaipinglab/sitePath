#' @importFrom utils tail

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
sitesMinEntropy <- function(x, ...)
    UseMethod("sitesMinEntropy")

#' @rdname sitesMinEntropy
#' @export
sitesMinEntropy.lineagePath <- function(x,
                                        minEffectiveSize = NULL,
                                        searchDepth = 1,
                                        method = c("compare", "insert", "delete"),
                                        ...) {
    paths <- .phyMSAmatch(x)
    tree <- attr(paths, "tree")
    # Set the minimal size of the group during the search
    if (is.null(minEffectiveSize)) {
        tipNum <- length(tree[["tip.label"]])
        pathNodeNum <- length(unique(unlist(paths)))
        minEffectiveSize <- ceiling(tipNum / pathNodeNum)
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
    nodeAlign <- .tipSeqsAlongPathNodes(paths, divNodes)
    # In case root node does not have any tips
    excludedNodes <- divNodes
    rootNode <- attr(paths, "rootNode")
    if (!rootNode %in% names(nodeAlign)) {
        excludedNodes <- c(rootNode, excludedNodes)
    }
    # Exclude the invariant sites
    loci <- attr(x, "loci")
    # Turn the site number into index for C++ code
    siteIndices <- attr(paths, "msaNumbering")[loci] - 1
    names(siteIndices) <- as.character(loci)
    # Group the result by path for all loci
    res <- lapply(paths, function(path) {
        path <- as.character(setdiff(path, excludedNodes))
        # Entropy minimization result for every locus
        lapply(
            X = siteIndices,
            FUN = .runEntropyMinimization,
            path = path,
            nodeAlign = nodeAlign,
            minimizeEntropy = minimizeEntropy,
            minEffectiveSize = minEffectiveSize,
            searchDepth = searchDepth
        )
    })
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
        .clusterByFixation(group)
    })
    clustersByPath <- .mergeClusters(clustersByPath)
    attr(res, "clustersByPath") <-
        .assignClusterNames(clustersByPath)
    attr(res, "paths") <- paths
    class(res) <- "sitesMinEntropy"
    return(res)
}

.runEntropyMinimization <- function(siteIndex,
                                    path,
                                    nodeAlign,
                                    minimizeEntropy,
                                    minEffectiveSize,
                                    searchDepth) {
    # Assign a variable to store the tip names and their info on amino acids.
    # They are the potential fixation segment
    nodeTips <- integer()
    previousAA <- NULL
    currentAA <- NULL
    previousNode <- NULL
    # The input for entropy minimization calculation
    nodeSummaries <- list()
    # Divergent nodes are not included anywhere in the result
    for (node in path) {
        # Get the related descendant tips and related sequences
        nodeTips <- as.integer(names(nodeAlign[[node]]))
        # Frequency of the amino acids at the locus
        aaSummary <- tableAA(nodeAlign[[node]], siteIndex)
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
        attr(seg[[1]], "AA") <-
            names(attr(nodeSummaries[[1]], "aaSummary"))
        attr(seg[[1]], "node") <- names(seg)
    }
    return(seg)
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
                # Once a potential divergent point having being found, safeguard
                # the truncated 'gp' (in 'res') to merge with have actual
                # overlap with the current non-truncated 'gp'
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
                    attr(sharedTips, "toMerge") <-
                        c(toMergeRefSites,
                          attr(sharedTips, "toMerge"))
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
