#' @rdname sitesMinEntropy
#' @name sitesMinEntropy
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
        minEffectiveSize <-
            length(tree[["tip.label"]]) / length(unique(unlist(paths)))
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    minEffectiveSize <- ceiling(minEffectiveSize)
    # Set the search depth for heuristic search
    if (searchDepth < 1) {
        stop("\"searchDepth\" should be at least 1")
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
    nodeAlign <- .tipSeqsAlongPathNodes(paths = paths,
                                        divNodes = divNodes)
    # Get the MSA numbering
    reference <- attr(paths, "msaNumbering")
    align <- attr(paths, "align")
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    # In case root node does not have any tips (because itself is a divergent
    # node)
    excludedNodes <- divNodes
    rootNode <- attr(paths, "rootNode")
    if (!rootNode %in% names(nodeAlign)) {
        excludedNodes <- c(rootNode, excludedNodes)
    }
    # Turn the site number into index for C++ code
    siteIndices <- reference[loci] - 1
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
    clustersByPath <- .clustersByPath(res)
    clustersByPath <- .mergeClusters(clustersByPath)
    attr(res, "clustersByPath") <-
        .assignClusterNames(clustersByPath)
    attr(res, "paths") <- paths
    class(res) <- "sitesMinEntropy"
    return(res)
}

.tipSeqsAlongPathNodes <- function(paths, divNodes) {
    tree <- attr(paths, "tree")
    align <- attr(paths, "align")
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
            childrenNode <- tree[["edge"]][which(tree[["edge"]][, 1] == n), 2]
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

.childrenTips <- function(tree, node) {
    maxTip <- length(tree[["tip.label"]])
    children <- integer()
    getChildren <- function(edges, parent) {
        children <<- c(children, parent[which(parent <= maxTip)])
        i <- which(edges[, 1] %in% parent)
        if (length(i) == 0L) {
            return(children)
        } else {
            parent <- edges[i, 2]
            return(getChildren(edges, parent))
        }
    }
    return(getChildren(tree[["edge"]], node))
}

.runEntropyMinimization <- function(siteIndex,
                                    path,
                                    nodeAlign,
                                    minimizeEntropy,
                                    minEffectiveSize,
                                    searchDepth) {
    # Assign a variable to store the tip names and their info on amino
    # acids. They are the potential fixation segment
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
        # Attach the node to the previous node if they're both purely fixed
        # and have the same AA fixed.
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
            # Oddly, R uses the name of the first variable when adding two
            # numeric vectors. So there is no need for names (AA) assignment
        }
        # Assign or re-assign the 'nodeTips' with 'aaSummary' to the
        # 'nodeSummaries'
        nodeSummaries[[node]] <- nodeTips
        previousAA <- currentAA
        previousNode <- node
    }
    # Return empty value if the site is purely fixed on the lineage
    seg <- list()
    if (length(nodeSummaries) >= 2) {
        seg <- minimizeEntropy(nodeSummaries,
                               minEffectiveSize,
                               searchDepth)
    }
    return(seg)
}

.clustersByPath <- function(fixations) {
    # Find the clustering for each lineage path
    res <- lapply(fixations, function(sp) {
        # Remove the site purely conserved on the lineage
        group <- list()
        for (site in names(sp)) {
            mp <- sp[[site]]
            if (length(mp) >= 2) {
                toAdd <- lapply(mp, function(tips) {
                    siteChar <- attr(tips, "AA")
                    attributes(tips) <- NULL
                    attr(tips, "site") <- siteChar
                    names(attr(tips, "site")) <- site
                    tips
                })
                group <- c(group, list(toAdd))
            }
        }
        if (length(group) == 0) {
            return(group)
        }
        # Group tips according to fixation points
        res <- group[[1]]
        for (p in group[-1]) {
            for (tips in p) {
                site <- attr(tips, "site")
                # Update grouping for each tips by growing a new list
                newGrouping <- list()
                for (i in seq_along(res)) {
                    gp <- res[[i]]
                    common <- sort(intersect(tips, gp))
                    # No new cluster when the coming tips have no overlap or are
                    # identical to tips in an existing cluster
                    if (length(common) == 0) {
                        newGrouping <- res[seq_len(i)]
                    } else if (identical(sort(gp), sort(tips))) {
                        attr(gp, "site") <- c(attr(gp, "site"), site)
                        if (i + 1 <= length(res)) {
                            trailing <- res[(i + 1):length(res)]
                        } else {
                            trailing <- list()
                        }
                        newGrouping <-
                            c(newGrouping, list(gp), trailing)
                        break
                    } else {
                        # A new cluster formed when there is overlapped between
                        # new coming tips and existing tips in a cluster
                        if (identical(sort(gp), common)) {
                            # The new coming tips includes the current group
                            # The extra tips stay for the next loop
                            tips <- setdiff(tips, gp)
                            # Update the SNP site info for the current group
                            attr(gp, "site") <-
                                c(attr(gp, "site"), site)
                            newGrouping <- c(newGrouping, list(gp))
                        } else if (identical(sort(tips), common)) {
                            # The new coming tips are included in the group
                            # (they are used up at this point)
                            separate <- setdiff(gp, tips)
                            attributes(separate) <- attributes(gp)
                            attr(tips, "site") <-
                                c(attr(gp, "site"), site)
                            if (i + 1 <= length(res)) {
                                trailing <- res[(i + 1):length(res)]
                            } else {
                                trailing <- list()
                            }
                            newGrouping <- c(newGrouping,
                                             list(tips),
                                             list(separate),
                                             trailing)
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
    })
    res <- res[which(lengths(res) != 0)]
    return(res)
}

.mergeClusters <- function(clustersByPath) {
    # Find the divergent point and remove the overlapped part
    res <- list(clustersByPath[[1]])
    # 'res' stores all the non-overlapped parts which means all the clusters are
    # unique
    for (gpIndex in seq_along(clustersByPath)[-1]) {
        # 'gp' is the complete path with overlapped parts
        gp <- clustersByPath[[gpIndex]]
        # Reset the variables to NULL for in case of no overlap
        t <- integer()
        # The index of 'res' which to merge with 'gp'
        toMergeIndex <- NULL
        # The index of 'gp' where the divergent point is. Each truncated 'gp' in
        # 'res' will have one but only the deepest will be used
        divergedIndex <- 0L
        # The number of shared tips at divergent point will be used to decide if
        # the two clusters are completely diverged or not
        sharedAtDiv <- integer()
        # Loop through 'res' to find the most related group
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
                # overlap with the 'gp'
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
        refSites <- attr(divergedTips, "site")
        # The non-shared part of the 'divergedTips'. This part will not be empty
        divergedTips <- setdiff(divergedTips,
                                unlist(clustersByPath[[toMergeIndex]]))
        attr(divergedTips, "site") <- refSites
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
                    # When 'sharedTips' is not empty, the fixation site should
                    # be the only info to give back to
                    attr(sharedTips, "site") <-
                        attr(toMerge[[i]], "site")
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
                if (identical(attr(grouping[[i]][[j]], "site"),
                              attr(grouping[[i]][[j + 1]], "site"))) {
                    nextMajor <- currMajor + 1
                } else {
                    # Create a new major number when encounter a divergent point
                    currMajor <- currMajor + 1L
                    # Reset mini number
                    currMini <- 1L
                    nextMajor <- currMajor
                }
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

#' @export
sitesMinEntropy <- function(x, ...) {
    UseMethod("sitesMinEntropy")
}

#' @export
print.sitesMinEntropy <- function(x, ...) {
    cat("This is a 'sitesMinEntropy' object.", "\n")
}
