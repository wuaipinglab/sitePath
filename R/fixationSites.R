#' @rdname fixationSites
#' @name fixationSites
#' @title Fixation sites prediction
#' @description After finding the \code{\link{lineagePath}} of a phylogenetic
#'   tree, \code{fixationSites} uses the result to find those sites that show
#'   fixation on some, if not all, of the lineages. Parallel evolution is
#'   relatively common in RNA virus. There is chance that some site be fixed in
#'   one lineage but does not show fixation because of different sequence
#'   context.
#' @param paths A \code{lineagePath} object returned from
#'   \code{\link{lineagePath}} function.
#' @param minEffectiveSize The minimum number of tips in a group.
#' @param searchDepth The function uses heuristic search but the termination of
#'   the search cannot be intrinsically decided. \code{searchDepth} is needed to
#'   tell the search when to stop.
#' @param method The strategy for predicting the fixation. The basic approach is
#'   entropy minimization and can be achieved by adding or removing fixation
#'   point, or by comparing the two.
#' @param ... further arguments passed to or from other methods.
#' @return A \code{fixationSites} object.
#' @seealso \code{\link{as.data.frame.fixationSites}}
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' fixationSites(lineagePath(tree))
fixationSites.lineagePath <- function(paths,
                                      minEffectiveSize = NULL,
                                      searchDepth = 1,
                                      method = c("compare", "insert", "delete"),
                                      ...) {
    paths <- .phyMSAmatch(paths)
    tree <- attr(paths, "tree")
    align <- attr(paths, "align")
    # Decide which miniminzing strategy
    minimizeEntropy <- switch(
        match.arg(method),
        "compare" = minEntropyByComparing,
        "insert" = minEntropyByInserting,
        "delete" = minEntropyByDeleting
    )
    # Set the 'minEffectiveSize'
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <-
            length(tree[["tip.label"]]) / length(unique(unlist(paths)))
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    minEffectiveSize <- ceiling(minEffectiveSize)
    # Set the 'searchDepth' for heuristic search
    if (searchDepth < 1) {
        stop("\"searchDepth\" should be at least 1")
    } else {
        searchDepth <- ceiling(searchDepth)
    }
    # Get the divergent nodes
    divNodes <- divergentNode(paths)
    # The tips and matching
    nodeAlign <- .tipSeqsAlongPathNodes(paths = paths,
                                        divNodes = divNodes)
    # Find fixation sites
    res <- .findFixationSite(
        paths = paths,
        nodeAlign = nodeAlign,
        divNodes = divNodes,
        minimizeEntropy = minimizeEntropy,
        minEffectiveSize = minEffectiveSize,
        searchDepth = searchDepth
    )
    # Cluster tips according to fixation sites
    clustersByPath <- .transitionClusters(res, paths, tree)
    clustersByPath <- .mergeClusters(clustersByPath)
    # Merge fixation sites on different paths if applicable
    res <- .combineFixations(res, tree, align)
    # Set 'paths' and 'clustersByPath' attributes
    attr(res, "paths") <- paths
    attr(res, "clustersByPath") <-
        .assignClusterNames(clustersByPath)
    class(res) <- "fixationSites"
    return(res)
}

.childrenTips <- function(tree, node) {
    maxTip <- length(tree[["tip.label"]])
    children <- integer(0)
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

.findFixationSite <- function(paths,
                              nodeAlign,
                              divNodes,
                              minimizeEntropy,
                              minEffectiveSize,
                              searchDepth) {
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
    # The variable to store the result from entropy minimization for
    # each path with those purely fixed excluded.
    res <- list()
    # Iterate each path
    rootNode <- attr(paths, "rootNode")
    if (!rootNode %in% names(nodeAlign)) {
        paths <- lapply(paths, "[", -1)
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
                    attr(res[[site]], "tree") <- attr(paths, "tree")
                    class(res[[site]]) <- "sitePath"
                }
            }
        }
    }
    return(res)
}

.transitionClusters <- function(fixations, paths, tree) {
    # Find the clustering for each lineage path
    groupByPath <- lapply(paths, function(p) {
        terminalTips <- .childrenTips(tree, p[length(p)])
        # Group fixation results by path rather than site
        group <- list()
        for (sp in fixations) {
            site <- attr(sp, "site")
            for (mp in sp) {
                tips <- mp[[length(mp)]]
                # Filter if not belong to the path
                if (all(terminalTips %in% tips)) {
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
    groupByPath <- groupByPath[which(lengths(groupByPath) != 0)]
    return(groupByPath)
}

.mergeClusters <- function(groupByPath) {
    # Find the divergent point and remove the overlapped part
    res <- list(groupByPath[[1]])
    # 'res' stores all the non-overlapped parts which means all the clusters are
    # unique
    for (gpIndex in seq_along(groupByPath)[-1]) {
        # 'gp' is the complete path with overlapped parts
        gp <- groupByPath[[gpIndex]]
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
            # All existing tips in the other 'gp' in 'groupByPath' to see if
            # overlapped with tips in the 'gp' to be merged
            allTips <- unlist(groupByPath[[i]])
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
                                unlist(groupByPath[[toMergeIndex]]))
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
        refSeqName <- attr(x, "reference")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}

#' @rdname plotFunctions
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidytree groupOTU
#' @importFrom ggplot2 scale_color_manual guides guide_legend
#' @importFrom ggrepel geom_label_repel
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' paths <- lineagePath(tree)
#' plot(paths)
#' fixations <- fixationSites(paths)
#' plot(fixations)
#' x <- fixationPath(fixations)
#' plot(x)
plot.fixationSites <- function(x,
                               y = TRUE,
                               ...) {
    tree <- as.treedata.fixationSites(x)
    grp <- levels(attr(tree@phylo, "Groups"))
    grp <- grp[which(grp != "0")]
    groupColors <-
        colorRampPalette(brewer.pal(9, "Set1"))(length(grp))
    names(groupColors) <- grp
    groupColors["0"] <- "black"

    p <- ggtree(tree, aes(color = Groups)) +
        scale_color_manual(values = groupColors) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme(legend.position = "left")
    if (y) {
        p <- p + geom_label_repel(
            aes(x = branch, label = SNPs),
            fill = "lightgreen",
            color = "black",
            min.segment.length = 0,
            na.rm = TRUE
        )
    }
    return(p)
}
