#' @rdname fixationSites
#' @name fixationSites
#' @title Fixation sites prediction
#' @description After finding the \code{\link{lineagePath}} of a phylogenetic
#'   tree, \code{fixationSites} uses the result to find those sites that show
#'   fixation on some, if not all, of the lineages. Parallel evolution is
#'   relatively common in RNA virus. There is chance that some site be fixed in
#'   one lineage but does not show fixation because of different sequence
#'   context.
#' @param paths a \code{lineagePath} object returned from
#'   \code{\link{lineagePath}} function or a \code{phylo} object after
#'   \code{\link{addMSA}}
#' @param minEffectiveSize A vector of two integers to specifiy minimum tree
#'   tips involved before/after mutation. Otherwise the mutation will not be
#'   counted into the return. If more than one number is given, the ancestral
#'   takes the first and descendant takes the second as the minimum. If only
#'   given one number, it's the minimum for both ancestral and descendant.
#' @param searchDepth The function uses heuristic search but the termination of
#'   the search cannot be intrinsically decided. \code{searchDepth} is needed to
#'   tell the search when to stop.
#' @param method The strategy for predicting the fixation. The basic approach is
#'   entropy minimization and can be achieved by adding or removing fixation
#'   point, or by comparing the two.
#' @param ... further arguments passed to or from other methods.
#' @return \code{fixationSites} returns a list of fixation mutations with names
#'   of the tips involved.
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
    nodeAlign <- .tipSeqsAlongPathNodes(
        paths = paths,
        divNodes = divNodes,
        tree = tree,
        align = align
    )
    res <- .findFixationSite(
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
    groupingByPath <- .transitionClusters(res, paths, tree)
    res <- .combineFixations(res, tree, align)
    attr(res, "paths") <- paths
    attr(res, "reference") <- reference
    attr(res, "clustersByPath") <- groupingByPath
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
                    attr(res[[site]], "tree") <- tree
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
    # Find the divergent point and remove the overlapped part
    grouping <- list(groupByPath[[1]])

    for (gpIndex in seq_along(groupByPath)[-1]) {
        gp <- groupByPath[[gpIndex]]
        toMergeIndex <- NULL
        divergedIndex <- 0L
        # Loop through to find the most related group
        for (i in seq_along(grouping)) {
            allTips <- unlist(grouping[[i]])
            for (j in seq_along(gp)[-1]) {
                if (all(!gp[[j]] %in% allTips)) {
                    m <- i
                    d <- j
                    break
                }
            }
            if (d > divergedIndex) {
                toMergeIndex <- m
                divergedIndex <- d
            }
        }
        # Find the tips before diverged
        sharedTips <- gp[[divergedIndex - 1]]
        refSites <- attr(sharedTips, "site")
        # The non-shared part
        divergedTips <- setdiff(sharedTips, allTips)
        attr(divergedTips, "site") <- refSites
        # Drop the overlapped part
        grouping[[gpIndex]] <-
            c(list(divergedTips), gp[divergedIndex:length(gp)])
        # Find the most related group
        toMerge <- grouping[[toMergeIndex]]
        # To determine where to add the new (truncated) group
        for (i in seq_along(toMerge)) {
            gpTips <- unlist(gp)
            if (all(!toMerge[[i]] %in% gpTips)) {
                # Find the tips before diverged
                sharedTips <- toMerge[[i - 1]]
                sites <- attr(sharedTips, "site")
                # The non-shared part
                divergedTips <- setdiff(sharedTips, gpTips)
                attr(divergedTips, "site") <- sites
                # The shared part
                sharedTips <- setdiff(sharedTips, divergedTips)
                attr(sharedTips, "site") <- sites
                attr(sharedTips, "toMerge") <- gpIndex
                attr(sharedTips, "toMergeRefSites") <- refSites
                # Reform
                if (i == 2) {
                    preTips <- list()
                } else {
                    preTips <- toMerge[seq_len(i - 2)]
                }
                grouping[[toMergeIndex]] <- c(preTips,
                                              list(sharedTips),
                                              list(divergedTips),
                                              toMerge[i:length(toMerge)])
                break
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

fixationSites.phylo <- function(paths, ...) {
    align <- attr(paths, "align")
    # Generate the site mapping from reference
    reference <- attr(paths, "reference")
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    res <- fixationSitesSearch(nodepath(paths), align, loci)
    res <- res[which(lengths(res) != 1)]
    return(res)
}

treemerBySite <- function(x, ...) {
    align <- attr(x, "align")
    # Generate the site mapping from reference
    reference <- attr(x, "reference")
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    res <- runTreemerBySite(nodepath(x), align, loci)
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
#' x <- sitewiseClusters(fixations)
#' plot(x)
plot.fixationSites <- function(x,
                               y = TRUE,
                               ...) {
    tree <- as.phylo.fixationSites(x)
    grp <- sitewiseClusters.fixationSites(x, minEffectiveSize = 0)
    grp <- as.list.sitewiseClusters(grp)
    clusterPaths <- list()
    rootNode <- getMRCA(tree, tree[["tip.label"]])
    for (cluster in names(grp)) {
        tips <- grp[[cluster]]
        ancestral <- getMRCA(tree, tips)
        if (is.null(ancestral)) {
            np <- nodepath(tree, rootNode, tips)
            clusterPaths[[cluster]] <- np[1:(length(np) - 1)]
        } else {
            clusterPaths[[cluster]] <- nodepath(tree, rootNode, ancestral)
        }
    }
    clusterInfo <- lapply(names(grp), function(g) {
        data.frame(row.names = grp[[g]],
                   "cluster" = rep(g, length(grp[[g]])))
    })
    clusterInfo <- do.call(rbind, clusterInfo)

    transMut <- list()
    for (sp in x) {
        site <- attr(sp, "site")
        for (mp in sp) {
            for (i in seq_along(mp)[-1]) {
                prevTips <- mp[[i - 1]]
                currTips <- mp[[i]]
                prevAA <- attr(prevTips, "AA")
                currAA <- attr(currTips, "AA")
                mutation <- paste0(prevAA, site, currAA)
                currCluster <-
                    unique(clusterInfo[as.character(currTips), ])
                names(currCluster) <- currCluster
                # Choose the most ancient cluster which first receive the
                # mutation
                curr <- names(which.min(lapply(
                    X = currCluster,
                    FUN = function(cluster) {
                        length(clusterPaths[[cluster]])
                    }
                )))
                # Find the transition node
                currPath <- clusterPaths[[curr]]
                trans <- as.character(currPath[length(currPath)])
                if (trans %in% names(transMut)) {
                    transMut[[trans]] <- c(transMut[[trans]], mutation)
                } else {
                    transMut[[trans]] <- mutation
                }
            }
        }
    }
    transMut <- lapply(transMut, unique)
    d <- as_tibble(t(vapply(
        X = names(transMut),
        FUN = function(trans) {
            snp <- transMut[[trans]]
            res <- character()
            snpNum <- length(snp)
            for (i in seq_len(snpNum)) {
                res <- paste0(res, snp[i])
                if (i < snpNum) {
                    if (i %% 4 == 0) {
                        res <- paste0(res, ",\n")
                    } else {
                        res <- paste0(res, ", ")
                    }
                }
            }
            res <- c(trans, res)
            names(res) <- c("node", "SNPs")
            return(res)
        },
        FUN.VALUE = character(2)
    )))
    d[["node"]] <- as.integer(d[["node"]])
    d <- full_join(as_tibble(tree), d, by = "node")
    tree <- as.treedata(d)

    groupColors <- colorRampPalette(brewer.pal(9, "Set1"))(length(grp))
    names(groupColors) <- names(grp)
    groupColors["0"] <- "black"

    p <- ggtree(groupOTU(tree, grp, group_name = "Groups"),
                aes(color = Groups)) +
        geom_label_repel(
            aes(x = branch, label = SNPs),
            fill = "lightgreen",
            color = "black",
            min.segment.length = 0,
            na.rm = TRUE
        ) +
        scale_color_manual(values = groupColors)
    if (y) {
        p <- p +
            guides(color = guide_legend(override.aes = list(size = 3), )) +
            theme(legend.position = "left")
    } else {
        p <- p + theme(legend.position = "none")
    }
    return(p)
}
