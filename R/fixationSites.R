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
    # Find fixation sites
    res <- .findFixationSite(
        paths = paths,
        minEffectiveSize = minEffectiveSize,
        searchDepth = searchDepth,
        method = method
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

.findFixationSite <- function(paths,
                              minEffectiveSize,
                              searchDepth,
                              method = c("compare", "insert", "delete")) {
    # Get the divergent nodes
    divNodes <- divergentNode(paths)
    # The tips and matching
    nodeAlign <- .tipSeqsAlongPathNodes(paths = paths,
                                        divNodes = divNodes)
    # Decide which minimizing strategy
    minimizeEntropy <- switch(
        match.arg(method),
        "compare" = minEntropyByComparing,
        "insert" = minEntropyByInserting,
        "delete" = minEntropyByDeleting
    )
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
    # In case root node does not have any tips (because itself is a divergent
    # node)
    rootNode <- attr(paths, "rootNode")
    if (!rootNode %in% names(nodeAlign)) {
        paths <- lapply(paths, "[", -1)
    }
    # Iterate each path
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
    res <- lapply(paths, function(p) {
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
    res <- res[which(lengths(res) != 0)]
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
    cat("This is a 'fixationSites' object.\n\nResult for",
        length(attr(x, "paths")),
        "paths:\n\n")
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
