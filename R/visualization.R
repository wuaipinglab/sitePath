#' @name visualization
#' @title Visualize results
#' @description Visualize \code{\link{lineagePath}} object. A tree diagram will
#'   be plotted and paths are black solid line while the trimmed nodes and tips
#'   will use grey dashed line.
#' @param x Could be a \code{\link{lineagePath}} object, a
#'   \code{\link{fixationSites}} object or a \code{sitePath} object.
#' @param y For \code{\link{lineagePath}} object, it is deprecated. For a
#'   \code{\link{fixationSites}} object, it is whether to show the fixation
#'   mutation between clusters. For a \code{sitePath} object, it can have more
#'   than one fixation path. This is to select which path to plot. The default
#'   is \code{NULL} which will plot all the paths.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @param ... Arguments in \code{plot.phylo} functions.
#' @return The function only makes plot and returns no value (It behaviors like
#'   the generic \code{\link{plot}} function).
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' paths <- lineagePath(tree)
#' plot(paths)
#' @importFrom graphics plot
#' @importFrom ape plot.phylo
#' @export
plot.lineagePath <- function(x, y = TRUE, showTips = FALSE, ...) {
    tree <- attr(x, "tree")
    tree <- ladderize(tree, right = FALSE)
    nEdges <- length(tree$edge.length)
    color <- rep("#d3d3d3", nEdges)
    lty <- rep(2, nEdges)
    width <- rep(1, nEdges)
    targetEdges <- which(tree$edge[, 2] %in% unique(unlist(x)))
    color[targetEdges] <- "#000000"
    lty[targetEdges] <- 1
    width[targetEdges] <- 2
    # TODO: Emphaszie the nodes along the lineagePath
    show.tip.label <- showTips
    plot.phylo(
        tree,
        edge.color = color,
        edge.lty = lty,
        edge.width = width,
        show.tip.label = show.tip.label,
        ...
    )
}

.transitionClusters <- function(fixations) {
    paths <- attr(fixations, "paths")
    tree <- attr(paths, "tree")
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
                    if (length(common) == 0) {
                        newGrouping <- res[1:i]
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
                            newGrouping <-
                                c(newGrouping,
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
                    preTips <- toMerge[1:(i - 2)]
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

.snpTracing <- function(grouping, minEffectiveSize) {
    grouping <- lapply(grouping, function(g) {
        g[which(vapply(
            X = g,
            FUN = function(tips) {
                length(tips) > minEffectiveSize || !is.null(attr(tips, "toMerge"))
            },
            FUN.VALUE = logical(1)
        ))]
    })
    # The existing tips in the tree
    tipClusters <- list()
    # The edges and its SNP of the tree
    parentNodes <- integer()
    childrenNodes <- integer()
    edgeSNPs <- list()

    # Keep track of the newly added internal node
    newParentNode <-
        length(unlist(grouping, recursive = FALSE)) + 1L
    parentNode <- newParentNode
    # Keep track of the newly added tip node
    tipNode <- 0L
    # A list to record the fixation sites of the parent nodes
    parentNodesSites <- list()

    for (gpIndex in seq_along(grouping)) {
        # The group to add onto the tree
        gp <- grouping[[gpIndex]]

        # The initial tips of the group
        currentTips <- gp[[1]]
        # Assume the initial reference site
        refSites <- attr(currentTips, "site")
        # Find where to merge and parent node, update reference site maybe
        for (i in seq_along(tipClusters)) {
            toMerge <- tipClusters[[i]]
            toMergeIndex <- attr(toMerge, "toMerge")
            if (!is.null(toMergeIndex) && toMergeIndex == gpIndex) {
                refSites <- attr(toMerge, "toMergeRefSites")
                parentNode <- parentNodes[which(childrenNodes == i)]
                break
            }
        }
        parentNodesSites[[as.character(parentNode)]] <- refSites

        # Track the tip and internal node
        tipNode <- tipNode + 1L
        startingNode <- tipNode
        tipClusters <- c(tipClusters, list(currentTips))
        # Define the initial edge of the group
        parentNodes <- c(parentNodes, parentNode)
        childrenNodes <- c(childrenNodes, tipNode)
        # SNP of the initial edge is set none
        edgeSNPs <- c(edgeSNPs, list(as.character(na.omit(
            vapply(
                X = names(refSites),
                FUN = function(site) {
                    ref <- refSites[site]
                    snp <- attr(currentTips, "site")[site]
                    if (ref == snp) {
                        return(NA_character_)
                    }
                    return(paste0(ref, site, snp))
                },
                FUN.VALUE = character(1)
            )
        ))))

        # Grow the tree
        for (tipIndex in seq_along(gp)[-1]) {
            tipNode <- tipNode + 1L
            currentTips <- gp[[tipIndex]]
            currentSites <- attr(currentTips, "site")
            # Attach the tip near the most related tips. Assume the reference
            # tips are the most related (least number of SNP)
            mostRelatedTipNode <- startingNode
            leastSNPnum <- sum(refSites != currentSites)
            # Loop through the rest existing tip clusters
            for (otherTipNode in seq_along(tipClusters)[-seq_len(startingNode)]) {
                otherSites <- attr(tipClusters[[otherTipNode]], "site")
                snpNum <- sum(otherSites != currentSites)
                if (snpNum < leastSNPnum) {
                    mostRelatedTipNode <- otherTipNode
                    leastSNPnum <- snpNum
                }
            }
            # Find the direct tree edge to the most related tips
            edgeIndex <- which(childrenNodes == mostRelatedTipNode)
            parentNode <- parentNodes[edgeIndex]
            # Tree growing differs according to the edge SNP
            parentSites <-
                parentNodesSites[[as.character(parentNode)]]
            snpSites <- as.character(na.omit(
                vapply(
                    X = names(parentSites),
                    FUN = function(site) {
                        ref <- parentSites[site]
                        snp <- currentSites[site]
                        if (ref == snp) {
                            return(NA_character_)
                        }
                        return(paste0(ref, site, snp))
                    },
                    FUN.VALUE = character(1)
                )
            ))
            edgeSNP <- edgeSNPs[[edgeIndex]]
            sharedWithEdgeSNP <- intersect(snpSites, edgeSNP)
            # A new internal node is needed when no SNP overlap
            if (length(sharedWithEdgeSNP) != 0) {
                newParentNode <- newParentNode + 1L
                # Insert the new internal node to the target edge
                parentNodes[edgeIndex] <- newParentNode
                parentNodes <- c(parentNodes, parentNode)
                childrenNodes <- c(childrenNodes, newParentNode)
                edgeSNPs <- c(edgeSNPs, list(sharedWithEdgeSNP))
                # Update SNP of the directly linked edge
                edgeSNPs[[edgeIndex]] <-
                    setdiff(edgeSNP, sharedWithEdgeSNP)
                # Calculate the site for the new internal node
                siteToChange <- substr(sharedWithEdgeSNP,
                                       2,
                                       nchar(sharedWithEdgeSNP) - 1)
                parentSites[siteToChange] <- substr(
                    sharedWithEdgeSNP,
                    nchar(sharedWithEdgeSNP),
                    nchar(sharedWithEdgeSNP)
                )
                parentNodesSites[[as.character(newParentNode)]] <-
                    parentSites
                # Update the parent node and edge SNP for the current tip node
                parentNode <- newParentNode
                snpSites <- setdiff(snpSites, sharedWithEdgeSNP)
            }
            # Add edge
            parentNodes <- c(parentNodes, parentNode)
            childrenNodes <- c(childrenNodes, tipNode)
            # Add edge SNP
            edgeSNPs <- c(edgeSNPs, list(snpSites))
            # Add the current tips
            tipClusters <- c(tipClusters, list(currentTips))
        }
    }
    res <- list(
        "edge" = cbind(parentNodes, childrenNodes),
        "edge.length" = lengths(edgeSNPs),
        "Nnode" = length(unique(parentNodes)),
        "tip.label" = as.character(lengths(unlist(
            grouping, recursive = FALSE
        )))
    )
    attr(res, "tipClusters") <- tipClusters
    attr(res, "edgeSNPs") <- edgeSNPs
    class(res) <- "phylo"
    return(res)
}

#' @name visualization
#' @description Visualize \code{\link{fixationSites}} object. The tips are
#'   clustered according to the fixation sites. The transition of fixation sites
#'   will be plotted as a phylogenetic tree. The length of each branch
#'   represents the number of fixation mutation between two clusters. The name
#'   of the tree tips indicate the number of sequences in the cluster.
#' @param recurringOnly Whether to plot recurring fixation mutation only. The
#'   default is FALSE.
#' @param minEffectiveSize The minimum size for a tip cluster in the plot
#' @examples
#' fixations <- fixationSites(paths)
#' plot(fixations)
#' @importFrom ape edgelabels
#' @importFrom ape axisPhylo
#' @export
plot.fixationSites <- function(x,
                               y = TRUE,
                               showTips = FALSE,
                               recurringOnly = FALSE,
                               minEffectiveSize = NULL,
                               ...) {
    grouping <- .transitionClusters(x)
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- mean(lengths(unlist(grouping, recursive = FALSE)))
    }
    snpTracing <- .snpTracing(grouping, minEffectiveSize)
    edgeSNPs <- attr(snpTracing, "edgeSNPs")
    if (recurringOnly) {
        allMutSites <- unlist(edgeSNPs)
        duplicatedSites <-
            unique(allMutSites[which(duplicated(allMutSites))])
        edgeSNPs <- lapply(edgeSNPs, function(sites) {
            sites[which(sites %in% duplicatedSites)]
        })
    }
    edge2show <- which(lengths(edgeSNPs) != 0)
    show.tip.label <- showTips
    plot.phylo(snpTracing, show.tip.label = show.tip.label, ...)
    axisPhylo(backward = FALSE)
    if (y) {
        edgelabels(
            text = vapply(
                edgeSNPs[edge2show],
                paste,
                collapse = ", ",
                FUN.VALUE = character(1)
            ),
            edge = edge2show
        )
    }
}

#' @name visualization
#' @description Visualize the \code{sitePath} object which can be extracted by
#'   using \code{\link{extractSite}} on the return of
#'   \code{\link{fixationSites}} and \code{\link{multiFixationSites}}.
#' @examples
#' sp <- extractSite(fixations, 139)
#' plot(sp)
#' @importFrom graphics title
#' @importFrom graphics legend
#' @seealso \code{\link{plotSingleSite}}, \code{\link{extractSite}}
#' @export
plot.sitePath <- function(x, y = NULL, showTips = FALSE, ...) {
    tree <- attr(x, "tree")
    # Prepare tree for plotting
    tree <- ladderize(tree, right = FALSE)
    rootNode <- getMRCA(tree, tree$tip.label)
    plotName <- character(0)
    nEdges <- length(tree$edge.length)
    color <- rep("#d3d3d3", nEdges)
    lty <- rep(2, nEdges)
    width <- rep(0.5, nEdges)
    AAnames <- character(0)
    if (is.null(y)) {
        sitePaths <- x[]
    } else {
        tryCatch(
            expr = sitePaths <- x[y],
            error = function(e) {
                stop("There are ",
                     length(x),
                     " in \"x\". ",
                     "The selection \"y\" is out of bounds.")
            }
        )
    }
    for (sp in sitePaths) {
        aaName <- character(0)
        for (tips in rev(sp)) {
            aa <- AA_FULL_NAMES[tolower(attr(tips, "AA"))]
            aaName <- c(aa, aaName)
            targetEdges <- tip2Edge(tree$edge, tips, rootNode)
            color[targetEdges] <- AA_COLORS[aa]
            lty[targetEdges] <- 1
            width[targetEdges] <- 2
        }
        AAnames <- c(AAnames, aaName)
        plotName <-
            c(plotName, paste0(AA_SHORT_NAMES[aaName], collapse = " -> "))
    }
    show.tip.label <- showTips
    plot.phylo(
        tree,
        show.tip.label = show.tip.label,
        edge.color = color,
        edge.lty = lty,
        edge.width = width,
        ...
    )
    sepChar <- "\n"
    if (sum(nchar(plotName) <= 18)) {
        sepChar <- ", "
    }
    title(main = attr(x, "site"),
          sub = paste(plotName, collapse = sepChar))
    legend(
        "topleft",
        title = "Amino acid",
        legend = AA_SHORT_NAMES[unique(AAnames)],
        fill = AA_COLORS[unique(AAnames)],
        box.lty = 0
    )
}

#' @rdname plotSingleSite
#' @name plotSingleSite
#' @title Color the tree by a single site
#' @description For \code{lineagePath}, the tree will be colored according to
#'   the amino acid of the site. The color scheme tries to assign
#'   distinguishable color for each amino acid.
#' @param x A \code{fixationSites} object from \code{\link{fixationSites}} or
#'   the return from \code{\link{addMSA}} function.
#' @param site One of the mutations in the \code{fixationSites} object. It
#'   should be from the \code{\link{names}} of the object. Or an integer to
#'   indicate a site could be provide. The numbering is consistent with the
#'   reference defined at \code{\link{fixationSites}}.
#' @param showPath If plot the lineage result from lineagePath.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @param ... Arguments in \code{plot.phylo} functions and other arguments.
#' @return The function only makes plot and returns no value (It behaviors like
#'   the generic \code{\link{plot}} function).
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' paths <- lineagePath(tree)
#' plotSingleSite(paths, 139)
#' @seealso \code{\link{plot.sitePath}}
#' @importFrom ape ladderize
#' @importFrom ape getMRCA
#' @export
plotSingleSite.lineagePath <- function(x,
                                       site,
                                       showPath = FALSE,
                                       showTips = FALSE,
                                       ...) {
    site <- .checkSite(site)
    tree <- attr(x, "tree")
    tree <- ladderize(tree, right = FALSE)
    align <- attr(x, "align")
    align <- strsplit(tolower(align), "")
    reference <- attr(x, "reference")
    tryCatch(
        expr = site <- match.arg(as.character(site), seq_along(reference)),
        error = function(e) {
            stop("\"site\": ", site, " is not within the length of reference")
        }
    )
    siteComp <- vapply(align,
                       FUN = "[[",
                       FUN.VALUE = character(1),
                       reference[site])
    nEdges <- length(tree$edge.length)
    color <- rep("#000000", nEdges)
    rootNode <- getMRCA(tree, tree$tip.label)
    group <- list()
    for (i in seq_along(siteComp)) {
        group[[siteComp[[i]]]] <- c(group[[siteComp[[i]]]], i)
    }
    AAnames <- AA_FULL_NAMES[names(group)]
    names(group) <- AA_COLORS[AAnames]
    for (g in names(group)) {
        tip2colorEdge(color, g, tree$edge, group[[g]], rootNode)
    }
    width <- rep(1, nEdges)
    if (showPath) {
        targetEdges <- which(tree$edge[, 2] %in% unique(unlist(x)))
        color[targetEdges] <- "#000000"
        width[targetEdges] <- 2
    }
    show.tip.label <- showTips
    plot.phylo(
        tree,
        show.tip.label = show.tip.label,
        edge.color = color,
        edge.width = width,
        main = site,
        ...
    )
    legend(
        "topleft",
        title = "Amino acid",
        legend = unique(AAnames),
        fill = AA_COLORS[unique(AAnames)],
        box.lty = 0
    )
}

#' @rdname plotSingleSite
#' @description For \code{fixationSites}, it will color the ancestral tips in
#'   red, descendant tips in blue and excluded tips in grey.
#' @param select Select which fixation path in to plot. The default is NULL
#'   which will plot all the fixations.
#' @examples
#' fixations <- fixationSites(paths)
#' plotSingleSite(fixations, 139)
#' @export
plotSingleSite.fixationSites <- function(x,
                                         site,
                                         select = NULL,
                                         ...) {
    site <- .checkSite(site)
    tryCatch(
        expr = site <- match.arg(as.character(site), choices = names(x)),
        error = function(e) {
            stop("\"site\": ", site, " is not a mutation of fixation")
        }
    )
    plot.sitePath(x = x[[site]], y = select, ...)
}

#' @rdname plotSingleSite
#' @description For \code{multiFixationSites}, it will color the tips which have
#'   their site fixed. The color will use the same amino acid color scheme as
#'   \code{plotSingleSite.lineagePath}
#' @examples
#' \dontrun{
#' multiFixations <- multiFixationSites(paths)
#' plotSingleSite(multiFixations, 1542)
#' }
#' @export
plotSingleSite.multiFixationSites <- function(x,
                                              site,
                                              select = NULL,
                                              ...) {
    site <- .checkSite(site)
    tryCatch(
        expr = site <- match.arg(as.character(site), choices = names(x)),
        error = function(e) {
            stop("\"site\": ", site, " is not a mutation of fixation")
        }
    )
    plot.sitePath(x = x[[site]], y = select, ...)
}

#' @export
plotSingleSite <- function(x, ...)
    UseMethod("plotSingleSite")

.checkSite <- function(site) {
    if (!is.numeric(site) ||
        any(site <= 0) || as.integer(site) != site) {
        stop("Please enter positive integer value for \"site\"")
    }
    if (length(site) != 1) {
        warning("\"site\" has more than one element, only the first ",
                site[1],
                " will be used.")
    }
    return(site[1])
}
