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
#' @importFrom ggtree ggtree aes theme scale_color_manual geom_tiplab
#' @export
plot.lineagePath <- function(x, y = TRUE, showTips = FALSE, ...) {
    tree <- attr(x, "tree")
    nNodes <- length(tree[["tip.label"]]) + tree[["Nnode"]]
    attr(tree, "group") <- rep(1, times = nNodes)
    attr(tree, "group")[unique(unlist(paths))] <- 0
    attr(tree, "group") <- factor(attr(tree, "group"))

    attr(tree, "size") <- rep(0.5, times = nNodes)
    attr(tree, "size")[unique(unlist(paths))] <- 1

    p <- ggtree(tree, aes(color = group, linetype = group, size = size)) +
        scale_color_manual(values = c("black", "grey")) +
        theme(legend.position = "none")
    if (showTips) {
        p <- p + geom_tiplab()
    }
    return(p)
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
#' @importFrom tidytree as_tibble
#' @importFrom ape edgelabels
#' @importFrom ape axisPhylo
#' @seealso \code{\link{as.phylo.fixationSites}}
#' @export
plot.fixationSites <- function(x,
                               y = TRUE,
                               showTips = FALSE,
                               recurringOnly = FALSE,
                               minEffectiveSize = NULL,
                               ...) {
    snpTracing <- as.phylo.fixationSites(x, minEffectiveSize)
    edgeSNPs <- attr(snpTracing, "edgeSNPs")
    if (recurringOnly) {
        allMutSites <- unlist(edgeSNPs)
        duplicatedSites <-
            unique(allMutSites[which(duplicated(allMutSites))])
        edgeSNPs <- lapply(edgeSNPs, function(sites) {
            res <- sites[which(sites %in% duplicatedSites)]
            attributes(res) <- attributes(sites)
            return(res)
        })
    }
    edge2show <- which(lengths(edgeSNPs) != 0)
    show.tip.label <- showTips
    plot.phylo(snpTracing, show.tip.label = show.tip.label, ...)
    axisPhylo(backward = FALSE)
    if (y) {
        edgelabels(
            text = vapply(
                X = edgeSNPs[edge2show],
                FUN = paste,
                collapse = ", ",
                FUN.VALUE = character(1)
            ),
            edge = vapply(
                X = edgeSNPs[edge2show],
                FUN = function(i) {
                    which(snpTracing[["edge"]][, 2] == attr(i, "edge")[2])
                },
                FUN.VALUE = integer(1)
            )
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
#' @importFrom graphics plot
#' @importFrom graphics title
#' @importFrom graphics legend
#' @importFrom ape plot.phylo
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
