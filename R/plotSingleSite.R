#' @importFrom ggplot2 aes ggtitle guides theme guide_legend
#' @importFrom ggplot2 scale_size scale_color_manual
#' @importFrom ggplot2 GeomSegment GeomText GeomPoint
#' @importFrom ggrepel geom_label_repel
#' @importFrom tidytree groupOTU
#' @importFrom ggtree ggtree geom_tippoint geom_tiplab

#' @rdname plotSingleSite
#' @title Color the tree by a single site
#' @description Plot and color the tree according to amino acid/nucleotide of
#'   the selected site. The color scheme depends on the \code{seqType} set in
#'   \code{\link{addMSA}} function.
#' @param x A \code{fixationSites} object from \code{\link{fixationSites}} or
#'   the return from \code{\link{lineagePath}} function.
#' @param site For \code{lineagePath}, it can be any site within sequence
#'   length. For \code{fixationSites} and \code{parallelSites}, it is restrained
#'   to a predicted fixation site. The numbering is consistent with the
#'   reference defined by \code{\link{setSiteNumbering}}.
#' @param ...  Other arguments. Since 1.5.4, the function uses
#'   \code{\link{ggtree}} as the base function to make plots so the arguments in
#'   \code{plot.phylo} will no longer work.
#' @return Since 1.5.4, the function returns a ggplot object so on longer
#'   behaviors like the generic \code{\link{plot}} function.
#' @seealso \code{\link{plot.sitePath}}
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' paths <- lineagePath(tree)
#' plotSingleSite(paths, 139)
plotSingleSite <- function(x, site, ...)
    UseMethod("plotSingleSite")

#' @rdname plotSingleSite
#' @description For \code{\link{lineagePath}}, the tree will be colored
#'   according to the amino acid of the site. The color scheme tries to assign
#'   distinguishable color for each amino acid.
#' @param showPath If plot the lineage result from \code{\link{lineagePath}}.
#'   The default is \code{TRUE}.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @export
plotSingleSite.lineagePath <- function(x,
                                       site,
                                       showPath = TRUE,
                                       showTips = FALSE,
                                       ...) {
    group <- extractTips.lineagePath(x, site)
    # Use different color scheme depending on the sequence type
    names(group) <- toupper(names(group))
    groupColors <- .siteColorScheme(attr(x, "seqType"))
    tree <- attr(x, "tree")
    group <- groupOTU(as_tibble(tree), group)
    group <- group[["group"]]
    size <- NULL
    sizeRange <- c(GeomSegment[["default_aes"]][["size"]], 1.5)
    # Set lineage nodes and non-lineage nodes as separate group
    if (showPath) {
        pathNodes <- unique(unlist(x))
        pathLabel <- ".lineage"
        # Color the path node black
        levels(group) <- c(levels(group), pathLabel)
        group[pathNodes] <- pathLabel
        lineageColor <- "black"
        names(lineageColor) <- pathLabel
        groupColors <- c(groupColors, lineageColor)
        # Set the size of the lineage nodes
        size <- rep(1, times = length(group))
        size[pathNodes] <- 2
    }
    p <- ggtree(tree, aes(color = group, size = size)) +
        scale_size(range = sizeRange, guide = FALSE) +
        scale_color_manual(values = groupColors) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme(legend.position = "left") +
        ggtitle(site)
    if (showTips) {
        p <- p + geom_tiplab()
    }
    return(p)
}

.siteColorScheme <- function(seqType) {
    if (seqType == "AA") {
        groupColors <- vapply(
            X = AA_FULL_NAMES,
            FUN = function(i) {
                AA_COLORS[[i]]
            },
            FUN.VALUE = character(1)
        )
    } else {
        groupColors <- NT_COLORS
    }
    names(groupColors) <- toupper(names(groupColors))
    return(groupColors)
}

#' @rdname plotSingleSite
#' @description For \code{\link{parallelSites}}, the tree will be colored
#'   according to the amino acid of the site if the mutation is not fixed.
#' @export
plotSingleSite.parallelSites <- function(x,
                                         site,
                                         showPath = TRUE,
                                         ...) {
    paths <- attr(x, "paths")
    tree <- attr(paths, "tree")
    tipNames <- tree[["tip.label"]]
    nNodes <- length(tipNames) + tree[["Nnode"]]
    parallelMut <- extractTips.parallelSites(x, site)
    fixationMut <- character()
    sporadicTip <- rep(FALSE, nNodes)
    for (node in names(parallelMut)) {
        tips <- parallelMut[[node]]
        if (attr(tips, "fixed")) {
            fixationMut[node] <- attr(tips, "mutName")[4]
        } else {
            sporadicTip[which(tipNames == node)] <- TRUE
        }
    }
    if (length(fixationMut) != 0) {
        attr(paths, "tree") <- .annotateSNPonTree(tree, fixationMut)
        p <- plotSingleSite.lineagePath(
            x = paths,
            site = site,
            showPath = showPath,
            showTips = FALSE
        ) +
            geom_label_repel(
                aes(x = branch, label = SNPs),
                fill = 'lightgreen',
                color = "black",
                min.segment.length = 0,
                na.rm = TRUE,
                size = GeomText[["default_aes"]][["size"]]
            )
    } else {
        p <- plotSingleSite.lineagePath(
            x = paths,
            site = site,
            showPath = showPath,
            showTips = FALSE
        )
    }
    if (any(sporadicTip)) {
        p <- p + geom_tippoint(aes(subset = sporadicTip,
                                   size = GeomPoint[["default_aes"]][["size"]]))
    }
    return(p)
}

#' @rdname plotSingleSite
#' @description For \code{\link{fixationSites}}, it will color the ancestral
#'   tips in red, descendant tips in blue and excluded tips in grey.
#' @param select Select which fixation path in to plot. The default is NULL
#'   which will plot all the fixations.
#' @export
#' @examples
#' fixations <- fixationSites(paths)
#' plotSingleSite(fixations, 139)
plotSingleSite.fixationSites <- function(x,
                                         site,
                                         select = NULL,
                                         ...) {
    plot.sitePath(x = .actualExtractSite(x, site), y = select, ...)
}
