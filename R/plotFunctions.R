#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 aes guides theme ggtitle guide_legend
#' @importFrom ggplot2 scale_size scale_color_manual
#' @importFrom tidytree groupOTU
#' @importFrom ggtree ggtree geom_tiplab theme_tree2
#' @importFrom ggrepel geom_label_repel

plot.phyMSAmatched <- function(x, y = TRUE) {
    p <- ggtree(as.phylo.phyMSAmatched(x))
    return(p)
}

#' @rdname plotFunctions
#' @title Visualize the results
#' @description The plot function to visualize the return of functions in the
#'   package. The underlying function applies \code{\link{ggplot2}}. The
#'   function name \code{plot} is used to keep the compatibility with previous
#'   versions, but they do not behave like the generic \code{\link{plot}}
#'   function since 1.5.4.
#' @param x Could be a \code{\link{lineagePath}} object,
#'   \code{\link{fixationSites}} object or \code{\link{fixationPath}} object.
#' @param y For \code{\link{lineagePath}} object, it is deprecated and no longer
#'   have effect since 1.5.4. For a \code{sitePath} object, it can have more
#'   than one fixation path. This is to select which path to plot. The default
#'   is \code{NULL} which will plot all the paths. It is the same as
#'   \code{select} in \code{plotSingleSite} For a \code{\link{fixationSites}} or
#'   a \code{\link{fixationPath}} object, it is whether to show the fixation
#'   mutation between clusters.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @param ... Other arguments. Since 1.5.4, the function uses
#'   \code{\link{ggtree}} as the base function to make plots so the arguments in
#'   \code{plot.phylo} will no longer work.
#' @return A ggplot object to make the plot. A \code{\link{lineagePath}} object
#'   will be plotted as a tree diagram will be plotted and paths are black solid
#'   line while the trimmed nodes and tips will use grey dashed line. A
#'   \code{\link{fixationSites}} object will be plotted as original phylogenetic
#'   tree marked with fixation substitutions. A \code{\link{fixationPath}}
#'   object will be plotted as a \code{phylo} object. The tips are clustered
#'   according to the fixation sites. The transition of fixation sites will be
#'   plotted as a phylogenetic tree. The length of each branch represents the
#'   number of fixation mutation between two clusters. The name of the tree tips
#'   indicate the number of sequences in the cluster.
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' paths <- lineagePath(tree)
#' plot(paths)
plot.lineagePath <- function(x,
                             y = TRUE,
                             showTips = FALSE,
                             ...) {
    tree <- attr(x, "tree")
    # Get number of ancestral nodes plus tip nodes
    nNodes <- length(tree[["tip.label"]]) + tree[["Nnode"]]
    # Set lineage nodes and non-lineage nodes as separate group
    group <- rep(1, times = nNodes)
    group[unique(unlist(x))] <- 0
    group <- factor(group)
    # Set line size
    size <- rep(1, times = nNodes)
    size[unique(unlist(x))] <- 2
    # Tree plot
    p <- ggtree(tree, aes(
        color = group,
        linetype = group,
        size = size
    )) +
        scale_size(range = c(GeomSegment[["default_aes"]][["size"]], 1.5)) +
        scale_color_manual(values = c("black", "gainsboro")) +
        theme(legend.position = "none")
    if (showTips) {
        p <- p + geom_tiplab()
    }
    return(p)
}

#' @rdname plotFunctions
#' @examples
#' fixations <- fixationSites(paths)
#' plot(fixations)
#' @export
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

#' @rdname plotFunctions
#' @description Visualize the \code{sitePath} object which can be extracted by
#'   using \code{\link{extractSite}} on the return of
#'   \code{\link{fixationSites}}.
#' @export
#' @examples
#' sp <- extractSite(fixations, 139)
#' plot(sp)
plot.sitePath <- function(x, y = NULL, showTips = FALSE, ...) {
    tree <- attr(x, "tree")
    if (is.null(y)) {
        sitePaths <- x[]
    } else {
        if (length(x) < y) {
            stop(
                "There are ",
                length(x),
                "lineages with fixation mutation. ",
                "The selection \"y\" is out of bounds."
            )
        }
        sitePaths <- x[y]
    }
    # Specify the color of mutations by pre-defined color set.
    groupColors <- .siteColorScheme(attr(x, "seqType"))
    # Collect the fixation mutation for each evolutionary pathway
    subtitle <- character()
    group <- list()
    for (sp in sitePaths) {
        aaName <- character()
        for (tips in sp) {
            aa <- toupper(attr(tips, "AA"))
            aaName <- c(aaName, aa)
            group[[aa]] <- c(group[[aa]], tips)
        }
        subtitle <- c(subtitle,
                      paste0(aaName, collapse = " -> "))
    }
    tree <- groupOTU(tree, group)
    # Color the excluded branch gray as the excluded lineage
    excludedLabel <- ".excluded"
    groupColors[[excludedLabel]] <- "#dcdcdc"
    levels(attr(tree, "group"))[which(levels(attr(tree, "group")) == "0")] <-
        excludedLabel
    # Use dashed line for excluded branches
    linetype <- as.integer(attr(tree, "group") == excludedLabel)
    linetype <- factor(linetype)
    # Just in case the fixation mutation name is too long
    sepChar <- "\n"
    if (sum(nchar(subtitle) <= 18)) {
        sepChar <- ", "
    }
    subtitle <- paste(subtitle, collapse = sepChar)
    # Make the plot
    p <- ggtree(tree, aes(color = group, linetype = linetype)) +
        scale_color_manual(values = groupColors) +
        guides(linetype = FALSE,
               color = guide_legend(override.aes = list(size = 3))) +
        theme(legend.position = "left") +
        ggtitle(label = attr(x, "site"), subtitle = subtitle)
    if (showTips) {
        p <- p + geom_tiplab()
    }
    return(p)
}

#' @rdname plotFunctions
#' @examples
#' x <- fixationPath(fixations)
#' plot(x)
#' @export
plot.fixationPath <- function(x,
                              y = TRUE,
                              ...) {
    tr <- attr(x, "SNPtracing")
    p <- ggtree(tr) +
        geom_tiplab(hjust = 0.5,
                    align = TRUE,
                    offset = 0.5) +
        theme_tree2()
    if (y) {
        p <- p + geom_label_repel(
            aes(x = branch, label = SNPs),
            fill = 'lightgreen',
            min.segment.length = 0,
            na.rm = TRUE
        )
    }
    return(p)
}
