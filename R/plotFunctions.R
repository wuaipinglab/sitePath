#' @importFrom graphics plot
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 aes guides theme ggtitle guide_legend
#' @importFrom ggplot2 scale_size scale_color_manual
#' @importFrom ggplot2 GeomText
#' @importFrom ape Nnode
#' @importFrom tidytree groupOTU
#' @importFrom ggtree ggtree geom_tiplab theme_tree2
#' @importFrom ggrepel geom_label_repel

#' @rdname plotFunctions
#' @title Visualize the results
#' @description The plot function to visualize the return of functions in the
#'   package. The underlying function applies \code{\link{ggplot2}}. The
#'   function name \code{plot} is used to keep the compatibility with previous
#'   versions, but they do not behave like the generic \code{\link{plot}}
#'   function since 1.5.4.
#' @description A \code{\link{phyMSAmatched}} object will be plotted as a tree
#'   diagram.
#' @param x The object to plot.
#' @param y Whether to show the fixation mutation between clusters. For
#'   \code{lineagePath} object and \code{sitePath} object, it is deprecated and
#'   no longer have effect since 1.5.4.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @param ... Other arguments. Since 1.5.4, the function uses
#'   \code{\link{ggtree}} as the base function to make plots so the arguments in
#'   \code{plot.phylo} will no longer work.
#' @return A ggplot object to make the plot.
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' plot(tree)
plot.phyMSAmatched <- function(x, y = TRUE, ...) {
    p <- ggtree(as.phylo.phyMSAmatched(x))
    return(p)
}

#' @rdname plotFunctions
#' @description A \code{\link{lineagePath}} object will be plotted as a tree
#'   diagram and paths are black solid line while the trimmed nodes and tips
#'   will use gray dashed line.
#' @examples
#' paths <- lineagePath(tree)
#' plot(paths)
#' @export
plot.lineagePath <- function(x,
                             y = TRUE,
                             showTips = FALSE,
                             ...) {
    tree <- attr(x, "tree")
    # Get number of ancestral nodes plus tip nodes
    nNodes <- Nnode(tree, internal.only = FALSE)
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
#' @description A \code{\link{fixationSites}} object will be plotted as original
#'   phylogenetic tree marked with fixation substitutions.
#' @param tipsGrouping A \code{list} to hold the grouping of tips for how the
#'   tree will be colored.
#' @examples
#' fixations <- fixationSites(paths)
#' plot(fixations)
#' @export
plot.fixationSites <- function(x,
                               y = TRUE,
                               tipsGrouping = NULL,
                               ...) {
    tree <- as.treedata.fixationSites(x, .node = tipsGrouping)
    grp <- levels(tree@data[["group"]])
    grp <- grp[which(grp != "0")]
    groupColors <-
        colorRampPalette(brewer.pal(9, "Set1"))(length(grp))
    names(groupColors) <- grp
    groupColors["0"] <- "black"

    p <- ggtree(tree, aes(color = group)) +
        scale_color_manual(values = groupColors) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme(legend.position = "left")
    if (y) {
        p <- p + geom_label_repel(
            aes(x = branch, label = SNPs),
            fill = "lightgreen",
            color = "black",
            min.segment.length = 0,
            na.rm = TRUE,
            ...
        )
    }
    return(p)
}

#' @rdname plotFunctions
#' @description A \code{sitePath} object can be extracted by using
#'   \code{\link{extractSite}} on the return of \code{\link{fixationSites}}.
#' @param select For a \code{sitePath} object, it can have result on more than
#'   one evolution pathway. This is to select which path to plot. The default is
#'   \code{NULL} which will plot all the paths. It is the same as \code{select}
#'   in \code{\link{plotSingleSite}}.
#' @export
#' @examples
#' sp <- extractSite(fixations, 139)
#' plot(sp)
plot.sitePath <- function(x,
                          y = NULL,
                          select = NULL,
                          showTips = FALSE,
                          ...) {
    if (is.null(select)) {
        sitePaths <- x[]
    } else {
        if (any(length(x) < select)) {
            stop("Some in 'select' is out of bounds for 'x'.")
        }
        sitePaths <- x[select]
    }
    tree <- attr(x, "tree")
    siteName <- attr(x, "site")
    # Specify the color of mutations by pre-defined color set.
    groupColors <- .siteColorScheme(attr(x, "seqType"))
    # Collect the fixation mutation for each evolutionary pathway
    subtitle <- character()
    fixationMut <- character()
    group <- list()
    for (sp in sitePaths) {
        # The first group is separately done
        tips <- sp[[1]]
        prevAA <- toupper(attr(tips, "AA"))
        aaName <- prevAA
        group[[prevAA]] <- c(group[[prevAA]], tips)
        # Do the rest
        for (i in seq_along(sp)[-1]) {
            tips <- sp[[i]]
            aa <- toupper(attr(tips, "AA"))
            aaName <- c(aaName, aa)
            group[[aa]] <- c(group[[aa]], tips)
            fixationMut[names(sp[i])] <- paste0(prevAA,
                                                siteName,
                                                aa)
            prevAA <- aa
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
    # Annotate the mutation on the tree
    if (is.null(y) || !y) {
        p <- ggtree(tree, aes(color = group, linetype = linetype))
    } else {
        tree <- .annotateSNPonTree(tree, fixationMut)
        p <- ggtree(tree, aes(color = group, linetype = linetype)) +
            geom_label_repel(
                aes(x = branch, label = SNPs),
                fill = 'lightgreen',
                color = "black",
                min.segment.length = 0,
                na.rm = TRUE,
                size = GeomText[["default_aes"]][["size"]]
            )
    }
    # Make the plot
    p <- p + scale_color_manual(values = groupColors) +
        guides(linetype = FALSE,
               color = guide_legend(override.aes = list(size = 3))) +
        theme(legend.position = "left") +
        ggtitle(label = siteName, subtitle = subtitle)
    if (showTips) {
        p <- p + geom_tiplab()
    }
    return(p)
}

#' @rdname plotFunctions
#' @description A \code{\link{fixationIndels}} object will be plotted as original
#'   phylogenetic tree marked with indel fixation.
plot.fixationIndels <- function(x, y = TRUE, ...) {
}

#' @rdname plotFunctions
#' @description A \code{\link{fixationPath}} object will be plotted as a
#'   \code{phylo} object. The tips are clustered according to the fixation
#'   sites. The transition of fixation sites will be plotted as a phylogenetic
#'   tree. The length of each branch represents the number of fixation mutation
#'   between two clusters.
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
