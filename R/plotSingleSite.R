#' @rdname plotSingleSite
#' @name plotSingleSite
#' @title Color the tree by a single site
#' @description For \code{\link{lineagePath}}, the tree will be colored
#'   according to the amino acid of the site. The color scheme tries to assign
#'   distinguishable color for each amino acid.
#' @param x A \code{fixationSites} object from \code{\link{fixationSites}} or
#'   the return from \code{\link{lineagePath}} function.
#' @param site For \code{lineagePath}, it can be any site within sequence
#'   length. For \code{fixationSites}, it is restrained to a predicted fixation
#'   site. The numbering is consistent with \code{reference} defined in the
#'   object.
#' @param showPath If plot the lineage result from \code{\link{lineagePath}}.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @param ... Arguments in \code{plot.phylo} functions and other arguments.
#' @return A ggplot object to make the plot. It does not behavior like the
#'   generic \code{\link{plot}} function.
#' @seealso \code{\link{plot.sitePath}}
#' @importFrom graphics plot legend
#' @importFrom ape ladderize getMRCA plot.phylo
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' paths <- lineagePath(tree)
#' plotSingleSite(paths, 139)
plotSingleSite.lineagePath <- function(x,
                                       site,
                                       showPath = FALSE,
                                       showTips = FALSE,
                                       ...) {
    group <- extractTips.lineagePath(x, site)
    # Use different color scheme depending on the sequence type
    seqType <- attr(x, "seqType")
    if (seqType == "AA") {
        names(group) <- AA_FULL_NAMES[names(group)]
        groupColors <- AA_COLORS
    } else {
        names(group) <- toupper(names(group))
        groupColors <- NT_COLORS
    }
    tree <- attr(x, "tree")
    group <- groupOTU(as_tibble(tree), group)
    group <- group[["group"]]
    # Set lineage nodes and non-lineage nodes as separate group
    if (showPath) {
        pathLabel <- ".lineage"
        # Color the path node black
        levels(group) <- c(levels(group), pathLabel)
        group[unique(unlist(x))] <- pathLabel
        lineageColor <- "black"
        names(lineageColor) <- pathLabel
        groupColors <- c(groupColors, lineageColor)
    }
    p <- ggtree(tree, aes(color = group)) +
        scale_color_manual(values = groupColors) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme(legend.position = "left") +
        ggtitle(site)
    if (showTips) {
        p <- p + geom_tiplab()
    }
    return(p)
}

#' @rdname plotSingleSite
#' @description For \code{\link{parallelSites}}, the tree will be colored
#'   according to the amino acid of the site if the mutation is not fixed.
#' @importFrom ggtree geom_tippoint
#' @export
plotSingleSite.parallelSites <- function(x,
                                         site,
                                         showPath = FALSE,
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
                na.rm = TRUE
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
        p <- p + geom_tippoint(mapping = aes(subset = sporadicTip))
    }
    return(p)
}

#' @rdname plotSingleSite
#' @description Visualize the \code{sitePath} object which can be extracted by
#'   using \code{\link{extractSite}} on the return of
#'   \code{\link{fixationSites}}.
#' @param y For a \code{sitePath} object, it can have more than one fixation
#'   path. This is to select which path to plot. The default is \code{NULL}
#'   which will plot all the paths. It is the same as \code{select} in
#'   \code{plotSingleSite}
#' @importFrom graphics title
#' @export
#' @seealso \code{\link{plotSingleSite}}, \code{\link{extractSite}}
#' @examples
#' fixations <- fixationSites(paths)
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
    if (attr(x, "seqType") == "AA") {
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

#' @rdname plotSingleSite
#' @description For \code{\link{fixationSites}}, it will color the ancestral
#'   tips in red, descendant tips in blue and excluded tips in grey.
#' @param select Select which fixation path in to plot. The default is NULL
#'   which will plot all the fixations.
#' @export
#' @examples
#' plotSingleSite(fixations, 139)
plotSingleSite.fixationSites <- function(x,
                                         site,
                                         select = NULL,
                                         ...) {
    plot.sitePath(x = .actualExtractSite(x, site), y = select, ...)
}

#' @rdname plotSingleSite
#' @export
plotSingleSite.multiFixationSites <- function(x,
                                              site,
                                              select = NULL,
                                              ...) {
    plot.sitePath(x = .actualExtractSite(x, site), y = select, ...)
}

#' @export
plotSingleSite <- function(x, ...)
    UseMethod("plotSingleSite")
