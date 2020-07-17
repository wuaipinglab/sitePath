#' @rdname plotSingleSite
#' @name plotSingleSite
#' @title Color the tree by a single site
#' @description For \code{lineagePath}, the tree will be colored according to
#'   the amino acid of the site. The color scheme tries to assign
#'   distinguishable color for each amino acid.
#' @param x A \code{fixationSites} object from \code{\link{fixationSites}} or
#'   the return from \code{\link{lineagePath}} function.
#' @param site For \code{lineagePath}, it can be any site within sequence
#'   length. For \code{fixationSites}, it is restrained to a predicted fixation
#'   site. The numbering is consistent with \code{reference} defined in the
#'   object.
#' @param showPath If plot the lineage result from lineagePath.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @param ... Arguments in \code{plot.phylo} functions and other arguments.
#' @return The function only makes plot and returns no value (It behaviors like
#'   the generic \code{\link{plot}} function).
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
    site <- .checkSite(site)
    tree <- attr(x, "tree")
    align <- attr(x, "align")
    align <- strsplit(tolower(align), "")
    reference <- attr(x, "msaNumbering")
    tryCatch(
        expr = site <- match.arg(as.character(site), seq_along(reference)),
        error = function(e) {
            stop("\"site\": ", site, " is not within the length of reference")
        }
    )
    siteComp <- vapply(
        X = align,
        FUN = "[[",
        FUN.VALUE = character(1),
        i = reference[site]
    )
    group <- list()
    for (i in seq_along(siteComp)) {
        group[[siteComp[[i]]]] <- c(group[[siteComp[[i]]]], i)
    }
    names(group) <- AA_FULL_NAMES[names(group)]
    groupColors <- AA_COLORS
    tree <- groupOTU(tree, group)
    # Set lineage nodes and non-lineage nodes as separate group
    if (showPath) {
        pathLabel <- ".Lineage"
        # Color the path node black
        levels(attr(tree, "group")) <-
            c(levels(attr(tree, "group")), pathLabel)
        attr(tree, "group")[unique(unlist(paths))] <- pathLabel
        lineageColor <- "black"
        names(lineageColor) <- pathLabel
        groupColors <- c(groupColors, lineageColor)
    }
    p <- ggtree(tree, aes(color = group)) +
        scale_color_manual(values = groupColors) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme(legend.position = "left")
    if (showTips) {
        p <- p + geom_tiplab()
    }
    return(p)
}

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
    subtitle <- character()
    group <- list()
    for (sp in sitePaths) {
        aaName <- character()
        for (tips in sp) {
            aa <- AA_FULL_NAMES[tolower(attr(tips, "AA"))]
            aaName <- c(aaName, aa)
            group[[aa]] <- c(group[[aa]], tips)
        }
        subtitle <-
            c(subtitle, paste0(AA_SHORT_NAMES[aaName], collapse = " -> "))
    }
    groupColors <- AA_COLORS
    excludedTips <- setdiff(seq_along(tree[["tip.label"]]),
                            unlist(group, use.names = FALSE))
    if (length(excludedTips) != 0) {
        excludedLabel <- ".excluded"
        group[[excludedLabel]] <- excludedTips
        excludedColor <- "#d3d3d3"
        names(excludedColor) <- excludedLabel
        groupColors <- c(groupColors, excludedColor)
    }
    groupColors[["0"]] <- "#d3d3d3"
    tree <- groupOTU(tree, group)
    sepChar <- "\n"
    if (sum(nchar(subtitle) <= 18)) {
        sepChar <- ", "
    }
    subtitle <- paste(subtitle, collapse = sepChar)
    p <- ggtree(tree, aes(color = group)) +
        scale_color_manual(values = groupColors) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
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
