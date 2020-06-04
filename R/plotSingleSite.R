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
