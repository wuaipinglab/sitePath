#' @rdname plotFixationSites
#' @title Plot the result of fixation sites
#' @description Visualize the results of \code{\link{paraFixSites}}
#' @param x return from \code{\link{paraFixSites}}
#' @param site the number of the site according to
#'   \code{\link{setSiteNumbering}}
#' @param ... further arguments passed to or from other methods.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' paraFix <- paraFixSites(zikv_tree_reduced, alignment = zikv_align_reduced)
#' plotFixationSites(paraFix)
plotFixationSites <- function(x, ...) {
    UseMethod("plotFixationSites")
}

#' @rdname plotFixationSites
#' @param site the number of the site according to
#'   \code{\link{setSiteNumbering}}. If not provided, all sites will be plotted
#'   as labels on the tree
#' @export
plotFixationSites.fixationSites <- function(x, site = NULL, ...) {
    if (length(x)) {
        if (is.null(site)) {
            p <- plot.fixationSites(x)
        } else {
            p <- plotSingleSite.fixationSites(x, site = site)
        }
    } else {
        message("There is no fixation sites detected")
        paths <- attr(x, "paths")
        p <- plot.lineagePath(paths)
    }
    return(p)
}

#' @rdname plotFixationSites
#' @export
plotFixationSites.paraFixSites <- function(x, site = NULL, ...) {
    fixations <- attr(x, "allFixSites")
    p <- plotFixationSites.fixationSites(x = fixations,
                                         site = site,
                                         ...)
    return(p)
}
