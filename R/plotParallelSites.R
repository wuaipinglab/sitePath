#' @rdname plotParallelSites
#' @title Plot the result of fixation sites
#' @description Visualize the results of \code{\link{paraFixSites}}
#' @param x return from \code{\link{paraFixSites}}
#' @param site the number of the site according to
#'   \code{\link{setSiteNumbering}}
#' @param ... further arguments passed to or from other methods.
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' paraFix <- paraFixSites(zikv_tree, alignment = zikv_align)
#' plotParallelSites(paraFix)
plotParallelSites <- function(x, ...) {
    UseMethod("plotParallelSites")
}

#' @rdname plotParallelSites
#' @export
plotParallelSites.parallelSites <- function(x, site = NULL, ...) {
    if (length(x)) {
        if (is.null(site)) {
            p <- plot.parallelSites(x)
        } else {
            p <- plotSingleSite.parallelSites(x, site = site)
        }
    } else {
        message("There is no parallel sites detected")
        paths <- attr(x, "paths")
        p <- plot.lineagePath(paths)
    }
    return(p)
}

#' @rdname plotParallelSites
#' @export
plotParallelSites.paraFixSites <- function(x, site = NULL, ...) {
    paraSites <- attr(x, "allParaSites")
    p <- plotParallelSites.parallelSites(x = paraSites,
                                         site = site,
                                         ...)
    return(p)
}
