#' @rdname sitePath-deprecated
#' @name sitePath-deprecated
#' @aliases multiFixationSites
#' @title Deprecated functions in package \sQuote{sitePath}
#' @description These functions are provided for compatibility with older
#'   versions of \sQuote{sitePath} only, and will be defunct at the next
#'   release.
#' @details The following functions are deprecated and will be made defunct; use
#'   the replacement indicated below: \itemize{ \item{multiFixationSites:
#'   \code{\link{fixationSites}}} }
#' @importFrom utils flush.console
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
NULL

#' @export
multiFixationSites.lineagePath <- function(paths,
                                           samplingSize = NULL,
                                           samplingTimes = 100,
                                           minEffectiveSize = 0,
                                           searchDepth = 1,
                                           method = c("compare", "insert", "delete"),
                                           ...) {
    .Deprecated("fixationSites")
    res <- fixationSites.lineagePath(
        paths = paths,
        minEffectiveSize = minEffectiveSize,
        searchDepth = searchDepth,
        method = method,
        ...
    )
    return(res)
}

#' @export
multiFixationSites <- function(paths, ...)
    UseMethod("multiFixationSites")
