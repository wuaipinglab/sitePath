#' @rdname extractTips
#' @name extractSite
#' @description The function \code{extractSite} can be used to extract the
#'   fixation info of a single site.
#' @return The function \code{extractSite} returns a \code{sitePath} object
#' @export
#' @examples
#' extractSite(mutations, 139)
extractSite.fixationSites <- function(x, site, ...) {
    return(.actualExtractSite(x, site))
}

.actualExtractSite <- function(x, site) {
    site <- .checkSite(site)
    tryCatch(
        expr = sp <- x[[as.character(site)]],
        error = function(e) {
            stop("\"site\": ", site, " is not found in \"x\".")
        }
    )
    return(sp)
}

#' @rdname extractTips
#' @name extractSite
#' @export
extractSite.multiFixationSites <- function(x, site, ...) {
    return(.actualExtractSite(x, site))
}

#' @export
extractSite <- function(x, site, ...)
    UseMethod("extractSite")
