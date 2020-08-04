#' @rdname extractSite
#' @name extractSite
#' @title Extract tips for a single site
#' @description The functions in \code{sitePath} usually include the results on
#'   more than one site. The function \code{extractSite} can be used to extract
#'   the predicted result on a single site.
#' @param x A \code{fixationSites} or a \code{parallelSites} object. More type
#'   will be supported in the later version.
#' @param site A site included in the result.
#' @param ... Other arguments
#' @return The predicted result of a single site
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' mutations <- fixationSites(lineagePath(tree))
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

#' @rdname extractSite
#' @name extractSite
#' @export
extractSite.multiFixationSites <- function(x, site, ...) {
    return(.actualExtractSite(x, site))
}

#' @export
extractSite <- function(x, site, ...)
    UseMethod("extractSite")
