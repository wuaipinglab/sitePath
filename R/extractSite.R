#' @rdname extractSite
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
extractSite <- function(x, site, ...) {
    UseMethod("extractSite")
}

#' @rdname extractSite
#' @export
extractSite.fixationSites <- function(x, site, ...) {
    return(.actualExtractSite(x, site))
}

.actualExtractSite <- function(x, site) {
    site <- .checkSite(site)
    sp <- x[[as.character(site)]]
    if (is.null(sp)) {
        stop("\"site\": ", site, " is not found in \"x\".")
    }
    return(sp)
}

.checkSite <- function(site) {
    if (!is.numeric(site) ||
        any(site <= 0) || as.integer(site) != site) {
        stop("Please enter positive integer value for \"site\"")
    }
    if (length(site) != 1) {
        site <- site[1]
        warning(
            "\"site\" has more than one element, ",
            "only the first element (",
            site,
            ") will be used."
        )
    }
    return(site)
}

#' @export
extractSite.parallelSites <- function(x, site, ...) {
    return(.actualExtractSite(x, site))
}
