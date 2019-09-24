#' @rdname extractTips
#' @name extractTips
#' @title Extract info for fixation of a single site
#' @description
#' The result of \code{\link{fixationSites}} contains all the possible sites
#' with fixation mutation. The function \code{extractTips} retrieves
#' the name of the tips involved in the fixation.
#' @param x
#' A \code{fixationSites} or a \code{multiFixationSites} or
#' a \code{sitePath} object.
#' @param site
#' A site predicted to experience fixation.
#' @param select
#' For a site, there theoretically might be more than one
#' fixation on different lineages. You may use this argument
#' to extract for a specific fixation of a site. The default
#' is the first fixation of the site.
#' @param ...
#' Other arguments
#' @return
#' The function \code{extractTips} returns the name of the
#' tips involved in the fixation.
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' mutations <- fixationSites(lineagePath(tree))
#' extractTips(mutations, 139)
#' @export
extractTips.fixationSites <- function(x,
                                      site,
                                      select = 1,
                                      ...) {
    sp <- extractSite(x, site)
    return(.actualExtractTips(sp, select))
}

#' @rdname extractTips
#' @export
extractTips.multiFixationSites <- function(x,
                                           site,
                                           select = 1,
                                           ...) {
    sp <- extractSite(x, site)
    return(.actualExtractTips(sp, select))
}

#' @rdname extractTips
#' @export
extractTips.sitePath <- function(x, select = 1, ...) {
    return(.actualExtractTips(x, select))
}

#' @export
extractTips <- function(x, ...)
    UseMethod("extractTips")

#' @rdname extractTips
#' @name extractSite
#' @description
#' The function \code{extractSite} can be used to extract the fixation info
#' of a single site.
#' @return
#' The function \code{extractSite} returns a \code{sitePath} object
#' @examples
#' extractSite(mutations, 139)
#' @export
extractSite.fixationSites <- function(x, site, ...) {
    return(.actualExtractSite(x, site))
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

.actualExtractSite <- function(x, site) {
    site <- .checkSite(site)
    paths <- attr(x, "paths")
    tree <- attr(paths, "tree")
    tryCatch(
        expr = sp <- x[[as.character(site)]],
        error = function(e) {
            stop("\"site\": ", site, " is not found in \"x\".")
        }
    )
    attr(sp, "tree") <- tree
    return(sp)
}

.actualExtractTips <- function(sp, select) {
    tree <- attr(sp, "tree")
    if (select <= 0 || as.integer(select) != select) {
        stop("Please enter a single positive integer for \"select\"")
    }
    tryCatch(
        expr = sp <- sp[[select]],
        error = function(e) {
            cat(e, "\n")
            if (length(select))
                stop(
                    "The site: ",
                    attr(sp, "site"),
                    " has ",
                    length(sp),
                    " fixation(s). Please choose a number from 1 to ",
                    length(sp),
                    " for \"select\"."
                )
        }
    )
    res <- list()
    for (i in sp) {
        aa <- attr(i, "AA")
        attributes(aa) <- NULL
        i <- tree$tip.label[i]
        attr(i, "AA") <- aa
        res <- c(res, list(i))
    }
    return(res)
}

as.data.frame.fixationSites <- function(x, ...) {
    cat("Under development\n")
}

as.data.frame.multiFixationSites <- function(x, ...) {
    cat("Under development\n")
}
