#' @rdname extractTips
#' @name extractTips
#' @title Extract sitePath for a single site
#' @description
#' Retrieve the name of the tips involved in the fixation
#' @param x
#' A \code{fixationSites} or a \code{multiFixationSites} object.
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
#' The name of the tips involved in the fixation
#' @export
extractTips.fixationSites <- function(x, site, select = 1, ...) {
    site <- checkSite(site)
    tryCatch(
        expr = res <- x[[as.character(site)]],
        error = function(e) {
            stop(paste("\"site\":", site, "is not found in \"x\"."))
        }
    )
    if (select <= 0 || as.integer(select) != select) {
        stop("Please enter a single positive integer for \"select\"")
    }
    tryCatch(
        expr = res <- res[[select]],
        error = function(e) {
            cat(e, "\n")
            if (length(select))
                stop(
                    paste(
                        "The site:",
                        site,
                        "has",
                        length(x),
                        "fixation(s). Please choose a number from 1 to",
                        length(x),
                        "for \"select\"."
                    )
                )
        }
    )
    return(res)
}

#' @rdname extractTips
#' @export
extractTips.multiFixationSites <-
    function(x, site, select = 1, ...) {
        site <- checkSite(site)
        tryCatch(
            expr = res <- x[[as.character(site)]],
            error = function(e) {
                stop(paste("\"site\":", site, "is not found in \"x\"."))
            }
        )
        if (select <= 0 || as.integer(select) != select) {
            stop("Please enter a single positive integer for \"select\"")
        }
        tryCatch(
            expr = res <- res[[select]],
            error = function(e) {
                cat(e, "\n")
                if (length(select))
                    stop(
                        paste(
                            "The site:",
                            site,
                            "has",
                            length(x),
                            "fixation(s). Please choose a integer from 1 to",
                            length(x),
                            "for \"select\"."
                        )
                    )
            }
        )
        return(res)
    }

#' @export
extractTips <- function(x, site, select, ...)
    UseMethod("extractTips")

#' @export
as.data.frame.fixationSites <- function(x, ...) {
    cat("Under development\n")
    paths <- attr(x, "paths")
    tree <- attr(paths, "tree")
    reference <- attr(paths, "reference")
    res <- as.data.frame(matrix(
        nrow = length(tree$tip.label),
        ncol = length(reference),
        dimnames = list(tree$tip.label, seq_along(reference))
    ))
    for (site in names(x)) {
        for (tips in x[[site]]) {
            fixedAA <- AA_FULL_NAMES[tolower(names(tips))]
            for (tip in tips[[2]]) {
                res[tip, site] <- fixedAA
            }
        }
    }
    whichNA <- is.na(res)
    res <-
        res[rowSums(whichNA) < ncol(res), colSums(whichNA) < nrow(res)]
    return(res)
}

#' @export
as.data.frame.multiFixationSites <- function(x, ...) {
    cat("Under development\n")
}
