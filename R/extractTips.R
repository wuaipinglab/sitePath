#' @rdname extractTips
#' @name extractTips
#' @title Extract grouped tips for a single site
#' @description The result of \code{\link{fixationSites}} and \code{sitePath}
#'   contains all the possible sites with fixation mutation. The function
#'   \code{extractTips} retrieves the name of the tips involved in the fixation.
#' @param x A \code{fixationSites} or a \code{sitePath} object.
#' @param site A site predicted to experience fixation.
#' @param select For a site, there theoretically might be more than one fixation
#'   on different lineages. You may use this argument to extract for a specific
#'   fixation of a site. The default is the first fixation of the site.
#' @param ... Other arguments
#' @return Tree tips grouped as \code{\link{list}}
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' mutations <- fixationSites(lineagePath(tree))
#' extractTips(mutations, 139)
extractTips.sitePath <- function(x, select = 1, ...) {
    tree <- attr(x, "tree")
    if (select <= 0 || as.integer(select) != select) {
        stop("Please enter a single positive integer for \"select\"")
    }
    tryCatch(
        expr = x <- x[[select]],
        error = function(e) {
            if (length(select))
                stop(
                    "The site: ",
                    attr(x, "site"),
                    " has ",
                    length(x),
                    " fixation(s). Please choose a number from 1 to ",
                    length(x),
                    " for \"select\"."
                )
        }
    )
    res <- list()
    for (i in x) {
        aa <- attr(i, "AA")
        attributes(aa) <- NULL
        i <- tree[["tip.label"]][i]
        attr(i, "AA") <- aa
        res <- c(res, list(i))
    }
    return(res)
}

#' @rdname extractTips
#' @export
extractTips.fixationSites <- function(x,
                                      site,
                                      select = 1,
                                      ...) {
    sp <- extractSite.fixationSites(x, site)
    return(extractTips.sitePath(sp, select))
}

#' @rdname extractTips
#' @description For \code{\link{lineagePath}}, the function \code{extractTips}
#'   groups all the tree tips according to the amino acid/nucleotide of the
#'   \code{site}.
#' @export
extractTips.lineagePath <- function(x, site, ...) {
    site <- .checkSite(site)
    align <- attr(x, "align")
    align <- strsplit(tolower(align), "")
    reference <- attr(x, "msaNumbering")
    # Get the site index of the alignment
    tryCatch(
        expr = site <- reference[[site]],
        error = function(e) {
            stop("The site: ",
                 attr(x, "site"),
                 "is not within the length of sequence alignment.")
        }
    )
    # Group the tree tips by amino acid/nucleotide of the site
    group <- list()
    for (tipName in names(align)) {
        siteChar <- align[[tipName]][site]
        group[[siteChar]] <- c(group[[siteChar]], tipName)
    }
    return(group)
}

#' @rdname extractTips
#' @description For \code{\link{parallelSites}} and \code{sitePara} object, the
#'   function \code{extractTips} retrieve all the tips with parallel mutation.
#' @export
extractTips.sitePara <- function(x, ...) {
    # The duplicated nodes should be removed because they are not truly parallel
    toRemoveNodes <- unique(unlist(lapply(
        X = x,
        FUN = function(mut) {
            nodes <- names(mut)
            nodes[duplicated(nodes)]
        }
    )))
    # The parallel nodes in one pair of paths might be duplicated for the other
    # pair. And they should be remove too
    res <- list()
    for (mut in x) {
        for (node in names(mut)) {
            if (!node %in% toRemoveNodes) {
                res[[node]] <- mut[[node]]
            }
        }
    }
    return(res)
}

#' @rdname extractTips
#' @export
extractTips.parallelSites <- function(x, site, ...) {
    parallelMut <- extractSite.parallelSites(x, site)
    return(extractTips.sitePara(parallelMut))
}

#' @export
extractTips <- function(x, ...)
    UseMethod("extractTips")
