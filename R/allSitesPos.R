#' @rdname allSitesPos
#' @title Retrieve position of all the sites
#' @description The function is a way to get position of the resulting sites
#'   from \code{\link{SNPsites}}, \code{\link{fixationSites}} and
#'   \code{\link{parallelSites}}. The numbering is consistent with what's being
#'   set by \code{\link{setSiteNumbering}}
#' @param x The object containing the sites from analysis
#' @param ... Other arguments
#' @return An integer vector for sites position
#' @export
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' tree <- addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
#' snp <- SNPsites(tree)
#' allSitesPos(snp)
allSitesPos <- function(x, ...) {
    UseMethod("allSitesPos")
}

#' @rdname allSitesPos
#' @export
allSitesPos.SNPsites <- function(x, ...) {
    as.integer(x)
}

#' @rdname allSitesPos
#' @export
allSitesPos.fixationSites <- function(x, ...) {
    as.integer(names(x))
}

#' @rdname allSitesPos
#' @export
allSitesPos.parallelSites <- function(x, ...) {
    as.integer(names(x))
}
