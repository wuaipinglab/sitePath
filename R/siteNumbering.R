#' @rdname setSiteNumbering
#' @name setSiteNumbering
#' @title Set site numbering to the reference sequence
#' @description A reference sequence can be used to define a global site
#'   numbering scheme for multiple sequence alignment. The gap in the reference
#'   will be skipped so the site ignored in numbering.
#' @param x The object to set site numbering. It could be a
#'   \code{\link{phyMSAmatched}} or a \code{lineagePath} object.
#' @param reference Name of reference for site numbering. The name has to be one
#'   of the sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar The character to indicate gap. The numbering will skip the
#'   gapChar for the reference sequence.
#' @param ... further arguments passed to or from other methods.
#' @return The input \code{x} with numbering mapped to \code{reference}.
#' @export
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' tree <- addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
#' setSiteNumbering(tree)
setSiteNumbering.phyMSAmatched <- function(x,
                                           reference = NULL,
                                           gapChar = "-",
                                           ...) {
    res <- .checkReference(x, reference, gapChar)
    return(res)
}

.checkReference <- function(x, reference, gapChar) {
    x <- .phyMSAmatch(x)
    align <- attr(x, "align")
    if (is.null(reference)) {
        # Use numbering of MSA when now reference provided
        attr(x, "msaNumbering") <- seq_len(nchar(align[1]))
    } else if (!is.character(gapChar) || length(gapChar) != 1  ||
               nchar(gapChar) != 1) {
        stop("\"gapChar\" only accepts one single character")
    } else {
        # The length of 'msaNumbering' is same as the reference sequence and the
        # numbers are the site indices of MSA
        attr(x, "msaNumbering") <-
            getReference(align[reference], gapChar)
        attr(x, "reference") <- reference
    }
    return(x)
}

.phyMSAmatch <- function(x) {
    # 'x' could be any object with 'tree' and 'align' attributes.
    # ('msaNumbering' and 'reference' attributes are optional)

    # Safeguard alignment sequence
    align <- attr(x, "align")
    if (is.null(align)) {
        stop("No alignment sequence found.")
    } else if (length(unique(nchar(align))) > 1) {
        stop("Sequence lengths are not the same in alignment.")
    }
    # Saftguard tree object
    tree <- attr(x, "tree")
    if (is.null(tree)) {
        stop("No phylogenetic tree found.")
    } else if (!inherits(tree, "phylo")) {
        stop("\"tree\" is not class phylo.")
    }
    # Map the names between tree and alignment
    m <- match(tree[["tip.label"]], names(align))
    if (any(is.na(m))) {
        stop("Tree tips and alignment names are not matched.")
    }
    # Update 'align' attribute as matched with tree tips
    attr(x, "align") <- align[m]
    return(x)
}

#' @rdname setSiteNumbering
#' @export
setSiteNumbering.lineagePath <- function(x,
                                         reference = NULL,
                                         gapChar = "-",
                                         ...) {
    res <- .checkReference(x, reference, gapChar)
    return(res)
}

setSiteNumbering.fixationSites <- function(x,
                                           reference = NULL,
                                           gapChar = '-',
                                           ...) {
    site2newRef <- cds2genome
    names(x) <- vapply(
        X = names(x),
        FUN = function(n) {
            site2newRef[[n]]
        },
        FUN.VALUE = integer(1)
    )
    # Update the 'site' attr in each 'sitePath'
    for (n in names(x)) {
        attr(x[[n]], "site") <- as.integer(n)
    }
    for (gpIndex in seq_along(attr(x, "clustersByPath"))) {
        for (i in seq_along(attr(x, "clustersByPath")[[gpIndex]])) {
            oldSiteName <-
                names(attr(attr(x, "clustersByPath")[[gpIndex]][[i]], "site"))
            names(attr(attr(x, "clustersByPath")[[gpIndex]][[i]], "site")) <-
                site2newRef[oldSiteName]
            toMerge <-
                attr(attr(x, "clustersByPath")[[gpIndex]][[i]], "toMerge")
            if (!is.null(toMerge)) {
                attr(attr(x, "clustersByPath")[[gpIndex]][[i]], "toMerge") <-
                    lapply(
                        X = toMerge,
                        FUN = function(sites) {
                            oldSiteName <- names(sites)
                            names(sites) <- site2newRef[oldSiteName]
                            return(sites)
                        }
                    )
            }
        }
    }
    return(x)
}

setSiteNumbering.fixationPath <- function(x,
                                          reference = NULL,
                                          gapChar = '-',
                                          ...) {
    site2newRef <- cds2genome
    attr(x, "SNPtracing")@data[["SNPs"]] <- vapply(
        X = strsplit(attr(x, "SNPtracing")@data[["SNPs"]], ", "),
        FUN = function(allSNPs) {
            if (is.na(allSNPs)) {
                return(allSNPs)
            } else {
                res <- vapply(
                    X = allSNPs,
                    FUN = function(snp) {
                        charLen <- nchar(snp)
                        preMut <- substr(snp, 1, 1)
                        postMut <- substr(snp, charLen, charLen)
                        site <-
                            site2newRef[[substr(snp, 2, charLen - 1)]]
                        return(paste0(preMut, site, postMut))
                    },
                    FUN.VALUE = character(1)
                )
                return(paste0(res, collapse = ", "))
            }
        },
        FUN.VALUE = character(1)
    )
    return(x)
}

#' @export
setSiteNumbering <- function(x, reference, gapChar, ...)
    UseMethod("setSiteNumbering")