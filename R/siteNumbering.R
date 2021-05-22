#' @rdname setSiteNumbering
#' @title Set site numbering to the reference sequence
#' @description A reference sequence can be used to define a global site
#'   numbering scheme for multiple sequence alignment. The gap in the reference
#'   sequence will be skipped for the numbering. Also, the site that is gap or
#'   amino acid/nucleotide for too many tips will be ignored but won't affect
#'   numbering.
#' @param x The object to set site numbering. It could be a
#'   \code{\link{phyMSAmatched}} or a \code{\link{lineagePath}} object.
#' @param reference Name of reference for site numbering. The name has to be one
#'   of the sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar The character to indicate gap. The numbering will skip the
#'   \code{gapChar} for the reference sequence.
#' @param minSkipSize The minimum number of tips to have gap or ambiguous amino
#'   acid/nucleotide for a site to be ignored in other analysis. This will not
#'   affect the numbering. The default is 0.8.
#' @param ... Further arguments passed to or from other methods.
#' @return The input \code{x} with numbering mapped to \code{reference}.
#' @export
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' tree <- addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
#' setSiteNumbering(tree)
setSiteNumbering <- function(x, reference, gapChar, ...) {
    UseMethod("setSiteNumbering")
}

#' @rdname setSiteNumbering
#' @export
setSiteNumbering.phyMSAmatched <- function(x,
                                           reference = NULL,
                                           gapChar = "-",
                                           minSkipSize = NULL,
                                           ...) {
    extraArgs <- list(...)
    usedSites <- extraArgs[["usedSites"]]
    x <- .phyMSAmatch(x)
    align <- attr(x, "align")
    if (is.null(reference)) {
        # Use numbering of MSA when now reference provided
        msaNumbering <- seq_len(nchar(align[1]))
    } else if (!is.character(gapChar) || length(gapChar) != 1 ||
               nchar(gapChar) != 1) {
        stop("\"gapChar\" only accepts one single character")
    } else {
        # The length of 'msaNumbering' is same as the reference sequence and the
        # numbers are the site indexes of MSA
        msaNumbering <- getReference(align[reference], gapChar)
        attr(x, "reference") <- reference
    }
    if (!is.null(usedSites)) {
        if (any(!is.numeric(usedSites)) || any(usedSites <= 0)) {
            stop("Please enter positive number for \"usedSites\"")
        } else if (any(usedSites > max(msaNumbering))) {
            stop("\"usedSites\" exceed length of the reference.")
        }
        msaNumbering <- msaNumbering[usedSites]
    }
    attr(x, "msaNumbering") <- msaNumbering
    attr(x, "gapChar") <- gapChar
    # Check and set the minimum size to remove a site
    if (is.null(minSkipSize)) {
        minSkipSize <- 0.8 * length(align)
    } else if (minSkipSize <= 0) {
        stop("'minSkipSize' cannot be zero or negative.")
    } else if (minSkipSize < 1 && minSkipSize > 0) {
        minSkipSize <- minSkipSize * length(align)
    }
    unambiguous <- .unambiguousChars(x)
    attr(x, "loci") <- which(vapply(
        X = attr(x, "msaNumbering") - 1,
        FUN = function(s) {
            # Summarize the number amino acid/nucleotide for the site
            siteSummary <- tableAA(align, s)
            # Drop the site if completely conserved
            if (length(siteSummary) == 1) {
                return(FALSE)
            } else {
                siteChars <- names(siteSummary)
                ambiguousChars <-
                    siteChars[which(!siteChars %in% unambiguous)]
                if (length(ambiguousChars) != 0) {
                    # Drop the site if number of tips having gap or ambiguous on
                    # this site is over the threshold
                    if (sum(siteSummary[ambiguousChars]) >= minSkipSize) {
                        return(FALSE)
                    }
                }
            }
            return(TRUE)
        },
        FUN.VALUE = logical(1)
    ))
    return(x)
}

.unambiguousChars <- function(x) {
    gapChar <- attr(x, "gapChar")
    if (attr(x, "seqType") == "AA") {
        res <- setdiff(AA_UNAMBIGUOUS, gapChar)
    } else {
        res <- setdiff(NT_UNAMBIGUOUS, gapChar)
    }
    return(res)
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
    # Safeguard tree object
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

setSiteNumbering.fixationSites <- function(x,
                                           reference = NULL,
                                           gapChar = '-',
                                           ...) {
    siteMapping <- list(...)[["siteMapping"]]
    names(x) <- vapply(
        X = names(x),
        FUN = function(n) {
            as.character(siteMapping[[n]])
        },
        FUN.VALUE = character(1)
    )
    # Update the 'site' attr in each 'sitePath'
    for (n in names(x)) {
        attr(x[[n]], "site") <- n
    }
    for (gpIndex in seq_along(attr(x, "clustersByPath"))) {
        for (i in seq_along(attr(x, "clustersByPath")[[gpIndex]])) {
            oldSiteName <-
                names(attr(attr(x, "clustersByPath")[[gpIndex]][[i]], "AA"))
            names(attr(attr(x, "clustersByPath")[[gpIndex]][[i]], "AA")) <-
                siteMapping[oldSiteName]
            toMerge <-
                attr(attr(x, "clustersByPath")[[gpIndex]][[i]], "toMerge")
            if (!is.null(toMerge)) {
                attr(attr(x, "clustersByPath")[[gpIndex]][[i]], "toMerge") <-
                    lapply(
                        X = toMerge,
                        FUN = function(sites) {
                            oldSiteName <- names(sites)
                            names(sites) <- siteMapping[oldSiteName]
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
    siteMapping <- list(...)[["siteMapping"]]
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
                            siteMapping[[substr(snp, 2, charLen - 1)]]
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

setSiteNumbering.paraFixSites <- function(x,
                                          reference = NULL,
                                          gapChar = '-',
                                          ...) {
    fixedSites <- attr(x, "fixSites")
    attr(x, "fixSites") <-
        setSiteNumbering.fixationSites(fixedSites,
                                       reference,
                                       gapChar,
                                       ...)
    return(x)
}
