#' @rdname addMSA
#' @name addMSA
#' @title Prepare data for sitePath analysis
#' @description sitePath requires both tree and sequence alignment to do the
#'   analysis. \code{addMSA} wraps \code{read.alignment} function in
#'   \code{seqinr} package and helps match names in tree and sequence alignment.
#'   Either provide the file path to an alignment file and its format or an
#'   alignment object from the return of \code{read.alignment} function. If both
#'   the file path and alignment object are given, the function will use the
#'   sequence in the alignment file.
#' @param tree a \code{phylo} object. This commonly can be from tree parsing
#'   function in \code{ape} or \code{ggtree}. All the \code{tip.label} should be
#'   found in the sequence alignment.
#' @param msaPath The file path to the multiple sequence alignment file
#' @param msaFormat The format of the multiple sequence alignment file
#' @param alignment an \code{alignment} object. This commonly can be from
#'   sequence parsing function in the \code{seqinr} package. Sequence names in
#'   the alignment should include all \code{tip.label} in the tree
#' @return \code{addMSA} returns a \code{phylo} object with matched multiple
#'   sequence alignment
#' @importFrom seqinr read.alignment
#' @importFrom methods is
#' @export
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
addMSA <- function(tree,
                   msaPath = "",
                   msaFormat = "",
                   alignment = NULL) {
    # Get/test the tree object
    if (is(tree, "treedata")) {
        tree <- tree@phylo
    }
    if (!inherits(tree, "phylo")) {
        stop("\"tree\" is not class phylo")
    }
    # Read alignment from the file if applicable
    if (file.exists(msaPath)) {
        alignment <- read.alignment(msaPath, msaFormat)
    } else if (is.null(alignment)) {
        stop("Alignment file \"", msaPath, "\" does not exist")
    }
    # Test the alignment object
    if (!is(alignment, "alignment")) {
        stop("\"alignment\" is not class alignment")
    }
    # Map the names between tree and alignment
    m <- match(tree$tip.label, alignment$nam)
    if (any(is.na(m))) {
        stop("Tree tips and alignment names are not matched")
    } else {
        align <- toupper(alignment$seq[m])
        if (length(unique(nchar(align))) > 1) {
            stop("Sequence lengths are not the same in alignment")
        }
    }
    attr(tree, "align") <- align
    # Use the numbering of MSA as the default site numbering
    attr(tree, "reference") <-
        .checkReference(tree, align, NULL, "-")
    return(tree)
}

.checkReference <- function(tree, align, reference, gapChar) {
    if (is.null(reference)) {
        reference <- seq_len(nchar(align[1]))
    } else {
        if (!is.character(gapChar) ||
            nchar(gapChar) != 1 || length(gapChar) != 1) {
            stop("\"gapChar\" only accepts one single character")
        }
        reference <-
            getReference(align[which(tree$tip.label == reference)], gapChar)
    }
    return(reference)
}

#' @rdname setSiteNumbering
#' @name setSiteNumbering
#' @title Set site numbering to the reference sequence
#' @description A reference sequence can be used to define a global site
#'   numbering scheme for multiple sequence alignment. The gap in the reference
#'   will be skipped so the site ignored in numbering.
#' @param x The object to set site numbering. It could be a \code{phylo} object
#'   after \code{\link{addMSA}} or a \code{lineagePath} object. The function for
#'   \code{fixaitonSites} and \code{multiFixationSites} will be added in later
#'   version.
#' @param reference Name of reference for site numbering. The name has to be one
#'   of the sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar The character to indicate gap. The numbering will skip the
#'   gapChar for the reference sequence.
#' @param ... further arguments passed to or from other methods.
#' @return A \code{phylo} object with site numbering mapped to reference
#'   sequence
#' @export
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' tree <- addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
#' setSiteNumbering(tree)
setSiteNumbering.phylo <- function(x,
                                   reference = NULL,
                                   gapChar = "-",
                                   ...) {
    align <- attr(x, "align")
    siteMapping <- .checkReference(x, align, reference, gapChar)
    attr(siteMapping, "refSeqName") <- reference
    attr(x, "reference") <- siteMapping
    return(x)
}

#' @rdname setSiteNumbering
#' @name setSiteNumbering
#' @export
setSiteNumbering.lineagePath <- function(x,
                                         reference = NULL,
                                         gapChar = "-",
                                         ...) {
    align <- attr(x, "align")
    tree <- attr(x, "tree")
    siteMapping <-
        .checkReference(tree, align, reference, gapChar)
    attr(siteMapping, "refSeqName") <- reference
    attr(x, "reference") <- siteMapping
    return(x)
}

#' @export
setSiteNumbering <- function(x, reference, gapChar, ...)
    UseMethod("setSiteNumbering")

# TODO: Need a class called "reference" for site numbering
