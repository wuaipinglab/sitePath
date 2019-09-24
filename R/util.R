#' @rdname addMSA
#' @name addMSA
#' @title Prepare data for sitePath analysis
#' @description
#' sitePath requires both tree and sequence alignment to do the analysis.
#' \code{addMSA} wraps \code{read.alignment} function in \code{seqinr} package
#' and helps match names in tree and sequence alignment. Either
#' provide the file path to an alignment file and its format or an alignment
#' object from the return of \code{read.alignment} function. If both the file
#' path and alignment object are given, the function will use the sequence
#' in the alignment file.
#' @param tree
#' a \code{phylo} object. This commonly can be from tree parsing function
#' in \code{ape} or \code{ggtree}. All the \code{tip.label} should be
#' found in the sequence alignment.
#' @param msaPath
#' The file path to the multiple sequence alignment file
#' @param msaFormat
#' The format of the multiple sequence alignment file
#' @param alignment
#' an \code{alignment} object. This commonly can be
#' from sequence parsing function in the \code{seqinr} package.
#' Sequence names in the alignment should include all \code{tip.label}
#' in the tree
#' @return
#' \code{addMSA} returns a \code{phylo} object with matched
#' multiple sequence alignment
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
#' @importFrom seqinr read.alignment
#' @importFrom methods is
#' @export
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
    # Similarity matrix
    sim <- getSimilarityMatrix(align)
    dimnames(sim) <- list(tree$tip.label, tree$tip.label)
    attr(tree, "simMatrix") <- sim
    return(tree)
}

#' @rdname setSiteNumbering
#' @name setSiteNumbering
#' @title Set site numbering to the reference sequence
#' @description
#' A reference sequence can be used to define a global site numbering
#' scheme for multiple sequence alignment. The gap in the reference
#' will be skipped so the site ignored in numbering.
#' @param x
#' The object to set site numbering. It could be a \code{phylo} object
#' after \code{\link{addMSA}} or a \code{lineagePath} object. The function
#' for \code{fixaitonSites} and \code{multiFixationSites} will be
#' added in later version.
#' @param reference
#' Name of reference for site numbering. The name has to be one of the
#' sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar
#' The character to indicate gap. The numbering will skip the gapChar
#' for the reference sequence.
#' @param ... further arguments passed to or from other methods.
#' @return
#' A \code{phylo} object with site numbering mapped to reference sequence
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' tree <- addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
#' setSiteNumbering(tree)
#' @export
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

#' @name zikv_align
#' @title Multiple sequence alignment of Zika virus polyprotein
#' @description
#' The raw protein sequences were downloaded from ViPR database
#' (\url{https://www.viprbrc.org/}) and aliged using MAFFT.
#' with default settings.
#' @format a \code{alignment} object
#' @usage data(zikv_align)
#' @docType data
"zikv_align"

#' @name zikv_tree
#' @title Phylogenetic tree of Zika virus polyprotein
#' @description
#' Tree was built from \code{\link{zikv_align}} using RAxML
#' with default settings.The tip ANK57896 was used as
#' outgroup to root the tree.
#' @format a \code{phylo} object
#' @usage data(zikv_tree)
#' @docType data
"zikv_tree"

#' @name h3n2_align
#' @title Multiple sequence alignment of H3N2's HA protein
#' @description
#' The raw protein sequences were downloaded from NCBI database.
#' @format a \code{alignment} object
#' @usage data(h3n2_align)
#' @docType data
"h3n2_align"

#' @name h3n2_tree
#' @title Phylogenetic tree of H3N2's HA protein
#' @description
#' Tree was built from \code{\link{h3n2_align}} using RAxML
#' with default settings.
#' @format a \code{phylo} object
#' @usage data(h3n2_tree)
#' @docType data
"h3n2_tree"

#' @name zikv_align_reduced
#' @title Truncated data for runnable example
#' @description
#' This is a truncated version of \code{\link{zikv_align}}
#' @format a \code{alignment} object
#' @usage data(zikv_align_reduced)
#' @docType data
"zikv_align_reduced"

#' @name zikv_tree_reduced
#' @title Truncated data for runnable example
#' @description
#' This is a truncated version of \code{\link{zikv_tree}}
#' @format a \code{phylo} object
#' @usage data(zikv_tree_reduced)
#' @docType data
"zikv_tree_reduced"

#' @name h3n2_align_reduced
#' @title Truncated data for runnable example
#' @description
#' This is a truncated version of \code{\link{h3n2_align}}
#' @format a \code{alignment} object
#' @usage data(h3n2_align_reduced)
#' @docType data
"h3n2_align_reduced"

#' @name h3n2_tree_reduced
#' @title Truncated data for runnable example
#' @description
#' This is a truncated version of \code{\link{h3n2_tree}}
#' @format a \code{phylo} object
#' @usage data(h3n2_tree_reduced)
#' @docType data
"h3n2_tree_reduced"

#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

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

.childrenTips <- function(tree, node) {
    maxTip <- length(tree$tip.label)
    children <- integer(0)
    getChildren <- function(edges, parent) {
        children <<- c(children, parent[which(parent <= maxTip)])
        i <- which(edges[, 1] %in% parent)
        if (length(i) == 0L) {
            return(children)
        } else {
            parent <- edges[i, 2]
            return(getChildren(edges, parent))
        }
    }
    return(getChildren(tree$edge, node))
}

.extendPaths <- function(paths, tree) {
    # Store the attributes of 'paths'
    attrs <- attributes(paths)
    # Extend each path
    paths <- lapply(paths, function(p) {
        # Starting node
        sn <- p[length(p)]
        # nodePath from starting node to each children tip
        extended <- lapply(.childrenTips(tree, sn), function(t) {
            # The staring node is excluded from the nodePath
            nodepath(tree, sn, t)[-1]
        })
        el <- lengths(extended)
        # The length of longest nodePath
        ml <- max(el)
        # The longest nodePath
        longest <- extended[which(el == ml)]
        # Zip up the longest nodePath until hitting diverging node
        extended <- lapply(
            X = seq_len(ml),
            FUN = function(i) {
                unique(vapply(
                    longest,
                    FUN = "[[",
                    FUN.VALUE = integer(1),
                    i
                ))
            }
        )
        # The length of the element in 'extended' would be 1 if the
        # node is shared by the longest nodePath
        c(p, unlist(extended[which(lengths(extended) == 1)]))
    })
    # Give back the attributes of 'paths'
    attributes(paths) <- attrs
    return(paths)
}
