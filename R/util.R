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
#' a \code{phylo} object. This commonly can be from tree paring function
#' in \code{ape} or \code{ggtree}. All the \code{tip.label} should be
#' found in the sequence alignment
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
#' msaPath <- system.file("extdata", "ZIKV.fasta", package = "sitePath")
#' addMSA(zikv_tree, msaPath = msaPath, msaFormat = "fasta")
#' @importFrom seqinr read.alignment
#' @importFrom methods is
#' @export
addMSA <- function(tree,
                   msaPath = "",
                   msaFormat = "",
                   alignment = NULL) {
    if (is(tree, "treedata")) {
        tree <- tree@phylo
    }
    if (!inherits(tree, "phylo")) {
        stop("\"tree\" is not class phylo")
    }
    if (file.exists(msaPath)) {
        alignment <- read.alignment(msaPath, msaFormat)
    } else if (is.null(alignment)) {
        stop(paste("There is no file in", msaPath))
    }
    attr(tree, "alignment") <- checkMatched(tree, alignment)
    return(tree)
}

#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

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
#' Tree was built from \code{zikv_align} using RAxML with default settings.
#' The tip “ANK57896” was used as outgroup to root the tree.
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
#' Tree was built from \code{h3n2_align} using RAxML with default settings.
#' @format a \code{phylo} object
#' @usage data(h3n2_tree)
#' @docType data
"h3n2_tree"

checkMatched <- function(tree, align) {
    if (!is(align, "alignment")) {
        stop("\"alignment\" is not class alignment")
    }
    m <- match(tree$tip.label, align$nam)
    if (any(is.na(m))) {
        stop("Tree tips and alignment names are not matched")
    } else {
        align <- align$seq[m]
        if (length(unique(nchar(align))) > 1) {
            stop("Sequence lengths are not the same in alignment")
        }
    }
    return(toupper(align))
}

checkReference <- function(tree, align, reference, gapChar) {
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

ChildrenTips <- function(tree, node) {
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

extendPaths <- function(paths, tree) {
    attrs <- attributes(paths)
    paths <- lapply(paths, function(p) {
        sn <- p[length(p)]
        extended <-
            lapply(ChildrenTips(tree, sn), function(t) {
                nodepath(tree, sn, t)[-1]
            })
        el <- lengths(extended)
        ml <- max(el)
        longest <- extended[which(el == ml)]
        extended <-
            lapply(seq_len(ml), function(i) {
                unique(vapply(longest, FUN = "[[", FUN.VALUE = 0, i))
            })
        c(p, unlist(extended[which(lengths(extended) == 1)]))
    })
    attributes(paths) <- attrs
    return(paths)
}
