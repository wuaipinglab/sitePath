#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

#' @importFrom seqinr read.alignment 
#' @importFrom ape read.tree
#' @importFrom ape root
NULL

#' @name zikv_align
#' @title Multiple sequence alignment of Zika virus polyprotein
#' @description 
#' The raw protein sequences were downloaded from ViPR database 
#' (\url{https://www.viprbrc.org/}) and aliged using MAFFT.
#' with default settings.
#' @format a \code{alignment} object
"zikv_align"

#' @name zikv_tree
#' @title Phylogenetic tree of Zika virus polyprotein
#' @description 
#' Tree was built from \code{zikv_align} using RAxML with default settings.
#' The tip “ANK57896” was used as outgroup to root the tree.
#' @format a \code{phylo} object
"zikv_tree"

#' @name zikv_paths
#' @title Phylogenetic lineages of Zika virus polyprotein
#' @description 
#' A list of node path to represent phylogenetic lineages of \code{zikv_tree}.
#' The minimal pairwise sequence similarity was set to 0.996.
#' @format a \code{sitePath} object
"zikv_paths"

#' @name zikv_fixations
#' @title Fixation sites found in Zika virus polyprotein
#' @description 
#' The fixation mutation predicted on \code{zikv_paths} by 
#' \code{\link{fixationSites}} with default setting.
#' @format a \code{fixationSites} object
"zikv_fixations"

#' @rdname treemer
#' @name treemer
#' @title Topology-dependent tree trimming
#' @description
#' \code{groupTips} uses sequence similarity to group tree tips.
#' Members in a group are always constrained to share the same 
#' ancestral node. Similarity between two tips is derived from 
#' their multiple sequence alignment. The site will not be counted
#' into total length if both are gap. Similarity is calculated
#' as number of matched divided by the corrected total length.
#' So far the detection of divergence is based on one simple rule:
#' the miminal pairwise similarity. The two branches are decided to
#' be divergent if the similarity is lower than the threshold. 
#' (Other more statistical approaches such as Kolmogorov-Smirnov 
#' Tests among pair-wise distance could be introduced in the future)
#' @param tree
#' a \code{phylo} object. This commonly can be from tree paring function
#' in \code{ape} or \code{ggtree}
#' All the \code{tip.label} should be found in the alignment
#' @param align
#' an \code{alignment} object. This commonly can be
#' from sequence parsing function in \code{ape} or \code{seqinr}
#' and many others. Sequence names in the alignment
#' should include all \code{tip.label} in the tree
#' @param similarity similarity threshold for tree trimming
#' @param simMatrix a diagonal matrix of similarity between sequences
#' @param forbidTrivial does not allow trivial trimming
#' @param tipnames if return as tipnames
#' @importFrom ape nodepath
#' @examples 
#' data("zikv_tree")
#' data("zikv_align")
#' groupTips(zikv_tree, zikv_align, 0.996)
#' @return grouping of tips
#' @export
groupTips <- function(tree,
                      align,
                      similarity,
                      simMatrix = NULL,
                      forbidTrivial = TRUE,
                      tipnames = TRUE) {
    simMatrix <- sortSimMatrix(tree, simMatrix)
    grouping <- trimTree(nodepath(tree),
                         checkMatched(tree, align),
                         simMatrix,
                         similarity,
                         TRUE)
    if (length(grouping) == 1 && forbidTrivial) {
        warning(paste0(
            similarity,
            " is too low of a cutoff resulting in trivial trimming"
        ))
    }
    if (tipnames) {
        return(lapply(grouping, function(g) {
            tree$tip.label[g]
        }))
    } else {
        return(grouping)
    }
}

#' @rdname treemer
#' @description
#' \code{sitePath} finds the lineages of a phylogenetic tree providing
#' the corresponding sequence alignment. This is done by trimming
#' the tree to the ancestor node of tips in each group and then find
#' the bifurcated terminals of the trimmed tree. The \code{\link{nodepath}}
#' between root node and the bifurcated terminals is the lineages.
#' In order to extend the search of mutational site. The lineages
#' will tag some of its trailing nodes. Here nodes up to
#' the ancestor of the tips with the longest \code{\link{nodepath}} 
#' are added.
#' @importFrom ape nodepath
#' @examples 
#' data("zikv_tree")
#' data("zikv_align")
#' sitePath(zikv_tree, zikv_align, 0.996)
#' @return path represent by node number
#' @export
sitePath <- function(tree,
                     align,
                     similarity,
                     simMatrix = NULL,
                     forbidTrivial = TRUE) {
    simMatrix <- sortSimMatrix(tree, simMatrix)
    align <- checkMatched(tree, align)
    # nodepath after trimming
    trimmedPaths <-
        unique(trimTree(nodepath(tree), align, simMatrix, similarity, FALSE))
    # get the bifurcated pre-terminal nodes and their path to the root
    # those paths are the so-called sitePaths (isolated)
    paths <- lapply(trimmedPaths, function(p)
        p[1:(length(p) - 1)])
    paths <-
        unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
    if (length(paths) == 0 && forbidTrivial) {
        warning(paste0(
            similarity,
            " is too low of a cutoff resulting in trivial trimming"
        ))
    }
    attr(paths, "tree") <- tree
    attr(paths, "align") <- align
    attr(paths, "class") <- "sitePath"
    return(paths)
}

# customGroupTips <- function(tree,
#                             align,
#                             pvalue = 0.05,
#                             forbidTrivial = TRUE,
#                             tipnames = TRUE) {
#     qualifyFunc <- function(x, y) {
#         return(ks.test(x, y)$p.value > pvalue)
#     }
#     grouping <- customTrimTree(
#         nodepath(tree),
#         checkMatched(tree, align),
#         cophenetic.phylo(tree),
#         tree$edge,
#         qualifyFunc,
#         TRUE
#     )
#     if (length(grouping) == 1 && forbidTrivial) {
#         warning(paste0(
#             pvalue,
#             " is too low of a cutoff resulting in trivial trimming"
#         ))
#     }
#     if (tipnames) {
#         return(lapply(grouping, function(g) {
#             tree$tip.label[g]
#         }))
#     } else {
#         return(grouping)
#     }
# }
# 
# customSitePath <- function(tree,
#                            align,
#                            pvalue = 0.05,
#                            forbidTrivial = TRUE) {
#     qualifyFunc <- function(x, y) {
#         return(ks.test(x, y)$p.value > pvalue)
#     }
#     # nodepath after trimming
#     trimmedPaths <- unique(
#         customTrimTree(
#             nodepath(tree),
#             align <- checkMatched(tree, align),
#             cophenetic.phylo(tree),
#             tree$edge,
#             qualifyFunc,
#             FALSE
#         )
#     )
#     # get the bifurcated pre-terminal nodes and their path to the root
#     # those paths are the so-called sitePaths (isolated)
#     paths <- lapply(trimmedPaths, function(p)
#         p[1:(length(p) - 1)])
#     paths <-
#         unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
#     if (length(paths) == 0 && forbidTrivial) {
#         warning(paste0(
#             pvalue,
#             " is too low of a cutoff resulting in trivial trimming"
#         ))
#     }
#     attr(paths, "tree") <- tree
#     attr(paths, "align") <- align
#     attr(paths, "class") <- "sitePath"
#     return(paths)
# }

#' @export
print.sitePath <- function(x, ...) {
    cat(length(x), "paths\n")
}

#' @importFrom methods is
checkMatched <- function(tree, align) {
    if (!is(align, "alignment")) {
        stop("align is not class alignment")
    }
    m <- match(tree$tip.label, align$nam)
    if (any(is.na(m))) {
        stop("tree tips and alignment names are not matched")
    } else {
        align <- align$seq[m]
        if (length(unique(nchar(align))) > 1) {
            stop("Sequence lengths are not the same in alignment")
        }
    }
    return(toupper(align))
}

sortSimMatrix <- function(tree, simMatrix) {
    if (!inherits(tree, "phylo")) {
        stop("tree is not class phylo")
    }
    colMatch <- match(tree$tip.label, colnames(simMatrix))
    rowMatch <- match(tree$tip.label, rownames(simMatrix))
    if (is.null(simMatrix)) {
        return(matrix(
            NA,
            ncol = length(tree$tip.label),
            nrow = length(tree$tip.label)
        ))
    } else {
        return(simMatrix[rowMatch, colMatch])
    }
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
