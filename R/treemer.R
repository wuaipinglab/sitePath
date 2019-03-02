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
#' @param tree The return from \code{\link{addMSA}} function
#' @param similarity
#' Similarity threshold for tree trimming. If not provided, an average
#' value of similarity among all sequences will be used.
#' @param simMatrix S diagonal matrix of similarity between sequences
#' @param forbidTrivial Does not allow trivial trimming
#' @param tipnames If return as tipnames
#' @importFrom ape nodepath
#' @examples
#' data("zikv_tree")
#' data("zikv_align")
#' tree <- addMSA(zikv_tree, seqs = zikv_align)
#' groupTips(tree, 0.996)
#' @return grouping of tips
#' @export
groupTips <- function(tree,
                      similarity = NULL,
                      simMatrix = NULL,
                      forbidTrivial = TRUE,
                      tipnames = TRUE) {
    if (is.null(similarity)) {
        simMatrix <- similarityMatrix(tree)
        similarity <- mean(simMatrix)
    } else {
        simMatrix <- sortSimMatrix(tree, simMatrix)
    }
    align <- attr(tree, "alignment")
    if (is.null(align)) {
        stop("No alignment found in \"tree\"")
    }
    grouping <- runTreemer(nodepath(tree),
                           align,
                           simMatrix,
                           similarity,
                           TRUE)
    if (length(grouping) == 1 && forbidTrivial) {
        warning(
            paste(
                "\"similarity\"",
                similarity,
                "is too low of a cutoff resulting in trivial trimming"
            )
        )
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
#' sitePath(tree, 0.996)
#' @return path represent by node number
#' @export
sitePath <- function(tree,
                     similarity = NULL,
                     simMatrix = NULL,
                     forbidTrivial = TRUE) {
    if (is.null(similarity)) {
        simMatrix <- similarityMatrix(tree)
        similarity <- mean(simMatrix)
    } else {
        simMatrix <- sortSimMatrix(tree, simMatrix)
    }
    align <- attr(tree, "alignment")
    if (is.null(align)) {
        stop("No alignment found in \"tree\"")
    }
    # nodepath after trimming
    trimmedPaths <-
        unique(runTreemer(nodepath(tree), align, simMatrix, similarity, FALSE))
    # get the bifurcated pre-terminal nodes and their path to the root
    # those paths are the so-called sitePaths (isolated)
    paths <- lapply(trimmedPaths, function(p)
        p[seq_len(length(p) - 1)])
    paths <-
        unique(paths[which(duplicated(paths) & lengths(paths) > 1)])
    if (length(paths) == 0 && forbidTrivial) {
        warning(
            paste(
                "\"similarity\"",
                similarity,
                "is too low of a cutoff resulting in trivial trimming"
            )
        )
    }
    attr(paths, "tree") <- tree
    attr(paths, "align") <- align
    class(paths) <- "sitePath"
    return(paths)
}

#' @export
print.sitePath <- function(x, ...) {
    cat(length(x), "paths\n")
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
