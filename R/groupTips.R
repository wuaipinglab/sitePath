#' @rdname similarityMatrix
#' @title Similarity between sequences
#' @description Get similarity between aligned sequences with gap ignored.
#' @param tree The return from \code{\link{addMSA}} function.
#' @return A diagonal matrix of similarity between sequences.
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' simMatrix <- similarityMatrix(tree)
similarityMatrix <- function(tree) {
    x <- .phyMSAmatch(tree)
    sim <- getSimilarityMatrix(attr(x, "align"))
    tree <- attr(x, "tree")
    dimnames(sim) <- list(tree[["tip.label"]], tree[["tip.label"]])
    return(sim)
}

#' @rdname groupTips
#' @name groupTips
#' @title The grouping of tree tips
#' @description \code{groupTips} uses sequence similarity to group tree tips.
#'   Members in a group are always constrained to share the same ancestral node.
#'   Similarity between two tips is derived from their multiple sequence
#'   alignment. The site will not be counted into total length if both are gap.
#'   Similarity is calculated as number of matched divided by the corrected
#'   total length. So far the detection of divergence is based on one simple
#'   rule: the minimal pairwise similarity. The two branches are decided to be
#'   divergent if the similarity is lower than the threshold.
#' @param tree The return from \code{\link{addMSA}}, \code{\link{lineagePath}}
#'   or other functions in the \code{sitePath} package.
#' @param similarity Similarity threshold for tree trimming in \code{groupTips}.
#'   If not provided, the mean similarity subtract standard deviation of all
#'   sequences will be used.
#' @param simMatrix A diagonal matrix of similarities for each pair of
#'   sequences.
#' @param forbidTrivial Does not allow trivial trimming.
#' @param tipnames If return tips as integer or tip names.
#' @param ... Other arguments.
#' @return \code{groupTips} returns grouping of tips.
#' @export
#' @importFrom ape nodepath
#' @importFrom stats sd
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' groupTips(tree)
groupTips.phyMSAmatched <- function(tree,
                                    similarity = NULL,
                                    simMatrix = NULL,
                                    forbidTrivial = TRUE,
                                    tipnames = TRUE,
                                    ...) {
    x <- .phyMSAmatch(tree)
    tree <- attr(x, "tree")
    # Set 'simMatrix'
    colMatch <- match(tree[["tip.label"]], colnames(simMatrix))
    rowMatch <- match(tree[["tip.label"]], rownames(simMatrix))
    if (is.null(simMatrix)) {
        simMatrix <- matrix(NA,
                            ncol = length(tree[["tip.label"]]),
                            nrow = length(tree[["tip.label"]]))
    } else {
        simMatrix <- simMatrix[rowMatch, colMatch]
    }
    # Set default 'similarity'
    if (is.null(similarity)) {
        simDist <- simMatrix[upper.tri(simMatrix)]
        similarity <- mean(simDist) - sd(simDist)
    }
    # Use treemer
    grouping <- runTreemer(
        tipPaths = nodepath(tree),
        alignedSeqs = attr(x, "align"),
        simMatrixInput = simMatrix,
        similarity = similarity,
        getTips = TRUE
    )
    # If the result is trivial
    if (length(grouping) == 1 && forbidTrivial) {
        warning(
            "The input \"similarity\" ",
            similarity,
            " is too low of a cutoff resulting in trivial trimming"
        )
    }
    # Return with tip names or not
    if (tipnames) {
        return(lapply(grouping, function(g) {
            tree[["tip.label"]][g]
        }))
    } else {
        return(grouping)
    }
}

#' @rdname groupTips
#' @export
groupTips.fixationSites <- function(tree, tipnames = TRUE, ...) {
    grp <- fixationPath.fixationSites(tree, minEffectiveSize = 0)
    res <- groupTips.fixationPath(grp, tipnames = tipnames)
    return(res)
}

#' @rdname groupTips
#' @export
groupTips.fixationPath <- function(tree, tipnames = TRUE, ...) {
    x <- tree
    groupName <- names(x)
    tree <- attr(x, "tree")
    attributes(x) <- NULL
    if (tipnames) {
        res <- lapply(x, function(tips) {
            attributes(tips) <- NULL
            return(tree[["tip.label"]][tips])
        })
    } else {
        res <- lapply(x, function(tips) {
            attributes(tips) <- NULL
            return(tips)
        })
    }
    names(res) <- groupName
    return(res)
}

#' @export
groupTips <- function(tree, ...) {
    UseMethod("groupTips")
}