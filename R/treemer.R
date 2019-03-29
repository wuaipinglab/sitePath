#' @rdname pre-assessment
#' @name pre-assessment
#' @title Things can be done before the analysis
#' @description
#' \code{similarityMatrix} calculates similarity between aligned sequences
#' The similarity matrix can be used in \code{\link{groupTips}}
#' or \code{\link{lineagePath}}
#' @param tree The return from \code{\link{addMSA}} function
#' @examples
#' data("zikv_tree")
#' data("zikv_align")
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' simMatrix <- similarityMatrix(tree)
#' @return
#' \code{similarityMatrix} returns a diagonal matrix of
#' similarity between sequences
#' @importFrom methods is
#' @export
similarityMatrix <- function(tree) {
    if (!inherits(tree, "phylo")) {
        stop("\"tree\" is not class phylo")
    } else if (is.null(attr(tree, "alignment"))) {
        stop("No alignment found in \"tree\"")
    }
    sim <- getSimilarityMatrix(attr(tree, "alignment"))
    dimnames(sim) <- list(tree$tip.label, tree$tip.label)
    return(sim)
}

sortSimMatrix <- function(tree, simMatrix) {
    if (!inherits(tree, "phylo")) {
        stop("\"tree\" is not class phylo")
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
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
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
#' \code{lineagePath} finds the lineages of a phylogenetic tree providing
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
#' lineagePath(tree, 0.996)
#' @return path represent by node number
#' @export
lineagePath <- function(tree,
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
    class(paths) <- "lineagePath"
    return(paths)
}

#' @export
print.lineagePath <- function(x, ...) {
    cat(length(x), "paths\n")
}

#' @rdname pre-assessment
#' @description
#' \code{sneakPeek} is intended to plot similarity as a threshold
#' against number of output lineagePath. This plot is intended to give user
#' a feel about how many sitePaths they should expect from
#' the similarity threshold. The number of lineagePath should not
#' be too many or too few. The result excludes where the number of lineagePath
#' is greater than number of tips divided by 20 or self-defined maxPath.
#' The zero lineagePath result will also be excluded
#' @param step
#' the similarity window for calculating and ploting. To better
#' see the impact of threshold on path number. This is preferably
#' specified. The default is one 50th of the difference between 1
#' and minimal pairwise sequence similarity.
#' @param maxPath
#' maximum number of path to return show in the plot. The number of path
#' in the raw tree can be far greater than trimmed tree. To better
#' see the impact of threshold on path number. This is preferably
#' specified. The default is one 20th of tree tip number.
#' @param minPath
#' minimum number of path to return show in the plot. To better
#' see the impact of threshold on path number. This is preferably
#' specified. The default is 1.
#' @param makePlot whether make a dot plot when return
#' @examples
#' sneakPeek(tree)
#' @return
#' \code{sneakPeek} return the similarity threhold against number of lineagePath.
#' There will be a simple dot plot between threshold and path number if
#' \code{makePlot} is TRUE.
#' @importFrom methods is
#' @importFrom graphics plot
#' @export
sneakPeek <- function(tree,
                      step = NULL,
                      maxPath = NULL,
                      minPath = 1,
                      makePlot = TRUE) {
    simMatrix <- similarityMatrix(tree)
    minSim <- min(simMatrix)
    if (is.null(step)) {
        step <- round(minSim - 1, 3) / 50
    }
    if (is.null(maxPath)) {
        maxPath <- length(tree$tip.label) / 20
    } else if (maxPath <= 0) {
        stop("Invalid \"maxPath\": less than or equal to zero")
    }
    if (minPath >= maxPath) {
        stop("Invalid \"minPath\": greater than \"maxPath\"")
    } else if (minPath < 0) {
        stop("Invalid \"minPath\": less than zero")
    }
    similarity <- numeric(0)
    pathNum <- integer(0)
    for (s in seq(1, minSim, step)) {
        paths <- lineagePath(
            tree,
            similarity = s,
            simMatrix = simMatrix,
            forbidTrivial = FALSE
        )
        if (maxPath < length(paths)) {
            next
        } else if (length(paths) <= minPath) {
            break
        }
        similarity <- c(similarity, s)
        pathNum <- c(pathNum, length(paths))
    }
    if (makePlot) {
        plot(similarity, pathNum)
    }
    return(data.frame(similarity, pathNum))
}
