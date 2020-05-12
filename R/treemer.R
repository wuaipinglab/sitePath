#' @rdname treemer
#' @name treemer
#' @title Topology-dependent tree trimming
#' @description \code{groupTips} uses sequence similarity to group tree tips.
#'   Members in a group are always constrained to share the same ancestral node.
#'   Similarity between two tips is derived from their multiple sequence
#'   alignment. The site will not be counted into total length if both are gap.
#'   Similarity is calculated as number of matched divided by the corrected
#'   total length. So far the detection of divergence is based on one simple
#'   rule: the miminal pairwise similarity. The two branches are decided to be
#'   divergent if the similarity is lower than the threshold. (Other more
#'   statistical approaches such as Kolmogorov-Smirnov Tests among pair-wise
#'   distance could be introduced in the future)
#' @param tree The return from \code{\link{addMSA}} function
#' @param similarity Similarity threshold for tree trimming in \code{groupTips}.
#'   If not provided, the mean similarity substract standard deviation of all
#'   sequences will be used. And for \code{lineagePath}, this decides how minor
#'   SNPs are to remove. If provided as fraction between 0 and 1, then the
#'   minimum number of SNP will be total tips times \code{similariy}. If
#'   provided as integer greater than 1, the minimum number will be
#'   \code{similariy}. The default \code{similariy} is 0.1 for
#'   \code{lineagePath}.
#' @param simMatrix A diagonal matrix of similarities for each pair of
#'   sequences. This parameter will not have effect in the function
#'   \code{lineagePath}.
#' @param forbidTrivial Does not allow trivial trimming
#' @param tipnames If return as tipnames
#' @importFrom ape nodepath
#' @importFrom stats sd
#' @examples
#' data('zikv_tree')
#' data('zikv_align')
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' groupTips(tree, 0.996)
#' @return grouping of tips
#' @export
groupTips <- function(tree,
                      similarity = NULL,
                      simMatrix = NULL,
                      forbidTrivial = TRUE,
                      tipnames = TRUE) {
    if (is.null(simMatrix)) {
        simMatrix <- similarityMatrix(tree)
    } else {
        colMatch <- match(tree$tip.label, colnames(simMatrix))
        rowMatch <- match(tree$tip.label, rownames(simMatrix))
        if (is.null(simMatrix)) {
            simMatrix <- matrix(
                NA,
                ncol = length(tree$tip.label),
                nrow = length(tree$tip.label)
            )
        } else {
            simMatrix <- simMatrix[rowMatch, colMatch]
        }
    }
    if (is.null(similarity)) {
        simDist <- simMatrix[upper.tri(simMatrix)]
        similarity <- mean(simDist) - sd(simDist)
    }
    align <- attr(tree, "align")
    if (is.null(align)) {
        stop("No alignment found in \"tree\"")
    }
    grouping <- runTreemer(
        tipPaths = nodepath(tree),
        alignedSeqs = align,
        simMatrixInput = simMatrix,
        similarity = similarity,
        getTips = TRUE
    )
    if (length(grouping) == 1 && forbidTrivial) {
        warning(
            "\"similarity\" ",
            similarity,
            " is too low of a cutoff resulting in trivial trimming"
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
#' @description \code{lineagePath} finds the lineages of a phylogenetic tree
#'   providing the corresponding sequence alignment. This is done by finding
#'   'major SNPs' which usually accumulate along the evolutionary pathways. are
#'   added.
#' @importFrom ape nodepath
#' @examples
#' lineagePath(tree)
#' @return path represent by node number
#' @export
lineagePath <- function(tree,
                        similarity = NULL,
                        simMatrix = NULL,
                        forbidTrivial = TRUE) {
    nTips <- length(tree[["tip.label"]])
    if (is.null(similarity)) {
        minSNP <- nTips * 0.1
    } else if (!is.numeric(similarity) || similarity <= 0) {
        stop("\"similarity\" only accepts positive numeric")
    } else if (similarity > 0 && similarity < 1) {
        minSNP <- nTips * similarity
    } else if (similarity > 1 && similarity < nTips) {
        minSNP <- ceiling(similarity)
    } else {
        stop(
            "\"similarity\" cannot be greater than total tips. ",
            "And better not be equal to 1."
        )
    }
    align <- attr(tree, "align")
    attr(tree, "align") <- NULL
    if (is.null(align)) {
        stop("No alignment found in \"tree\"")
    }
    paths <- mergePaths(lapply(
        X = majorSNPtips(align, minSNP),
        FUN = function(tips) {
            nodepath(tree, from = nTips + 1, to = getMRCA(tree, tips))
        }
    ))
    if (length(paths) == 0 && forbidTrivial) {
        warning(
            "\"lineagePath\" now uses 'major SNP' to find lineages ",
            "rather than sequence similarities. And the parameter ",
            "\"similarity\" decides the least percentage/number of ",
            "'major SNPs'. The input \"similarity\" of value ",
            similarity,
            "is too high leading to no usable 'major SNP' and ",
            "hence no lineage path found."
        )
    }
    attr(paths, "reference") <- attr(tree, "reference")
    attr(tree, "reference") <- NULL
    attr(paths, "tree") <- tree
    attr(paths, "align") <- align
    attr(paths, "rootNode") <- getMRCA(tree, tree[["tip.label"]])
    class(paths) <- "lineagePath"
    return(paths)
}

#' @export
print.lineagePath <- function(x, ...) {
    cat(length(x), "paths\n")
}

#' @rdname pre-assessment
#' @name pre-assessment
#' @title Things can be done before the analysis
#' @description \code{similarityMatrix} calculates similarity between aligned
#'   sequences The similarity matrix can be used in \code{\link{groupTips}} or
#'   \code{\link{lineagePath}}
#' @param tree The return from \code{\link{addMSA}} function
#' @examples
#' data('zikv_tree')
#' data('zikv_align')
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' simMatrix <- similarityMatrix(tree)
#' @return \code{similarityMatrix} returns a diagonal matrix of similarity
#'   between sequences
#' @importFrom methods is
#' @export
similarityMatrix <- function(tree) {
    if (!inherits(tree, "phylo")) {
        stop("\"tree\" is not class phylo")
    } else if (is.null(attr(tree, "align"))) {
        stop("No alignment found in \"tree\"")
    }
    sim <- getSimilarityMatrix(attr(tree, "align"))
    dimnames(sim) <- list(tree$tip.label, tree$tip.label)
    return(sim)
}

#' @rdname pre-assessment
#' @description \code{sneakPeek} is intended to plot 'similarity' (actually the
#'   least percentage of 'major SNP') as a threshold against number of output
#'   lineagePath. This plot is intended to give user a rought view about how
#'   many lineages they could expect from the 'similarity' threshold in the
#'   function \code{\link{lineagePath}}. The number of lineagePath is preferably
#'   not be too many or too few. The result excludes where the number of
#'   lineagePath is greater than number of tips divided by 20 or user-defined
#'   maxPath. The zero lineagePath result will also be excluded.
#' @param step the 'similarity' window for calculating and ploting. To better
#'   see the impact of threshold on path number. The default is 10.
#' @param maxPath maximum number of path to return show in the plot. The number
#'   of path in the raw tree can be far greater than trimmed tree. To better see
#'   the impact of threshold on path number. This is preferably specified. The
#'   default is one 20th of tree tip number.
#' @param minPath minimum number of path to return show in the plot. To better
#'   see the impact of threshold on path number. The default is 1.
#' @param makePlot whether make a dot plot when return
#' @examples
#' sneakPeek(tree)
#' @return \code{sneakPeek} return the similarity threhold against number of
#'   lineagePath. There will be a simple dot plot between threshold and path
#'   number if \code{makePlot} is TRUE.
#' @importFrom methods is
#' @importFrom graphics plot
#' @export
sneakPeek <- function(tree,
                      step = 10,
                      maxPath = NULL,
                      minPath = 1,
                      makePlot = FALSE) {
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
    for (s in seq(from = 0.3, to = 0.01, length.out = step)) {
        paths <- lineagePath(tree,
                             similarity = s,
                             forbidTrivial = FALSE)
        if (length(paths) > maxPath) {
            break
        } else if (length(paths) <= minPath) {
            similarity <- c(similarity, s)
            pathNum <- c(pathNum, length(paths))
            next
        }
        similarity <- c(similarity, s)
        pathNum <- c(pathNum, length(paths))
    }
    if (makePlot) {
        plot(similarity, pathNum)
    }
    return(data.frame(similarity, pathNum))
}
