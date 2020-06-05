#' @rdname lineagePath
#' @name lineagePath
#' @title Resolving lineage paths using SNP
#' @description \code{lineagePath} finds the lineages of a phylogenetic tree
#'   providing the corresponding sequence alignment. This is done by finding
#'   'major SNPs' which usually accumulate along the evolutionary pathways.
#' @param tree The return from \code{\link{addMSA}} function
#' @param similarity This decides how minor SNPs are to remove. If provided as
#'   fraction between 0 and 1, then the minimum number of SNP will be total tips
#'   times \code{similariy}. If provided as integer greater than 1, the minimum
#'   number will be \code{similariy}. The default \code{similariy} is 0.05 for
#'   \code{lineagePath}.
#' @param simMatrix Deprecated and will not have effect.
#' @param forbidTrivial Does not allow trivial trimming.
#' @return path represent by node number
#' @importFrom ape nodepath
#' @export
#' @examples
#' data('zikv_tree')
#' data('zikv_align')
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' lineagePath(tree)
lineagePath <- function(tree,
                        similarity = NULL,
                        simMatrix = NULL,
                        forbidTrivial = TRUE) {
    nTips <- length(tree[["tip.label"]])
    if (is.null(similarity)) {
        minSNP <- nTips * 0.05
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

#' @rdname lineagePath
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
#' @return \code{sneakPeek} return the similarity threhold against number of
#'   lineagePath. There will be a simple dot plot between threshold and path
#'   number if \code{makePlot} is TRUE.
#' @export
#' @examples
#' sneakPeek(tree)
sneakPeek <- function(tree,
                      step = 10,
                      maxPath = NULL,
                      minPath = 1,
                      makePlot = TRUE) {
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
    similarity <- numeric()
    pathNum <- integer()
    allPaths <- list()
    for (s in seq(from = 0.05, to = 0.01, length.out = step)) {
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
        allPaths[[as.character(similarity)]] <- paths
    }
    res <- data.frame(similarity, pathNum)
    attr(res, "allPaths") <- allPaths
    class(res) <- "lineagePathSP"
    if (makePlot) {
        plot(similarity, pathNum)
    }
    return(res)
}
