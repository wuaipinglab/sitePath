#' @importFrom stats sd quantile
#' @importFrom ape nodepath getMRCA node.depth.edgelength Ntip
#' @importFrom gridExtra arrangeGrob grid.arrange

#' @rdname lineagePath
#' @title Resolving lineage paths using SNP
#' @description \code{lineagePath} finds the lineages of a phylogenetic tree
#'   providing the corresponding sequence alignment. This is done by finding
#'   'major SNPs' which usually accumulate along the evolutionary pathways.
#' @param tree The return from \code{\link{addMSA}} or \code{sneakPeek}
#'   function.
#' @param similarity This decides how minor SNPs are to remove. If provided as
#'   fraction between 0 and 1, then the minimum number of SNP will be total tips
#'   times \code{similariy}. If provided as integer greater than 1, the minimum
#'   number will be \code{similariy}. The default \code{similariy} is 0.05 for
#'   \code{lineagePath}.
#' @param simMatrix Deprecated and will not have effect.
#' @param forbidTrivial Does not allow trivial trimming.
#' @param ... Other arguments.
#' @return Lineage path represent by node number.
#' @export
#' @examples
#' data('zikv_tree')
#' data('zikv_align')
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' lineagePath(tree)
lineagePath <- function(tree, similarity, ...) {
    UseMethod("lineagePath")
}

#' @rdname lineagePath
#' @export
lineagePath.phyMSAmatched <- function(tree,
                                      similarity = NULL,
                                      simMatrix = NULL,
                                      forbidTrivial = TRUE,
                                      ...) {
    x <- .phyMSAmatch(tree)
    tree <- attr(x, "tree")
    # Find total number of tree tips
    nTips <- Ntip(tree)
    # The range of results using different 'minSize'
    rangeOfResults <- attr(x, "rangeOfResults")
    maxSize <-
        max(vapply(rangeOfResults, attr, integer(1), "minSize"))
    # Set the number of lineage paths
    if (is.null(similarity)) {
        pathNums <- lengths(rangeOfResults)
        pathNums <- pathNums[which(pathNums != 1)]
        fourthQuantile <- quantile(pathNums)[[4]]
        paths <-
            rangeOfResults[[which.min(abs(pathNums - fourthQuantile))]]
        minSNP <- attr(paths, "minSize")
        similarity <- minSNP / nTips
    } else {
        if (!is.numeric(similarity) || similarity <= 0) {
            stop("\"similarity\" only accepts positive numeric")
        } else if (similarity > 1) {
            minSNP <- ceiling(similarity)
        } else {
            minSNP <- nTips * similarity
        }
        if (minSNP > maxSize) {
            if (forbidTrivial) {
                minSNP <- maxSize
            } else {
                stop(
                    "\"similarity\" cannot be greater than total tips. ",
                    "And better not be equal to 1."
                )
            }
        }
        similarity <- minSNP / nTips
        paths <- rangeOfResults[[which.min(vapply(
            X = rangeOfResults,
            FUN = function(ps) {
                abs(attr(ps, "minSize") - minSNP)
            },
            FUN.VALUE = numeric(1)
        ))]]
    }
    # Transfer attributes
    attributes(paths) <- attributes(x)
    # Set attributes 'similarity' and 'minSNP'
    attr(paths, "similarity") <- similarity
    attr(paths, "minSize") <- minSNP
    class(paths) <- c("lineagePath", "phyMSAmatched")
    snp <- SNPsites(paths)
    attr(paths, "allSNP") <- attr(snp, "allSNP")
    return(paths)
}

.dropLessDiverged <- function(paths, tree) {
    allTipsNum <- Ntip(tree)
    toKeep <- vapply(
        X = seq_along(paths),
        FUN = function(i) {
            refPath <- paths[[i]]
            # All tips at the terminal
            endTipsNum <-
                length(.childrenTips(tree, tail(refPath, 1)))
            # The diverged part from all other lineage
            divPath <- setdiff(refPath, unlist(paths[-i]))
            divTipsNum <- length(.childrenTips(tree, divPath[1]))
            # If the diverged part (terminal tips excluded) has larger
            # proportion
            div <- (divTipsNum - endTipsNum) / divTipsNum
            all <- endTipsNum / allTipsNum
            div > all
        },
        FUN.VALUE = logical(1)
    )
    paths <- paths[which(toKeep)]
    return(paths)
}

.mergeDivergedPaths <- function(tree, zValue, rootNode) {
    edgeLength <- node.depth.edgelength(tree)
    res <- nodepath(tree)
    # Iterate each path to see if diverged from all other paths
    i <- 1
    while (i <= length(res)) {
        # The target lineage path excluding the terminal node
        refPath <- res[[i]]
        endNode <- refPath[length(refPath)]
        # Assume the current terminal tips are to keep
        toDrop <- FALSE
        # Compare with rest of the paths to decide if to keep
        for (p in res[-i]) {
            if (edgeLength[endNode] > edgeLength[p[length(p)]]) {
                # Always keep the longer path
                next
            }
            # Pair up the two paths for comparison
            pathPair <- list(refPath, p)
            # The intersection of the two paths
            commonNodes <- intersect(refPath, p)
            if (length(commonNodes) == 1 &&
                commonNodes == rootNode) {
                # The two paths diverged at the root
                divNodes <- rootNode
            } else {
                divNodes <- divergentNode(pathPair)
            }
            attr(pathPair, "tree") <- tree
            # Fake the 'align' attribute, so the function will work even without
            # the actual sequence
            attr(pathPair, "align") <- tree[["tip.label"]]
            # The tips and the corresponding ancestral node. Because the tree
            # was forced to be bifurcated, there won't be any tips on the
            # divergent node
            pathNodeTips <-
                .tipSeqsAlongPathNodes(pathPair, divNodes)
            # The genetic distance of all tips to the path
            allDist <- vapply(
                X = names(pathNodeTips),
                FUN = function(n) {
                    an <- as.integer(n)
                    tips <- pathNodeTips[[n]]
                    tipsDist <- edgeLength[tips] - edgeLength[an]
                    mean(tipsDist)
                },
                FUN.VALUE = numeric(1)
            )
            divThreshold <- mean(allDist) + sd(allDist) * zValue
            # The diverged part of the current path compared with the other path
            divDist <- edgeLength[endNode] - edgeLength[divNodes]
            if (divDist < divThreshold) {
                toDrop <- TRUE
                break
            }
        }
        if (toDrop) {
            res <- res[-i]
        } else {
            i <- i + 1
        }
    }
    return(res)
}

#' @rdname lineagePath
#' @description \code{sneakPeek} is intended to plot 'similarity' (actually the
#'   least percentage of 'major SNP') as a threshold against number of output
#'   lineagePath. This plot is intended to give user a rough view about how many
#'   lineages they could expect from the 'similarity' threshold in the function
#'   \code{\link{lineagePath}}. The number of lineagePath is preferably not be
#'   too many or too few. The result excludes where the number of lineagePath is
#'   greater than number of tips divided by 20 or user-defined maxPath. The zero
#'   lineagePath result will also be excluded.
#' @param step the 'similarity' window for calculating and plotting. To better
#'   see the impact of threshold on path number. The default is 10.
#' @param maxPath maximum number of path to return show in the plot. The number
#'   of path in the raw tree can be far greater than trimmed tree. To better see
#'   the impact of threshold on path number. This is preferably specified. The
#'   default is one 20th of tree tip number.
#' @param minPath minimum number of path to return show in the plot. To better
#'   see the impact of threshold on path number. The default is 1.
#' @param makePlot Whether make a plot when return.
#' @return \code{sneakPeek} return the similarity threhold against number of
#'   lineagePath. There will be a simple dot plot between threshold and path
#'   number if \code{makePlot} is TRUE.
#' @export
#' @examples
#' sneakPeek(tree, step = 3)
sneakPeek <- function(tree,
                      step = 9,
                      maxPath = NULL,
                      minPath = 0,
                      makePlot = TRUE) {
    x <- .phyMSAmatch(tree)
    tree <- attr(x, "tree")
    nTips <- Ntip(tree)
    # Check 'maxPath'
    if (is.null(maxPath)) {
        maxPath <- nTips
    } else if (maxPath <= 0) {
        stop("Invalid \"maxPath\": less than or equal to zero")
    }
    # Check 'minPath'
    if (minPath >= maxPath) {
        stop("Invalid \"minPath\": greater than \"maxPath\"")
    } else if (minPath < 0) {
        stop("Invalid \"minPath\": less than zero")
    }
    # The range of results using different 'minSize'
    rangeOfResults <- attr(x, "rangeOfResults")
    rangeOfResults <-
        rangeOfResults[which(lengths(rangeOfResults) != 1)]
    sizeRange <-
        range(vapply(rangeOfResults, attr, integer(1), "minSize"))
    # Try every 'similarity'
    similarity <- numeric()
    pathNum <- integer()
    allPaths <- list()
    for (s in seq(from = sizeRange[1],
                  to = sizeRange[2],
                  length.out = step)) {
        paths <- lineagePath(x,
                             similarity = s,
                             forbidTrivial = FALSE)
        if (length(paths) > maxPath) {
            # Terminate when hit 'maxPath'
            break
        } else if (length(paths) <= minPath) {
            # Go to next 'similarity' when less than 'minPath'
            next
        }
        similarity <- c(similarity, s)
        pathNum <- c(pathNum, length(paths))
        allPaths[[as.character(s)]] <- paths
    }
    res <- as.data.frame(t(vapply(
        X = attr(x, "rangeOfResults"),
        FUN = function(i) {
            c("similarity" = attr(i, "minSize") / nTips,
              "pathNum" = length(i))
        },
        FUN.VALUE = numeric(2)
    )))
    attr(res, "phyMSAmatched") <- x
    class(res) <- c(class(res), "sneakPeekedPaths")
    # Combine all plots of the lineages
    if (makePlot) {
        g <- lapply(allPaths, function(path) {
            plot(path) + ggtitle(
                label = paste0(
                    attr(path, "minSize"),
                    " (",
                    round(attr(path, "similarity") * 100, 2),
                    "%)"
                ),
                subtitle = paste0(length(path), " path(s)")
            )
        })
        grid.arrange(arrangeGrob(grobs = g))
    }
    return(res)
}

#' @rdname lineagePath
#' @description When used on the return of \code{sneakPeek}, a
#'   \code{lineagePath} with the closest \code{similarity} will be retrieved
#'   from the returned value.
#' @export
#' @examples
#' x <- sneakPeek(tree, step = 3)
#' lineagePath(x, similarity = 0.05)
lineagePath.sneakPeekedPaths <- function(tree, similarity, ...) {
    tr <- attr(tree, "phyMSAmatched")
    lineagePath.phyMSAmatched(tr, similarity, ...)
}
