#' @rdname lineagePath
#' @name lineagePath
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
#' @importFrom ape nodepath getMRCA
#' @export
#' @examples
#' data('zikv_tree')
#' data('zikv_align')
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' lineagePath(tree)
lineagePath.phyMSAmatched <- function(tree,
                                      similarity = NULL,
                                      simMatrix = NULL,
                                      forbidTrivial = TRUE,
                                      ...) {
    x <- .phyMSAmatch(tree)
    tree <- attr(x, "tree")
    # Find total number of tree tips
    nTips <- length(tree[["tip.label"]])
    # Set 'minSNP' with 'nTips' and 'similarity'
    if (is.null(similarity)) {
        minSNP <- nTips * 0.05
        similarity <- 0.05
    } else if (!is.numeric(similarity) || similarity <= 0) {
        stop("\"similarity\" only accepts positive numeric")
    } else if (similarity > 0 && similarity < 0.5) {
        minSNP <- nTips * similarity
    } else if (similarity >= 0.5 && similarity < 1) {
        minSNP <- nTips * (1 - similarity)
        # The reminder of using major SNP to find lineages
        if (forbidTrivial) {
            warning(
                "\"lineagePath\" now uses 'major SNP' to find lineages ",
                "rather than sequence similarities. And the parameter ",
                "\"similarity\" decides the least percentage/number of ",
                "'major SNPs'. The input \"similarity\" of value ",
                similarity,
                "is over 0.5 and the resulting lineage path will ",
                "probably be trivial for other analysis."
            )
        }
    } else if (similarity > 1 && similarity < nTips) {
        minSNP <- ceiling(similarity)
        similarity <- minSNP / nTips
    } else {
        stop(
            "\"similarity\" cannot be greater than total tips. ",
            "And better not be equal to 1."
        )
    }
    align <- attr(x, "align")
    rootNode <- getMRCA(tree, tree[["tip.label"]])
    # Get all lineages using the terminal node found by SNP
    paths <- mergePaths(lapply(
        X = majorSNPtips(align, minSNP),
        FUN = function(tips) {
            nodepath(tree, from = rootNode, to = getMRCA(tree, tips))
        }
    ))
    # Transfer attributes
    paths <- .phyMSAtransfer(paths, x)
    # Set attributes 'similarity' and 'rootNode'
    attr(paths, "similarity") <- similarity
    attr(paths, "rootNode") <- rootNode
    class(paths) <- "lineagePath"
    return(paths)
}

treemerBySite <- function(x, ...) {
    tree <- attr(x, "tree")
    align <- attr(x, "align")
    # Generate the site mapping from reference
    # Exclude the invariant sites
    loci <- attr(x, "loci")
    res <- runTreemerBySite(nodepath(tree), align, loci)
    return(res)
}


#' @export
lineagePath <- function(tree, similarity, ...)
    UseMethod("lineagePath")

#' @export
print.lineagePath <- function(x, ...) {
    cat(
        length(x),
        " lineage paths using ",
        attr(x, "similarity") * 100,
        "% as \"major SNP\" threshold \n",
        sep = ""
    )
}

#' @rdname plotFunctions
#' @title Visualize the results
#' @description The plot function to visualize the return of functions in the
#'   package. Though the function name \code{plot} is used, the plot functions
#'   here do not behave like the generic \code{\link{plot}} function. The
#'   underlying function applies \code{\link{ggplot2}}.
#' @param x Could be a \code{\link{lineagePath}} object,
#'   \code{\link{fixationSites}} object or \code{\link{fixationPath}} object.
#' @param y For \code{\link{lineagePath}} object, it is whether to emphasize the
#'   lineage branches by using thicker line. For a \code{\link{fixationSites}}
#'   or a \code{\link{fixationPath}} object, it is whether to show the fixation
#'   mutation between clusters.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @param ... Other arguments. Since 1.5.4, the function uses
#'   \code{\link{ggtree}} as the base function to make plots so the arguments in
#'   \code{plot.phylo} will no longer work.
#' @return A ggplot object to make the plot. A \code{\link{lineagePath}} object
#'   will be plotted as a tree diagram will be plotted and paths are black solid
#'   line while the trimmed nodes and tips will use grey dashed line. A
#'   \code{\link{fixationSites}} object will be plotted as original phylogenetic
#'   tree marked with fixation substitutions. A \code{\link{fixationPath}}
#'   object will be plotted as a \code{phylo} object. The tips are clustered
#'   according to the fixation sites. The transition of fixation sites will be
#'   plotted as a phylogenetic tree. The length of each branch represents the
#'   number of fixation mutation between two clusters. The name of the tree tips
#'   indicate the number of sequences in the cluster.
#' @importFrom ggtree ggtree
#' @importFrom ggplot2 aes theme ggtitle scale_color_manual scale_size
#' @export
plot.lineagePath <- function(x,
                             y = TRUE,
                             showTips = FALSE,
                             ...) {
    tree <- attr(x, "tree")
    # Get number of ancestral nodes plus tip nodes
    nNodes <- length(tree[["tip.label"]]) + tree[["Nnode"]]
    # Set lineage nodes and non-lineage nodes as separate group
    group <- rep(1, times = nNodes)
    group[unique(unlist(x))] <- 0
    group <- factor(group)
    # Set line size
    if (y) {
        size <- rep(1, times = nNodes)
        size[unique(unlist(x))] <- 2
    } else {
        size <- NULL
    }
    # Tree plot
    p <- ggtree(tree, aes(
        color = group,
        linetype = group,
        size = size
    )) +
        scale_size(range = c(GeomSegment[["default_aes"]][["size"]], 1.5)) +
        scale_color_manual(values = c("black", "gainsboro")) +
        theme(legend.position = "none") +
        ggtitle(attr(x, "similarity"))
    if (showTips) {
        p <- p + geom_tiplab()
    }
    return(p)
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
#' @param makePlot Whether make a plot when return.
#' @return \code{sneakPeek} return the similarity threhold against number of
#'   lineagePath. There will be a simple dot plot between threshold and path
#'   number if \code{makePlot} is TRUE.
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @export
#' @examples
#' sneakPeek(tree, step = 3)
sneakPeek <- function(tree,
                      step = 9,
                      maxPath = NULL,
                      minPath = 1,
                      makePlot = FALSE) {
    x <- .phyMSAmatch(tree)
    tree <- attr(x, "tree")
    # Check 'maxPath'
    if (is.null(maxPath)) {
        maxPath <- length(tree[["tip.label"]]) / 20
    } else if (maxPath <= 0) {
        stop("Invalid \"maxPath\": less than or equal to zero")
    }
    # Check 'minPath'
    if (minPath >= maxPath) {
        stop("Invalid \"minPath\": greater than \"maxPath\"")
    } else if (minPath < 0) {
        stop("Invalid \"minPath\": less than zero")
    }
    # Try every 'similarity'
    similarity <- numeric()
    pathNum <- integer()
    allPaths <- list()
    for (s in seq(from = 0.05,
                  to = 0.01,
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
    res <- data.frame(similarity, pathNum)
    attr(res, "allPaths") <- allPaths
    class(res) <- c(class(res), "sneakPeekedPaths")
    # Combine all plots of the lineages
    if (makePlot) {
        g <- lapply(allPaths, plot, y = FALSE)
        grid.arrange(arrangeGrob(grobs = g))
    }
    return(res)
}

#' @rdname lineagePath
#' @description When used on the return of \code{sneakPeek}, a
#'   \code{lineagePath} with the closest \code{similarity} will be retrived from
#'   the returned value.
#' @export
#' @examples
#' x <- sneakPeek(tree, step = 3)
#' lineagePath(x, similarity = 0.05)
lineagePath.sneakPeekedPaths <- function(tree, similarity, ...) {
    allPaths <- attr(tree, "allPaths")
    allPaths[[which.min(abs(as.numeric(names(allPaths)) - similarity))]]
}
