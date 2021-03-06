#' @importFrom ggplot2 ggplot geom_point theme element_blank element_rect
#' @importFrom aplot insert_left
#' @importFrom ggtree geom_tiplab

#' @rdname plotMutSites
#' @title Plot tree and mutation sites
#' @description The mutated sites for each tip in a phylogenetic tree will be
#'   represented as colored dots positioned by their site number.
#' @param x An \code{\link{SNPsites}} object.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
#' @param ... Other arguments
#' @return A tree plot with SNP as dots for each tip.
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' plotMutSites(SNPsites(tree))
plotMutSites <- function(x, ...) {
    UseMethod("plotMutSites")
}

#' @rdname plotMutSites
#' @export
plotMutSites.SNPsites <- function(x, showTips = FALSE, ...) {
    phyMSAmatched <- attr(x, "phyMSAmatched")
    # Use 'ggtree' to make tree plot
    treePlot <- ggtree(as.phylo(phyMSAmatched))
    if (showTips) {
        treePlot <- treePlot + geom_tiplab()
    }
    allSNP <- attr(x, "allSNP")
    seqType <- attr(phyMSAmatched, "seqType")
    snpPlot <- .createSNPplot(allSNP, seqType)
    # Combine the two plots and return
    return(insert_left(snpPlot, treePlot, 2))
}

.createSNPplot <- function(allSNP, seqType) {
    # Specify the color of mutations by pre-defined color set.
    snpColors <- .siteColorScheme(seqType)
    # Use 'ggplot' to make SNP plot as dots
    snpPlot <- ggplot(allSNP, aes(x = Pos,
                                  y = Accession,
                                  fill = SNP)) +
        geom_point(
            shape = 23,
            size = 1,
            stroke = 0,
            color = snpColors[allSNP[, "SNP"]]
        ) +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = "white"),
            legend.position = "none"
        )
    return(snpPlot)
}

#' @rdname plotMutSites
#' @export
plotMutSites.lineagePath <- function(x, ...) {
    snpPlot <- .createSNPplot(allSNP = attr(x, "allSNP"),
                              seqType = attr(x, "seqType"))
    treePlot <- plot(x, ...)
    return(insert_left(snpPlot, treePlot, 2))
}

#' @rdname plotMutSites
#' @export
plotMutSites.fixationSites <- function(x, ...) {
    paths <- attr(x, "paths")
    snpPlot <- .createSNPplot(allSNP = attr(paths, "allSNP"),
                              seqType = attr(paths, "seqType"))
    treePlot <- plot(x, ...)
    return(insert_left(snpPlot, treePlot, 2))
}
