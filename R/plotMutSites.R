#' @rdname plotMutSites
#' @name plotMutSites
#' @title Plot tree and mutation sites
#' @description The mutated sites for each tip in a phylogenetic tree will be
#'   represented as colored dots positioned by their site number.
#' @param x An \code{\link{SNPsites}} object.
#' @param ... Other arguments
#' @return A tree plot with SNP as dots for each tip.
#' @importFrom ggplot2 ggplot geom_point element_blank element_rect
#' @importFrom aplot insert_left
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' plotMutSites(SNPsites(tree))
plotMutSites.SNPsites <- function(x, ...) {
    allSNP <- attr(x, "allSNP")
    snpColors <- vapply(
        X = AA_FULL_NAMES,
        FUN = function(i) {
            AA_COLORS[[i]]
        },
        FUN.VALUE = character(1)
    )
    names(snpColors) <- toupper(names(snpColors))
    p <- ggplot(allSNP, aes(
        x = Pos,
        y = Accession,
        fill = SNP
    )) +
        geom_point(
            shape = 23,
            size = 1,
            stroke = 0,
            color = as.character(snpColors[allSNP[, "SNP"]])
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
    tree <- attr(x, "tree")
    return(insert_left(p, ggtree(tree), 2))
}

#' @export
plotMutSites <- function(x, ...) {
    UseMethod("plotMutSites")
}
