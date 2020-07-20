#' @rdname SNPsites
#' @title Finding sites with variation
#' @description Single nucleotide polymorphism (SNP) in the whole package refers
#'   to variation of amino acid. \code{SNPsite} will try to find SNP in the
#'   multiple sequence alignment. A reference sequence and gap character may be
#'   specified to number the site.
#' @param tree A \code{\link{phyMSAmatched}} object.
#' @param minSNP Minimum number of a mutation to be a SNP. The default is 10th
#'   of the total tree tips.
#' @return A \code{SNPsites} object.
#' @importFrom stats complete.cases
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' SNPsites(tree)
SNPsites <- function(tree, minSNP = NULL) {
    x <- .phyMSAmatch(tree)
    tree <- attr(x, "tree")
    # Set default 'minSNP' value
    if (is.null(minSNP)) {
        minSNP <- length(tree[["tip.label"]]) / 10
    }
    align <- attr(x, "align")
    msaNumbering <- attr(x, "msaNumbering")
    allSNP <- .findSNPsites(align, msaNumbering, NULL)
    # Calculate the frequency of each mutation/SNP
    snpSummary <- as.data.frame(table(allSNP[["Pos"]],
                                      allSNP[["SNP"]]))
    allSNP <- merge(
        x = allSNP,
        y = snpSummary,
        by.x = c("Pos", "SNP"),
        by.y = c("Var1", "Var2"),
        all.x = TRUE
    )
    # Filter out low frequency mutation/SNP
    allSNP <- allSNP[which(allSNP[, "Freq"] >= minSNP),
                     c("Accession", "Pos", "SNP")]
    rownames(allSNP) <- NULL
    # Extract all the qualified sites as 'res' to be compatible with the return
    # of previous version
    res <- sort(unique(allSNP[["Pos"]]))
    attr(res, "allSNP") <- allSNP
    # Transfer attributes
    res <- .phyMSAtransfer(res, x)
    class(res) <- "SNPsites"
    return(res)
}

.findSNPsites <- function(align, msaNumbering, refSeqName) {
    # Find SNP for each tree tip by comparing with the consensus sequence
    if (is.null(refSeqName)) {
        # Find the major SNP of each site as the consensus sequence
        referenceSeq <- vapply(
            X = msaNumbering,
            FUN = function(s) {
                aaSummary <- tableAA(align, s - 1)
                # The amino acid/nucleotide having the most appearance
                names(aaSummary)[which.max(aaSummary)]
            },
            FUN.VALUE = character(1)
        )
        align <- strsplit(x = align, split = "")
    } else {
        align <- strsplit(x = align, split = "")
        referenceSeq <- align[[refSeqName]]
    }
    allSNP <- lapply(seq_along(referenceSeq), function(site) {
        snp <- vapply(
            X = align,
            FUN = "[[",
            i = msaNumbering[site],
            FUN.VALUE = character(1)
        )
        snp <- snp[which(snp != referenceSeq[[site]] & snp != '-')]
        res <- data.frame(
            "Accession" = names(snp),
            "Pos" = rep(site, length(snp)),
            "SNP" = snp
        )
        return(res)
    })
    allSNP <- do.call(rbind, allSNP)
    return(allSNP)
}

.phyMSAtransfer <- function(receive, give) {
    attr(receive, "align") <- attr(give, "align")
    attr(receive, "tree") <- attr(give, "tree")
    attr(receive, "msaNumbering") <- attr(give, "msaNumbering")
    attr(receive, "reference") <- attr(give, "reference")
    return(receive)
}

#' @export
print.SNPsites <- function(x, ...) {
    x <- as.integer(x)
    print(x)
}

#' @rdname plotMutSites
#' @name plotMutSites
#' @title Plot tree and mutation sites
#' @description The mutated sites for each tip in a phylogenetic tree will be
#'   represented as colored dots positioned by their site number.
#' @param x An \code{\link{SNPsites}} object.
#' @param showTips Whether to plot the tip labels. The default is \code{FALSE}.
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
plotMutSites.SNPsites <- function(x, showTips = FALSE, ...) {
    allSNP <- attr(x, "allSNP")
    # Specify the color of mutations by pre-defined color set.
    snpColors <- vapply(
        X = AA_FULL_NAMES,
        FUN = function(i) {
            AA_COLORS[[i]]
        },
        FUN.VALUE = character(1)
    )
    names(snpColors) <- toupper(names(snpColors))
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
    # Use 'ggtree' to make tree plot
    treePlot <- ggtree(attr(x, "tree"))
    if (showTips) {
        treePlot <- treePlot + geom_tiplab()
    }
    # Combine the two plots and return
    return(insert_left(snpPlot, treePlot, 2))
}

#' @export
plotMutSites <- function(x, ...) {
    UseMethod("plotMutSites")
}
