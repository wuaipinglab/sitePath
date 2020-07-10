#' @rdname parallelSites
#' @title Amino acid in parallel lineages
#' @description A site might be mutated to the same amino acid in multiple
#'   lineages.
#' @param x A \code{\link{lineagePath}} object. And amino acid sequence must be
#'   provided.
#' @param minFreq Minimum number of an amino acid to be reported and contribute
#'   to make a parallel site. The default is 1.
#' @return A \code{parallelSites} objects.
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' paths <- lineagePath(tree)
#' parallelSites(paths)
parallelSites <- function(x,
                          minFreq = NULL,
                          refSeqName = NULL) {
    if (length(x) == 1) {
        stop("There should be at least two lineages but \"x\" only has 1.")
    }
    x <- .phyMSAmatch(x)
    # Set default 'minSNP' value
    if (is.null(minFreq)) {
        minFreq <- 1
    }
    align <- attr(x, "align")
    msaNumbering <- attr(x, "msaNumbering")
    allSNP <- .findSNPsites(align, msaNumbering, refSeqName)
    # Get the divergent nodes
    divNodes <- divergentNode(paths)
    # The tips and matching sequence
    nodeAlign <- .tipSeqsAlongPathNodes(paths = paths,
                                        divNodes = divNodes)
    # The tips in all lineages
    tipNames <- as.phylo.phyMSAmatched(x)[["tip.label"]]
    pathsTips <- lapply(x, function(p) {
        lapply(as.character(p[which(!p %in% divNodes)]), function(n) {
            tipNames[as.integer(names(nodeAlign[[n]]))]
        })
    })
    pathsTips <- .mergeClusters(pathsTips)
    pathsTips <- lapply(pathsTips, unlist)
    # The parallel site should at least appears in two different lineages
    allSNP <- split(x = allSNP,
                    f = allSNP[, c("Pos", "SNP")],
                    drop = TRUE)
    isParallel <- which(vapply(
        X = allSNP,
        FUN = function(snp) {
            tips <- snp[, "Accession"]
            sum(vapply(
                X = pathsTips,
                FUN = function(p) {
                    sum(tips %in% p) >= minFreq
                },
                FUN.VALUE = logical(1)
            )) >= 2
        },
        FUN.VALUE = logical(1)
    ))
    if (length(isParallel) == 0) {
        stop("No parallel site found using \"minFreq\" ", minFreq)
    }
    allSNP <- do.call(rbind, allSNP[isParallel])
    rownames(allSNP) <- NULL
    res <- sort(unique(allSNP[["Pos"]]))
    attr(res, "allSNP") <- allSNP
    # Transfer attributes
    res <- .phyMSAtransfer(res, x)
    class(res) <- "parallelSites"
    return(res)
}

#' @export
print.parallelSites <- function(x, ...) {
    x <- as.integer(x)
    print(x)
}
