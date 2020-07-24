#' @rdname parallelSites
#' @name parallelSites
#' @title Mutation in multiple lineages
#' @description A site may have mutated on parallel lineages.
#' @param x A \code{\link{sitesMinEntropy}} object.
#' @param minSNP Minimum number of tips to have mutations on a lineage.
#' @return A \code{sitesMinEntropy} object
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' paths <- lineagePath(tree)
#' x <- sitesMinEntropy(paths)
#' parallelSites(x)
parallelSites.sitesMinEntropy <- function(x, minSNP, ...) {
    paths <- attr(x, "paths")
    align <- attr(paths, "align")
    reference <- attr(paths, "msaNumbering")
    hasParallelMut <- Reduce("+", lapply(x, lengths))
    hasParallelMut <- names(which(hasParallelMut > 1))
    if (length(hasParallelMut) == 0) {
        stop("There doesn't seem to have any mutation in parallel lineages")
    }
    fixationParallel <- list()
    sporadicParallel <- list()
    for (segs in x) {
        fixationMut <- list()
        sporadicMut <- list()
        for (siteName in hasParallelMut) {
            site <- reference[as.integer(siteName)]
            seg <- segs[siteName]
            # Collect all non-fixation mutations
            mut <- list()
            for (tips in seg) {
                fixedAA <- attr(tips, "AA")
                tipsAA <- substr(x = align[tips],
                                 start = site,
                                 stop = site)
                mutTips <-
                    names(which(tipsAA != fixedAA & tipsAA != '-'))
                mutName <- c(fixedAA, siteName, tipsAA[mutTips])
                mut <- c(mut, lapply(mutName, function(m) {
                    attr(m, "tip") <- t
                    return(m)
                }))
            }
            sporadicMut <- c(sporadicMut, mut)
            # Collect all fixation mutations
            if (length(seg) >= 2) {
                mut <- list()
                for (i in seq_along(seg)[-1]) {
                    prevTips <- seg[i - 1]
                    currTips <- seg[i]
                    prevAA <- attr(prevTips, "AA")
                    currAA <- attr(currTips, "AA")
                    mutNode <- attr(currTips, "node")
                    mutName <- c(prevAA, siteName, currAA)
                    currTipsAA <- substr(x = align[currTips],
                                         start = site,
                                         stop = site)
                    tips <- names(which(currTipsAA == currAA))
                    attr(mutName, "tips") <- tips
                    mut[[mutNode]] <- mutName
                }
                fixationMut <- c(fixationMut, mut)
            }
        }
    }
    allParallel <- data.frame(
        "Accession" = character(),
        "Pos" = integer(),
        "From" = character(),
        "To" = character(),
        "Fixed" = logical()
    )
    res <- sort(unique(allParallel[["Pos"]]))
    attr(res, "allParallel") <- allParallel
    attr(res, "paths") <- paths
    class(res) <- "parallelSites"
    return(res)
}

#' @export
parallelSites <- function(x, ...) {
    UseMethod("parallelSites")
}

#' @export
print.parallelSites <- function(x, ...) {
    cat("This is a 'parallelSites' object.\n")
}
