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
    # There must be at least two lineages to have mutations
    hasParallelMut <- Reduce("+", lapply(x, lengths))
    hasParallelMut <- names(which(hasParallelMut > 1))
    if (length(hasParallelMut) == 0) {
        stop("There doesn't seem to have any mutation in parallel lineages")
    }
    # Collect the sporadic and fixation mutation on each lineage
    sporadicParallel <- list()
    fixationParallel <- list()
    # Iterate entropy minimization result for each lineage
    for (segs in x) {
        # To collect mutation on the current lineage
        sporadicMut <- list()
        fixationMut <- list()
        # The site have mutated on at least two lineages
        for (siteName in hasParallelMut) {
            # Convert the site to index of multiple sequence alignment
            site <- reference[as.integer(siteName)]
            # The entropy minimization of a single site on the lineage
            seg <- segs[[siteName]]
            # Find all sporadic mutations by comparing AA/nucleotide of each tip
            # with the fixed one
            for (tips in seg) {
                # The fixed AA/nucleotide of the group
                fixedAA <- attr(tips, "AA")
                # The real AA/nucleotide of each tip named with tip name
                tipsAA <- substr(x = align[tips],
                                 start = site,
                                 stop = site)
                # The tips with AA/nucleotide different from the fixed one
                mut <- lapply(
                    X = which(tipsAA != fixedAA & tipsAA != '-'),
                    FUN = function(i) {
                        # Return the mutation info
                        mutName <- c(fixedAA, siteName, tipsAA[[i]])
                        # Add the tip name
                        attr(mutName, "tips") <- names(tipsAA[i])
                        return(mutName)
                    }
                )
                # Add the mutation info to all mutation collection of the
                # lineage
                sporadicMut <- c(sporadicMut, mut)
            }
            # Collect all fixation mutations if any for the site
            if (length(seg) >= 2) {
                mut <- list()
                # Compare the fixed AA/nucleotide between two adjacent groups
                for (i in seq_along(seg)[-1]) {
                    prevTips <- seg[[i - 1]]
                    currTips <- seg[[i]]
                    prevAA <- attr(prevTips, "AA")
                    currAA <- attr(currTips, "AA")
                    mutNode <- attr(currTips, "node")
                    mutName <- c(prevAA, siteName, currAA)
                    # The real AA/nucleotide of each tip named with tip name
                    currTipsAA <- substr(x = align[currTips],
                                         start = site,
                                         stop = site)
                    # Find the tips actually have the same AA/nucleotide as the
                    # fixed one
                    tips <- names(which(currTipsAA == currAA))
                    # Add the tip name
                    attr(mutName, "tips") <- tips
                    mut[[mutNode]] <- mutName
                }
                fixationMut <- c(fixationMut, mut)
            }
        }
        sporadicParallel <- c(sporadicParallel, list(sporadicMut))
        fixationParallel <- c(fixationParallel, list(fixationMut))
    }
    # allParallel <- data.frame(
    #     "Accession" = character(),
    #     "Pos" = integer(),
    #     "From" = character(),
    #     "To" = character(),
    #     "Fixed" = logical()
    # )
    # res <- sort(unique(allParallel[["Pos"]]))
    # attr(res, "allParallel") <- allParallel
    # attr(res, "paths") <- paths
    # class(res) <- "parallelSites"
    # return(res)
}

#' @export
parallelSites <- function(x, ...) {
    UseMethod("parallelSites")
}

#' @export
print.parallelSites <- function(x, ...) {
    cat("This is a 'parallelSites' object.\n")
}
