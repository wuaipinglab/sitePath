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
    # The site with mutation on each lineage
    sporadicParallelSites <- list()
    fixationParallelSites <- list()
    # Iterate entropy minimization result for each lineage. This part is to
    # remove the duplicate mutations on the overlapped part of the lineages
    for (segs in x) {
        # To collect mutation on the current lineage
        sporadicMut <- list()
        fixationMut <- list()
        # To collect all sites with mutation on the current lineage
        sporadicMutSites <- character()
        fixationMutSites <- character()
        # The existing mutation tip/node of the previous lineages
        existingSporadic <- unlist(lapply(sporadicParallel, names))
        existingFixation <- unlist(lapply(fixationParallel, names))
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
                # Add the mutation info to sporadic mutation collection of the
                # lineage if the the tip is not found in the previous lineages
                mut <- mut[setdiff(names(mut), existingSporadic)]
                # Add the mutation info to sporadic mutation collection of the
                # lineage
                for (mutNode in names(mut)) {
                    mutName <- list(mut[[mutNode]])
                    if (mutNode %in% names(sporadicMut)) {
                        sporadicMut[[mutNode]] <- c(sporadicMut[[mutNode]],
                                                    mutName)
                    } else {
                        sporadicMut[[mutNode]] <- mutName
                    }
                    sporadicMutSites <- c(sporadicMutSites,
                                          siteName)
                }
            }
            # Collect all fixation mutations if any for the site
            if (length(seg) >= 2) {
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
                    # Add the mutation info to fixation mutation collection of
                    # the lineage if the the same group of tips is not found in
                    # the previous lineages
                    if (!mutNode %in% existingFixation) {
                        mutName <- list(mutName)
                        if (mutNode %in% names(fixationMut)) {
                            fixationMut[[mutNode]] <- c(fixationMut[[mutNode]],
                                                        mutName)
                        } else {
                            fixationMut[[mutNode]] <- mutName
                        }
                        fixationMutSites <- c(fixationMutSites,
                                              siteName)
                    }
                }
            }
        }
        sporadicParallel <- c(sporadicParallel, list(sporadicMut))
        fixationParallel <- c(fixationParallel, list(fixationMut))
        sporadicParallelSites <- c(sporadicParallelSites,
                                   list(unique(sporadicMutSites)))
        fixationParallelSites <- c(fixationParallelSites,
                                   list(unique(fixationMutSites)))
    }
    return(fixationParallelSites)
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
