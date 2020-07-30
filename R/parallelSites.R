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
    # Iterate entropy minimization result for each lineage. This part is to
    # remove the duplicate mutations on the overlapped part of the lineages
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
                        mutTips <- names(tipsAA[i])
                        # Add the tip name
                        attr(mutTips, "mutName") <- c(fixedAA,
                                                      siteName,
                                                      tipsAA[[i]])
                        return(mutTips)
                    }
                )
                # Add the mutation info to sporadic mutation collection of the
                # lineage
                for (mutNode in names(mut)) {
                    mutTips <- list(mut[[mutNode]])
                    if (mutNode %in% names(sporadicMut)) {
                        sporadicMut[[mutNode]] <- c(sporadicMut[[mutNode]],
                                                    mutTips)
                    } else {
                        sporadicMut[[mutNode]] <- mutTips
                    }
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
                    mutTips <- names(which(currTipsAA == currAA))
                    # Add the tip name
                    attr(mutTips, "mutName") <- mutName
                    # Add the mutation info to fixation mutation collection
                    mutTips <- list(mutTips)
                    if (mutNode %in% names(fixationMut)) {
                        fixationMut[[mutNode]] <- c(fixationMut[[mutNode]],
                                                    mutTips)
                    } else {
                        fixationMut[[mutNode]] <- mutTips
                    }
                }
            }
        }
        sporadicParallel <- c(sporadicParallel, list(sporadicMut))
        fixationParallel <- c(fixationParallel, list(fixationMut))
    }
    res <- integer()
    pathsNum <- length(paths)
    # Compare fixation result on each pair of lineages
    for (i in seq_len(pathsNum)[-pathsNum]) {
        sporadicRef <- sporadicParallel[[i]]
        fixationRef <- fixationParallel[[i]]
        for (j in seq(i + 1, pathsNum)) {
            # Both lineage should each have their own unique sporadic mutations
            sporadicMut <- sporadicParallel[[j]]
            res <- c(res, as.integer(
                .parallelSitesOnTwoPaths(sporadicRef, sporadicMut)
            ))
            # Both lineage should each have their own unique fixation mutations
            fixationMut <- fixationParallel[[j]]
            # res <- c(res, as.integer(
            #     .parallelSitesOnTwoPaths(fixationRef, fixationMut)
            # ))
        }
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
    attr(res, "paths") <- paths
    class(res) <- "parallelSites"
    return(res)
}

.parallelSitesOnTwoPaths <- function(mutatNodes, otherNodes) {
    res <- NULL
    mutatDiff <- setdiff(names(mutatNodes), names(otherNodes))
    otherDiff <- setdiff(names(otherNodes), names(mutatNodes))
    if (length(mutatDiff) != 0 && length(otherDiff) != 0) {
        mutatSites <- .groupMutationsBySites(mutatNodes, mutatDiff)
        otherSites <- .groupMutationsBySites(otherNodes, otherDiff)
        candidateSites <- intersect(names(mutatSites),
                                    names(otherSites))
        # Check if the candidate sites are qualified because the uniqueness
        # judged by ancestral node could be incorrect if entropy minimization
        # result is off on the two lineages
        res <- character()
        for (site in candidateSites) {
            # All candidate groups of tip(s) on the two lineages for the site
            mutat <- mutatSites[[site]]
            other <- otherSites[[site]]
            # Assume all tip groups in one set has actual parallel site
            mutatQualifed <- rep(TRUE, length(mutat))
            # Iterate one set of groups and compare with the other group set
            for (i in seq_along(mutat)) {
                mTips <- mutat[[i]]
                for (j in seq_along(other)) {
                    oTips <- other[[j]]
                    if (length(intersect(mTips, oTips)) != 0) {
                        # Disqualify the group if it has overlap with the group
                        # in another set
                        mutatQualifed[i] <- FALSE
                        # The group with overlap in another set is dropped
                        other <- other[-j]
                        # There is no point looking at the rest
                        break
                    }
                }
            }
            # The site is qualified only when both sets pass the check
            if (any(mutatQualifed) && length(other) != 0) {
                res <- c(res, site)
            }
        }
    }
    return(res)
}

.groupMutationsBySites <- function(nodeGrouped, siteNames) {
    nodeGrouped <- nodeGrouped[siteNames]
    res <- list()
    for (node in nodeGrouped) {
        for (mut in node) {
            site <- attr(mut, "mutName")[2]
            res[[site]] <- c(res[[site]], list(mut))
        }
    }
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
