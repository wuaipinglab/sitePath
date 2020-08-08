#' @rdname parallelSites
#' @name parallelSites
#' @title Mutation in multiple lineages
#' @description A site may have mutated on parallel lineages.
#' @param x A \code{\link{sitesMinEntropy}} object.
#' @param minEffectiveSize The minimum number of mutations to be qualified as
#'   parallel on at least two lineages. The default is 1.
#' @param method The strategy for finding parallel site. The default is to
#'   consider any mutation regardless of the amino acid/nucleotide before and
#'   after mutation; Or the exact same mutation; Or the mutation having amino
#'   acid/nucleotide before or after mutation.
#' @param ... Other arguments.
#' @return A \code{sitesMinEntropy} object
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' paths <- lineagePath(tree)
#' x <- sitesMinEntropy(paths)
#' parallelSites(x)
parallelSites.sitesMinEntropy <- function(x,
                                          minEffectiveSize = NULL,
                                          method = c("all", "exact",
                                                     "pre", "post"),
                                          ...) {
    paths <- attr(x, "paths")
    align <- attr(paths, "align")
    reference <- attr(paths, "msaNumbering")
    # There must be at least two lineages to have mutations
    hasParallelMut <- Reduce("+", lapply(x, lengths))
    hasParallelMut <- names(which(hasParallelMut > 1))
    if (length(hasParallelMut) == 0) {
        stop("There doesn't seem to have any mutation in parallel lineages")
    }
    # Set the 'minEffectiveSize'
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- 1
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    minEffectiveSize <- ceiling(minEffectiveSize)
    # The function to the corresponding method
    mutSelectFunc <- .parallelMutSelectFunc(method)
    # To collect the result by site
    res <- list()
    # Collect the sporadic and fixation mutation on each lineage
    sporadicParallel <- list()
    fixationParallel <- list()
    mixedParallel <- list()
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
                        attr(mutTips, "fixed") <- FALSE
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
                    # Get the tip names
                    mutTips <- names(align[currTips])
                    # Add the tip name
                    attr(mutTips, "mutName") <- mutName
                    attr(mutTips, "fixed") <- TRUE
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
        mixedRef <- c(sporadicMut, fixationMut)
        # Compare mutation result on each pair of lineages
        for (mixedMut in mixedParallel) {
            res <- .collectParallelSites(mixedRef,
                                         mixedMut,
                                         res,
                                         mutSelectFunc,
                                         minEffectiveSize)
        }
        # Add the mutation result to the collection
        mixedParallel <- c(mixedParallel, list(mixedRef))
        sporadicParallel <- c(sporadicParallel, list(sporadicMut))
        fixationParallel <- c(fixationParallel, list(fixationMut))
    }
    attr(res, "paths") <- paths
    class(res) <- "parallelSites"
    return(res)
}

.parallelMutSelectFunc <- function(method = c("all", "exact",
                                              "pre", "post")) {
    resFunc <- switch(
        match.arg(method),
        "all" = function(mutat,
                         other,
                         mutatFix,
                         otherFix,
                         minEffectiveSize) {
            res <- character()
            # The name of the fixation mutation
            mutatFixMut <- .getMutationNames(mutatFix)
            otherFixMut <- .getMutationNames(otherFix)
            # Apply the 'minEffectiveSize' constrain to see if any mutation
            # still qualified on the two lineages
            if ((length(mutat) >= minEffectiveSize ||
                 length(mutatFixMut) > 0) &&
                (length(other) >= minEffectiveSize ||
                 length(otherFixMut) > 0)) {
                # Collect the names of all true parallel mutations
                mutNames <- .getMutationNames(c(mutat, other))
                res <- unique(mutNames)
            }
            return(res)
        },
        "exact" = function(mutat,
                           other,
                           mutatFix,
                           otherFix,
                           minEffectiveSize) {
            # The name of the fixation mutation
            mutatFixMut <- .getMutationNames(mutatFix)
            otherFixMut <- .getMutationNames(otherFix)
            # Summarize the names of parallel mutations on the two lineages
            mutatSummary <- table(.getMutationNames(mutat))
            mutatMutNames <- c(mutatFixMut,
                               names(which(mutatSummary >= minEffectiveSize)))
            otherSummary <- table(.getMutationNames(other))
            otherMutNames <- c(otherFixMut,
                               names(which(otherSummary >= minEffectiveSize)))
            # Apply the 'minEffectiveSize' constrain: a mutation has to be
            # qualified on both lineages
            res <- intersect(mutatMutNames, otherMutNames)
            return(res)
        },
        "pre" = function(mutat,
                         other,
                         mutatFix,
                         otherFix,
                         minEffectiveSize) {
            # The qualified amino acid/nucleotide before the mutation on the two
            # lineages
            qualifedAA <- intersect(
                .qualifiedMutAA(mutat, mutatFix, 1, minEffectiveSize),
                .qualifiedMutAA(other, otherFix, 1, minEffectiveSize)
            )
            res <- c(.selectMutByAA(mutat, 1, qualifedAA),
                     .selectMutByAA(other, 1, qualifedAA))
            return(res)
        },
        "post" = function(mutat,
                          other,
                          mutatFix,
                          otherFix,
                          minEffectiveSize) {
            # The qualified amino acid/nucleotide before the mutation on the two
            # lineages
            qualifedAA <- intersect(
                .qualifiedMutAA(mutat, mutatFix, 3, minEffectiveSize),
                .qualifiedMutAA(other, otherFix, 3, minEffectiveSize)
            )
            res <- c(.selectMutByAA(mutat, 3, qualifedAA),
                     .selectMutByAA(other, 3, qualifedAA))
            return(res)
        }
    )
    return(resFunc)
}

.getMutationNames <- function(mutat) {
    vapply(
        X = mutat,
        FUN = function(mut) {
            paste0(attr(mut, "mutName"), collapse = "")
        },
        FUN.VALUE = character(1)
    )
}

.qualifiedMutAA <- function(mutat, mutatFix, i, minEffectiveSize) {
    mutAAsummary <- table(vapply(
        X = mutat,
        FUN = function(mut) {
            attr(mut, "mutName")[i]
        },
        FUN.VALUE = character(1)
    ))
    fixedMutAA <- vapply(
        X = mutatFix,
        FUN = function(mut) {
            attr(mut, "mutName")[i]
        },
        FUN.VALUE = character(1)
    )
    c(fixedMutAA, names(which(mutAAsummary >= minEffectiveSize)))
}

.selectMutByAA <- function(mutat, i, qualifiedAA) {
    res <- vapply(
        X = mutat,
        FUN = function(mut) {
            mutName <- attr(mut, "mutName")
            if (mutName[i] %in% qualifiedAA) {
                return(paste0(mutName, collapse = ""))
            }
            return(NA_character_)
        },
        FUN.VALUE = character(1)
    )
    res[which(!is.na(res))]
}

.collectParallelSites <- function(mutatNodes,
                                  otherNodes,
                                  res,
                                  mutSelectFunc,
                                  minEffectiveSize) {
    allMutatNodes <- names(mutatNodes)
    allOtherNodes <- names(otherNodes)
    mutatSitesFull <- .groupMutationsBySites(mutatNodes,
                                             allMutatNodes)
    otherSitesFull <- .groupMutationsBySites(otherNodes,
                                             allOtherNodes)
    # Use the ancestral node to to remove overlapped part
    mutatDiff <- setdiff(allMutatNodes, allOtherNodes)
    otherDiff <- setdiff(allOtherNodes, allMutatNodes)
    # There has to have non-overlap part for both lineage
    if (length(mutatDiff) != 0 && length(otherDiff) != 0) {
        mutatSites <- .groupMutationsBySites(mutatNodes, mutatDiff)
        otherSites <- .groupMutationsBySites(otherNodes, otherDiff)
        candidateSites <- intersect(names(mutatSites),
                                    names(otherSites))
        # Check if the candidate sites are qualified because the uniqueness
        # judged by ancestral node could be incorrect if entropy minimization
        # result is off on the two lineages
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
                # Apply the constrains to get the mutations that meet the
                # constrain of "minEffectiveSize" and filtering "method"
                mutat <- mutat[which(mutatQualifed)]
                # Fixation mutation will ignore 'minEffectiveSize' constrain
                mutatFix <- Filter(function(mut) {
                    attr(mut, "fixed")
                }, mutat)
                otherFix <- Filter(function(mut) {
                    attr(mut, "fixed")
                }, other)
                qualifiedMut <- mutSelectFunc(mutat,
                                              other,
                                              mutatFix,
                                              otherFix,
                                              minEffectiveSize)
                if (length(qualifiedMut) > 0) {
                    toAdd <- c(
                        Filter(function(mut) {
                            paste0(attr(mut, "mutName"),
                                   collapse = "") %in% qualifiedMut
                        }, mutatSitesFull[[site]]),
                        Filter(function(mut) {
                            paste0(attr(mut, "mutName"),
                                   collapse = "") %in% qualifiedMut
                        }, otherSitesFull[[site]])
                    )
                    res[[site]] <- c(res[[site]], list(toAdd))
                    class(res[[site]]) <- "sitePara"
                }
            }
        }
    }
    return(res)
}

.groupMutationsBySites <- function(nodeGrouped, nodeNames) {
    nodeGrouped <- nodeGrouped[nodeNames]
    res <- list()
    for (node in names(nodeGrouped)) {
        for (mut in nodeGrouped[[node]]) {
            site <- attr(mut, "mutName")[2]
            mut <- list(mut)
            names(mut) <- node
            res[[site]] <- c(res[[site]], mut)
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
    cat("This is a 'parallelSites' object.\n\nResult for",
        length(attr(x, "paths")),
        "paths:\n\n")
    if (length(x) == 0) {
        cat("No parallel site found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(x, "reference")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}

#' @export
print.sitePara <- function(x, ...) {
    class(x) <- NULL
    print(x)
}
