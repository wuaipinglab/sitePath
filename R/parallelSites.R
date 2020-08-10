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
    # The index of 'mutName' attribute to be used for applying constrains
    mutNameIndex <- switch(
        match.arg(method),
        "all" = 2,
        "exact" = 4,
        "pre" = 1,
        "post" = 3
    )
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
                        # The tip name to return
                        mutTips <- names(tipsAA[i])
                        # The mutation info used later
                        mutName <- c(fixedAA, siteName, tipsAA[[i]])
                        mutName <- c(mutName,
                                     paste0(mutName, collapse = ""))
                        # Add the tip name
                        attr(mutTips, "mutName") <- mutName
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
                    # The mutation info used later
                    mutName <- c(prevAA, siteName, currAA)
                    mutName <- c(mutName,
                                 paste0(mutName, collapse = ""))
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
                                         mutNameIndex,
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

.qualifiedMutAA <- function(mutat, i, minEffectiveSize) {
    # Fixation mutation will ignore 'minEffectiveSize' constrain
    fixationMutAA <- character()
    sporadicMutAA <- character()
    for (mutTips in mutat) {
        # Get the indexed 'mutName' for deciding parallelity
        mutAA <- attr(mutTips, "mutName")[i]
        if (attr(mutTips, "fixed")) {
            fixationMutAA <- c(fixationMutAA, mutAA)
        } else {
            sporadicMutAA <- c(sporadicMutAA, mutAA)
        }
    }
    # Summarize the number of each amino acid/nucleotide of each sporadic
    # mutation
    c(fixationMutAA, names(which(
        table(sporadicMutAA) >= minEffectiveSize
    )))
}

.selectMutByAA <- function(mutat, i, qualifiedAA) {
    # Select the mutation name according to the result of qualification
    res <- vapply(
        X = mutat,
        FUN = function(mut) {
            mutName <- attr(mut, "mutName")
            if (mutName[i] %in% qualifiedAA) {
                return(mutName[4])
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
                                  mutNameIndex,
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
                # The qualified 'mutName' (indexed according to the constrain
                # mode) on the two lineages
                qualifedAA <- intersect(
                    .qualifiedMutAA(mutat, mutNameIndex, minEffectiveSize),
                    .qualifiedMutAA(other, mutNameIndex, minEffectiveSize)
                )
                # Collect all mutation names that qualify the constrain
                qualifiedMut <- c(
                    .selectMutByAA(mutat, mutNameIndex, qualifedAA),
                    .selectMutByAA(other, mutNameIndex, qualifedAA)
                )
                if (length(qualifiedMut) > 0) {
                    # Only the tips with qualified mutation will be added to the
                    # result. The result still keeps the mutation separate of
                    # the two lineages
                    toAdd <- c(
                        Filter(function(mut) {
                            attr(mut, "mutName")[4] %in% qualifiedMut
                        }, mutatSitesFull[[site]]),
                        Filter(function(mut) {
                            attr(mut, "mutName")[4] %in% qualifiedMut
                        }, otherSitesFull[[site]])
                    )
                    # The result of a site will contain qualified tips on the
                    # pair of lineages
                    res[[site]] <- c(res[[site]], list(toAdd))
                    # Re-assign or assign the class
                    attr(res[[site]], "site") <- site
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
    cat(
        "This is a 'sitePara' object.\n\nSite",
        attr(x, "site"),
        "may have parallel mutation on",
        length(x),
        "pair of paths:\n\n"
    )
    mutSummary <- table(vapply(
        X = extractTips.sitePara(x),
        FUN = function(mutTips) {
            attr(mutTips, "mutName")[4]
        },
        FUN.VALUE = character(1)
    ))
    mutInfo <- character()
    for (mutName in names(mutSummary)) {
        mutInfo <-
            c(mutInfo, paste0(mutName, "(", mutSummary[[mutName]], ")"))
    }
    cat(
        paste0(mutInfo, collapse = ", "),
        "\n\nIn the bracket are the number of tips",
        "involved in the mutation\n"
    )
}
