#' @importFrom ape Ntip

#' @rdname parallelSites
#' @title Mutation across multiple phylogenetic lineages
#' @description A site may have mutated on parallel lineages. Mutation can occur
#'   on the same site across the phylogenetic lineages solved by
#'   \code{\link{lineagePath}}. The site will be considered mutated in parallel
#'   if the mutation occurs on the non-overlap part of more than two lineages.
#'   The amino acid/nucleotide before and after the mutation can be allowed
#'   different on different lineages or only the exact same mutations are
#'   considered.
#' @param x A \code{\link{lineagePath}} or a \code{\link{sitesMinEntropy}}
#'   object.
#' @param minSNP The minimum number of mutations to be qualified as parallel on
#'   at least two lineages. The default is 1.
#' @param mutMode The strategy for finding parallel site. The default \code{all}
#'   is to consider any mutation regardless of the amino acid/nucleotide before
#'   and after mutation; Or \code{exact} to force mutation to be the same; Or
#'   \code{pre}/\code{post} to select the site having amino acid/nucleotide
#'   before/after mutation.
#' @param ... The arguments in \code{\link{sitesMinEntropy}}.
#' @return A \code{parallelSites} object
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' paths <- lineagePath(tree)
#' x <- sitesMinEntropy(paths)
#' parallelSites(x)
parallelSites <- function(x, ...) {
    UseMethod("parallelSites")
}

#' @rdname parallelSites
#' @export
parallelSites.lineagePath <- function(x,
                                      minSNP = NULL,
                                      mutMode = c("all", "exact",
                                                  "pre", "post"),
                                      ...) {
    minEntropy <- sitesMinEntropy.lineagePath(x, ...)
    res <- parallelSites.sitesMinEntropy(x = minEntropy,
                                         minSNP = minSNP,
                                         mutMode = mutMode)
    return(res)
}

#' @rdname parallelSites
#' @export
parallelSites.sitesMinEntropy <- function(x,
                                          minSNP = NULL,
                                          mutMode = c("all", "exact",
                                                      "pre", "post"),
                                          ...) {
    paths <- attr(x, "paths")
    align <- attr(paths, "align")
    reference <- attr(paths, "msaNumbering")
    unambiguous <- .unambiguousChars(paths)
    # Set the 'minSNP' constrain
    if (is.null(minSNP)) {
        minSNP <- attr(attr(x, "paths"), "minSize")
    } else if (minSNP != 1) {
        tree <- as.phylo.phyMSAmatched(paths)
        nTips <- Ntip(tree)
        minSNP <- .checkMinEffectiveSize(
            x = minSNP,
            varName = "minSNP",
            totalSize = nTips,
            maxSize = nTips
        )
    }
    # The index of 'mutName' attribute to be used for applying constrains
    mutNameIndex <- switch(
        match.arg(mutMode),
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
        for (siteName in names(segs)) {
            # Convert the site to index of multiple sequence alignment
            site <- reference[as.integer(siteName)]
            # The entropy minimization of a single site on the lineage
            seg <- segs[[siteName]]
            # Find all sporadic mutations by comparing AA/nucleotide of each tip
            # with the fixed one
            for (tips in seg) {
                # The fixed AA/nucleotide of the group
                fixedAA <- attr(tips, "AA")
                # Skip the group if the fixed AA/nucleotide is ambiguous or gap
                if (!fixedAA %in% unambiguous ||
                    length(attr(tips, "aaSummary")) < 2)
                    next
                # The real AA/nucleotide of each tip named with tip name
                tipsAA <- substr(x = align[tips],
                                 start = site,
                                 stop = site)
                # The tips with AA/nucleotide different from the fixed one
                mut <- lapply(
                    X = which(tipsAA != fixedAA & tipsAA %in% unambiguous),
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
            # There has to be at least one fixation on the lineage and at least
            # two of the mutation is neither gap nor ambiguous character
            # Collect all fixation mutations if any for the site
            if (.qualifiedFixation(seg, unambiguous)) {
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
                                         minSNP)
        }
        # Add the mutation result to the collection
        mixedParallel <- c(mixedParallel, list(mixedRef))
        sporadicParallel <- c(sporadicParallel, list(sporadicMut))
        fixationParallel <- c(fixationParallel, list(fixationMut))
    }
    attr(res, "paths") <- paths
    attr(res, "clustersByPath") <- attr(x, "clustersByPath")
    class(res) <- "parallelSites"
    return(res)
}

.collectParallelSites <- function(mutatNodes,
                                  otherNodes,
                                  res,
                                  mutNameIndex,
                                  minSNP) {
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
                # constrain of "minSNP" and filtering "mutMode"
                mutat <- mutat[which(mutatQualifed)]
                # The qualified 'mutName' (indexed according to the constrain
                # mode) on the two lineages
                qualifedAA <- intersect(
                    .qualifiedMutAA(mutat, mutNameIndex, minSNP),
                    .qualifiedMutAA(other, mutNameIndex, minSNP)
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

.qualifiedMutAA <- function(mutat, i, minSNP) {
    # Fixation mutation will ignore 'minSNP' constrain
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
    c(fixationMutAA, names(which(table(sporadicMutAA) >= minSNP)))
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

#' @rdname parallelSites
#' @export
parallelSites.paraFixSites <- function(x, ...) {
    res <- attr(x, "allParaSites")
    return(res)
}
