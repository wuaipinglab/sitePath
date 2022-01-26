#' @rdname paraFixSites
#' @title The fixation sites with mutation on parallel lineage
#' @description The operation between the results of \code{\link{fixationSites}}
#'   and \code{\link{parallelSites}}.
#' @param x A \code{lineagePath} object returned from \code{\link{lineagePath}}
#'   function.
#' @param ... further arguments passed to or from other methods.
#' @return A \code{paraFixSites} object.
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' paraFixSites(zikv_tree_reduced, alignment = zikv_align_reduced)
paraFixSites <- function(x, ...) {
    UseMethod("paraFixSites")
}

#' @rdname paraFixSites
#' @param alignment An \code{alignment} object. This commonly can be from
#'   sequence parsing function in the \code{\link{seqinr}} package. Sequence
#'   names in the alignment should include all \code{tip.label} in the tree
#' @param seqType The type of the sequence in the alignment file. The default is
#'   "AA" for amino acid. The other options are "DNA" and "RNA".
#' @param Nmin The parameter for identifying phylogenetic pathway using SNP. If
#'   provided as fraction between 0 and 1, then the minimum number of SNP will
#'   be total tips times \code{Nmin}. If provided as integer greater than 1, the
#'   minimum number will be \code{Nmin}.
#' @param reference Name of reference for site numbering. The name has to be one
#'   of the sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar The character to indicate gap. The numbering will skip the
#'   \code{gapChar} for the reference sequence.
#' @param minSkipSize The minimum number of tips to have gap or ambiguous amino
#'   acid/nucleotide for a site to be ignored in other analysis. This will not
#'   affect the numbering. The default is 0.8.
#' @export
paraFixSites.phylo <- function(x,
                               alignment = NULL,
                               seqType = c("AA", "DNA", "RNA"),
                               Nmin = NULL,
                               reference = NULL,
                               gapChar = "-",
                               minSkipSize = NULL,
                               ...) {
    paths <- lineagePath.phylo(
        tree = x,
        alignment = alignment,
        seqType = seqType,
        Nmin = Nmin,
        reference = reference,
        gapChar = gapChar,
        minSkipSize = minSkipSize
    )
    res <- paraFixSites.lineagePath(paths, ...)
    return(res)
}

#' @rdname paraFixSites
#' @export
paraFixSites.treedata <- function(x, ...) {
    tree <- as.phylo(x)
    res <- paraFixSites.phylo(tree)
    return(res)
}

#' @rdname paraFixSites
#' @param minEffectiveSize The minimum number of tips in a group.
#' @param searchDepth The function uses heuristic search but the termination of
#'   the search cannot be intrinsically decided. \code{searchDepth} is needed to
#'   tell the search when to stop.
#' @param method The strategy for predicting the fixation. The basic approach is
#'   entropy minimization and can be achieved by adding or removing fixation
#'   point, or by comparing the two.
#' @export
paraFixSites.lineagePath <- function(x,
                                     minEffectiveSize = NULL,
                                     searchDepth = 1,
                                     method = c("compare", "insert", "delete"),
                                     ...) {
    minEntropy <- sitesMinEntropy.lineagePath(x,
                                              minEffectiveSize,
                                              searchDepth,
                                              method)
    res <- paraFixSites.sitesMinEntropy(minEntropy, ...)
    return(res)
}

#' @rdname paraFixSites
#' @param category Could be \code{parallelOnly}, \code{fixationOnly},
#'   \code{intersect} or \code{union}.
#' @param minSNP The minimum number of mutations to be qualified as parallel on
#'   at least two lineages. The default is 1.
#' @param mutMode The strategy for finding parallel site. The default \code{all}
#'   is to consider any mutation regardless of the amino acid/nucleotide before
#'   and after mutation; Or \code{exact} to force mutation to be the same; Or
#'   \code{pre}/\code{post} to select the site having amino acid/nucleotide
#'   before/after mutation.
#' @export
paraFixSites.sitesMinEntropy <- function(x,
                                         category = c("intersect",
                                                      "union",
                                                      "parallelOnly",
                                                      "fixationOnly"),
                                         minSNP = NULL,
                                         mutMode = c("all", "exact",
                                                     "pre", "post"),
                                         ...) {
    # All fixation sites and site names
    fixSites <- fixationSites.sitesMinEntropy(x)
    fixSiteNames <- allSitesName.fixationSites(fixSites)
    # All parallel sites and site names
    paraSites <- parallelSites.sitesMinEntropy(x, minSNP, mutMode)
    paraSiteNames <- allSitesName.parallelSites(paraSites)
    # Derive the sites based on the specified category
    sites <- switch(
        match.arg(category),
        "intersect" = intersect(fixSiteNames, paraSiteNames),
        "union" = union(fixSiteNames, paraSiteNames),
        "parallelOnly" = setdiff(paraSiteNames, fixSiteNames),
        "fixationOnly" = setdiff(fixSiteNames, paraSiteNames)
    )
    # Site names as result
    res <- sort(as.integer(sites))
    attr(res, "allFixSites") <- fixSites
    attr(res, "allParaSites") <- paraSites
    # Subset of the fixation and parallel sites
    fixSites <- fixSites[intersect(sites, fixSiteNames)]
    paraSites <- paraSites[intersect(sites, paraSiteNames)]
    # Set fixation sites attribute if any
    if (length(fixSites)) {
        # Set 'paths' and 'clustersByPath' attributes
        attr(fixSites, "paths") <- attr(x, "paths")
        attr(fixSites, "clustersByPath") <-
            attr(x, "clustersByPath")
        class(fixSites) <- "fixationSites"
        attr(res, "fixSites") <- fixSites
    }
    # Set 'allSNP' (to represent parallel sites) attribute if any
    if (length(paraSites)) {
        allSNP <- unlist(paraSites,
                         recursive = FALSE,
                         use.names = FALSE)
        allSNP <- unlist(allSNP,
                         recursive = FALSE,
                         use.names = FALSE)
        allSNP <- do.call(rbind, lapply(allSNP, function(tips) {
            Pos <- as.integer(rep(attr(tips, "mutName")[2], length(tips)))
            SNP <- rep(attr(tips, "mutName")[3], length(tips))
            data.frame("Accession" = tips,
                       "Pos" = Pos,
                       "SNP" = SNP)
        }))
        attr(res, "allSNP") <- allSNP
    }
    attr(res, "paths") <- attr(x, "paths")
    class(res) <- "paraFixSites"
    return(res)
}
