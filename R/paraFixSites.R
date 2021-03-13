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
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' paraFixSites(lineagePath(tree))
paraFixSites <- function(x, ...) {
    UseMethod("paraFixSites")
}

#' @rdname paraFixSites
#' @param category Could be \code{parallelOnly}, \code{fixationOnly},
#'   \code{intersect} or \code{union}.
#' @param minEffectiveSize The minimum number of tips in a group.
#' @param searchDepth The function uses heuristic search but the termination of
#'   the search cannot be intrinsically decided. \code{searchDepth} is needed to
#'   tell the search when to stop.
#' @param method The strategy for predicting the fixation. The basic approach is
#'   entropy minimization and can be achieved by adding or removing fixation
#'   point, or by comparing the two.
#' @export
paraFixSites.lineagePath <- function(x,
                                     category = c("intersect", "union",
                                                  "parallelOnly", "fixationOnly"),
                                     minEffectiveSize = NULL,
                                     searchDepth = 1,
                                     method = c("compare", "insert", "delete"),
                                     minSNP = NULL,
                                     mutMode = c("all", "exact",
                                                 "pre", "post"),
                                     ...) {
    minEntropy <- sitesMinEntropy.lineagePath(paths,
                                              minEffectiveSize,
                                              searchDepth,
                                              method)
    res <- paraFixSites.sitesMinEntropy(minEntropy,
                                        category,
                                        minSNP,
                                        mutMode)
    return(res)
}

#' @rdname paraFixSites
#' @export
paraFixSites.sitesMinEntropy <- function(x,
                                         category = c("intersect", "union",
                                                      "parallelOnly", "fixationOnly"),
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
    # Subset of the fixation and parallel sites
    fixSites <- fixSites[intersect(sites, fixSiteNames)]
    paraSites <- paraSites[intersect(sites, paraSiteNames)]
    # Site names as result
    res <- sort(as.integer(sites))
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
            data.frame(
                "Accession" = tips,
                "Pos" = Pos,
                "SNP" = SNP
            )
        }))
        attr(res, "allSNP") <- allSNP
    }
    attr(res, "paths") <- attr(x, "paths")
    class(res) <- "paraFixSites"
    return(res)
}
