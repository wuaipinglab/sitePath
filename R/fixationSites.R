#' @rdname fixationSites
#' @title Fixation sites prediction
#' @description After finding the \code{\link{lineagePath}} of a phylogenetic
#'   tree, \code{fixationSites} uses the result to find those sites that show
#'   fixation on some, if not all, of the lineages. The number of tips before
#'   and after the fixation mutation is expected to be more than
#'   \code{minEffectiveSize}. Also, the fixation will be skipped if the amino
#'   acid/nucleotide is gap or ambiguous character. A lineage has to have at
#'   least one fixation mutation to be reported.
#' @param paths A \code{lineagePath} object returned from
#'   \code{\link{lineagePath}} function.
#' @param minEffectiveSize The minimum number of tips in a group.
#' @param searchDepth The function uses heuristic search but the termination of
#'   the search cannot be intrinsically decided. \code{searchDepth} is needed to
#'   tell the search when to stop.
#' @param method The strategy for predicting the fixation. The basic approach is
#'   entropy minimization and can be achieved by adding or removing fixation
#'   point, or by comparing the two.
#' @param ... further arguments passed to or from other methods.
#' @return A \code{fixationSites} object.
#' @seealso \code{\link{as.data.frame.fixationSites}}
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' fixationSites(lineagePath(tree))
fixationSites <- function(paths, ...) {
    UseMethod("fixationSites")
}

#' @rdname fixationSites
#' @export
fixationSites.lineagePath <- function(paths,
                                      minEffectiveSize = NULL,
                                      searchDepth = 1,
                                      method = c("compare", "insert", "delete"),
                                      ...) {
    minEntropy <- sitesMinEntropy.lineagePath(paths,
                                              minEffectiveSize,
                                              searchDepth,
                                              method)
    res <- fixationSites.sitesMinEntropy(minEntropy, ...)
    return(res)
}

#' @rdname fixationSites
#' @export
fixationSites.sitesMinEntropy <- function(paths, ...) {
    x <- paths
    paths <- attr(x, "paths")
    tree <- attr(paths, "tree")
    align <- attr(paths, "align")
    seqType <- attr(paths, "seqType")
    unambiguous <- .unambiguousChars(paths)
    # 'res' is going to be the return of this function. Each entry in the list
    # is the 'sitePath' for a site. Each site ('sitePath') consists of 'mutPath'
    # that is named by the starting node name. The fixed AA and number of
    # non-dominant AA is also stored.
    res <- list()
    for (segs in x) {
        for (site in names(segs)) {
            seg <- segs[[site]]
            # There has to be at least one fixation on the lineage and at least
            # two of the mutation is neither gap nor ambiguous character
            if (.qualifiedFixation(seg, unambiguous)) {
                i <- as.integer(site)
                # Test if the slot for the site is empty
                if (is.null(res[[site]])) {
                    # Initiate the first 'mutPath' for the site
                    res[[site]][[1]] <- lapply(seg, function(tips) {
                        attr(tips, "tipsAA") <- substr(
                            x = align[tips],
                            start = i,
                            stop = i
                        )
                        return(tips)
                    })
                    attr(res[[site]], "site") <- i
                    attr(res[[site]], "tree") <- tree
                    attr(res[[site]], "seqType") <- seqType
                    class(res[[site]]) <- "sitePath"
                } else {
                    # Assume a new 'mutPath' is to add (not combined by default)
                    targetIndex <- length(res[[site]]) + 1
                    # The index to extract the terminal tips of the 'mutPath'
                    endIndex <- length(seg)
                    finalAA <- attr(seg[[endIndex]], "AA")
                    # The following is to decide if any 'mutPath' can be
                    # combined
                    existPaths <- res[[site]]
                    # Combine the two 'mutPath' when the fixation tips before
                    # the terminal tips are identical and the final fixation
                    # mutation are the same
                    toCombine <- vapply(
                        X = existPaths,
                        FUN = function(ep) {
                            identical(lapply(seg, c)[-endIndex],
                                      lapply(ep, c)[-endIndex]) &&
                                finalAA == attr(ep[[endIndex]], "AA")
                        },
                        FUN.VALUE = logical(1)
                    )
                    if (any(toCombine)) {
                        existIndex <- which(toCombine)
                        # These are the candidates to combine. The additional
                        # condition be all the descendant tips are included.
                        toCombine <- unique(unlist(lapply(
                            X = c(res[[site]][existIndex], list(seg)),
                            FUN = "[[",
                            i = endIndex
                        )))
                        # The all descendant tips of the node where 'mutPath's
                        # are combined
                        allTips <-
                            .childrenTips(tree, getMRCA(tree, toCombine))
                        if (all(allTips %in% toCombine)) {
                            # Create the newly combined 'mutPath'
                            seg[[endIndex]] <- toCombine
                            attr(seg[[endIndex]], "AA") <- finalAA
                            # Remove the combined 'mutPath'
                            res[[site]] <- res[[site]][-existIndex]
                            targetIndex <- length(res[[site]]) + 1
                        } else {
                            # Add new fixation for the site if no existing
                            # mutation path can be combined with
                            targetIndex <- length(existPaths) + 1
                        }
                    }
                    seg <- lapply(seg, function(tips) {
                        attr(tips, "tipsAA") <- substr(
                            x = align[tips],
                            start = i,
                            stop = i
                        )
                        return(tips)
                    })
                    res[[site]][[targetIndex]] <- seg
                    attr(res[[site]], "site") <- i
                    attr(res[[site]], "tree") <- tree
                    attr(res[[site]], "seqType") <- seqType
                    class(res[[site]]) <- "sitePath"
                }
            }
        }
    }
    # Set 'paths' and 'clustersByPath' attributes
    attr(res, "paths") <- paths
    attr(res, "clustersByPath") <- attr(x, "clustersByPath")
    class(res) <- "fixationSites"
    return(res)
}

.qualifiedFixation <- function(seg, unambiguous) {
    siteChars <- unique(vapply(
        X = seg,
        FUN = attr,
        FUN.VALUE = character(1),
        which = "AA"
    ))
    sum(siteChars %in% unambiguous) >= 2
}
