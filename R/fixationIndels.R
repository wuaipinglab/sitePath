#' @rdname fixationIndels
#' @title Fixation indels prediction
#' @description The fixation of insertions of deletions.
#' @param x The return from \code{\link{sitesMinEntropy}} function.
#' @param ... Other arguments.
#' @return A \code{fixationIndels} object.
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' fixationIndels(sitesMinEntropy(tree))
fixationIndels <- function(x, ...) {
    UseMethod("fixationIndels")
}

#' @rdname fixationIndels
#' @export
fixationIndels.sitesMinEntropy <- function(x, ...) {
    paths <- attr(x, "paths")
    tree <- attr(paths, "tree")
    seqType <- attr(paths, "seqType")
    gapChar <- attr(paths, "gapChar")
    minSize <- attr(paths, "minSize")
    # 'res' is going to be the return of this function. Each entry in the list
    # is the 'indelPath' for a fragment of sequence.
    res <- list()
    for (segs in x) {
        pathNodeTips <- lapply(attr(segs, "pathNodeTips"), as.integer)
        prevSite <- -1
        currIndels <- list()
        for (site in names(segs)) {
            seg <- segs[[site]]
            # Find the tips having 'gapChar' at the site
            siteChars <- vapply(
                X = seg,
                FUN = attr,
                FUN.VALUE = character(1),
                which = "AA"
            )
            tipsWithDeletion <- seg[which(siteChars == gapChar)]
            if (length(tipsWithDeletion)) {
                currSite <- as.integer(site)
                # Test the continuity of the deletion
                if (currSite - prevSite == 1) {
                    # Find the overlapping tips to further ensure the continuity
                    for (iter in seq_along(currIndels)) {
                        # Existing tips with continuing deletion
                        refTips <- currIndels[[iter]]
                        indelSites <- c(attr(refTips, "indelSites"),
                                        currSite)
                        for (tips in tipsWithDeletion) {
                            continued <- intersect(refTips, tips)
                            # The deletion of the tips is ended if the current
                            # site is not gap
                            ended <- setdiff(refTips, continued)
                            # A new deletion is started if a new group of tips
                            # are gap at the current site
                            started <- setdiff(tips, continued)
                            if (length(continued)) {
                                continued <- .findAncestralNode(continued,
                                                                pathNodeTips,
                                                                indelSites)
                                currIndels[iter] <- continued
                                if (length(ended)) {
                                    ended <- .findAncestralNode(ended,
                                                                pathNodeTips,
                                                                indelSites)
                                    currIndels <-
                                        c(currIndels, ended)
                                }
                            } else {
                                if (length(ended)) {
                                    ended <- .findAncestralNode(ended,
                                                                pathNodeTips,
                                                                indelSites)
                                    currIndels[iter] <- ended
                                }
                            }
                            if (length(started)) {
                                started <- .findAncestralNode(started,
                                                              pathNodeTips,
                                                              currSite)
                                currIndels <- c(currIndels, started)
                            }
                        }
                    }
                } else {
                    # Initiate the first deletion fragment or re-initiate new
                    # deletion fragment if the gap can't be extended due to
                    # discontunity of the site
                    currIndels <-
                        lapply(tipsWithDeletion, function(tips) {
                            attr(tips, "indelSites") <- currSite
                            return(tips)
                        })
                }
                # Update the 'prevSite' only when the site is a gap
                prevSite <- currSite
            }
        }
        # All indel for the current path
        for (tips in currIndels) {
            if (length(tips) >= minSize) {
                indelSites <- attr(tips, "indelSites")
                if (length(indelSites) > 1) {
                    indelSites <- range(indelSites)
                    indelSites <- paste0(indelSites, collapse = "-")
                }
                node <- attr(tips, "node")
                res[[indelSites]][[node]] <- tips
                attr(res[[indelSites]], "indelSites") <- indelSites
                attr(res[[indelSites]], "tree") <- tree
                attr(res[[indelSites]], "seqType") <- seqType
                class(res[[indelSites]]) <- "indelPath"
            }
        }
    }
    # Set 'paths' and 'clustersByPath' attributes
    attr(res, "paths") <- paths
    class(res) <- "fixationIndels"
    return(res)
}

.findAncestralNode <- function(tipsWithGap,
                               pathNodeTips,
                               indelSites) {
    res <- list()
    # The tips to be grouped
    currTips <- integer()
    ancestralNode <- names(pathNodeTips[1])
    # Iterate the tips along the path
    for (node in names(pathNodeTips)) {
        tips <- pathNodeTips[[node]]
        # To find the tips that are grouped in the 'tipsWithGap'
        if (any(tips %in% tipsWithGap)) {
            # The ancestral node is from its starting tips in 'pathNodeTips'
            if (length(currTips) == 0) {
                ancestralNode <- node
            }
            # Accumulating the tips (it's assumed the all 'tips' are in
            # 'tipsWithGap' if any)
            currTips <- c(currTips, tips)
            tipsWithGap <- setdiff(tipsWithGap, tips)
        } else {
            if (length(currTips)) {
                # The continuity stopped and new tip group formed
                attr(currTips, "node") <- ancestralNode
                res[[ancestralNode]] <- currTips
            }
            # Reset the tips to be grouped
            currTips <- integer()
        }
    }
    if (length(currTips)) {
        # The continuity stopped and new tip group formed
        attr(currTips, "node") <- ancestralNode
        attr(currTips, "indelSites") <- indelSites
        res[[ancestralNode]] <- currTips
    }
    return(res)
}
