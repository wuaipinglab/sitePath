#' @rdname as.data.frame
#' @title Convert results to Data Frame
#' @description Convert return of functions in \code{sitePath} package to a
#'   \code{\link{data.frame}} so can be better worked with.
#' @param x A \code{\link{fixationSites}} object.
#' @param row.names NULL or a character vector giving the row names for the data
#'   frame. Missing values are not allowed.
#' @param optional Unimplemented.
#' @param tipname Logical: if the return give fixation mutations betweem tip
#'   groups or tip names within a group. The default is \code{FALSE}.
#' @param ... Other arguments.
#' @return A \code{\link{data.frame}} object.
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' fixations <- fixationSites(lineagePath(tree))
#' as.data.frame(fixations)
as.data.frame.fixationSites <- function(x,
                                        row.names = NULL,
                                        optional = FALSE,
                                        tipname = FALSE,
                                        ...) {
    tree <- as.phylo.fixationSites(x)
    grp <-
        sitewiseClusters.fixationSites(x, minEffectiveSize = 0)
    grp <- as.list.sitewiseClusters(grp)
    if (tipname) {
        res <- lapply(names(grp), function(n) {
            tips <- grp[[n]]
            data.frame(
                row.names = as.character(tips),
                "tipname" = tree[["tip.label"]][tips],
                "group" = rep(n, length(tips))
            )
        })
        res <- do.call(rbind, res)

    } else {
        clusterPaths <- list()
        rootNode <- getMRCA(tree, tree[["tip.label"]])

        for (cluster in names(grp)) {
            tips <- grp[[cluster]]
            ancestral <- getMRCA(tree, tips)
            if (is.null(ancestral)) {
                np <- nodepath(tree, rootNode, tips)
                clusterPaths[[cluster]] <-
                    np[1:(length(np) - 1)]
            } else {
                clusterPaths[[cluster]] <- nodepath(tree, rootNode, ancestral)
            }
        }
        clusterInfo <- lapply(names(grp), function(g) {
            data.frame(row.names = grp[[g]],
                       "cluster" = rep(g, length(grp[[g]])))
        })
        clusterInfo <- do.call(rbind, clusterInfo)
        transMut <- list()
        for (sp in x) {
            site <- attr(sp, "site")
            for (mp in sp) {
                for (i in seq_along(mp)[-1]) {
                    prevTips <- mp[[i - 1]]
                    currTips <- mp[[i]]
                    prevAA <- attr(prevTips, "AA")
                    currAA <- attr(currTips, "AA")
                    mutation <- paste0(prevAA, site, currAA)
                    prevCluster <-
                        unique(clusterInfo[as.character(prevTips), ])
                    names(prevCluster) <- prevCluster
                    # Choose the most recent cluster to stay un-mutated
                    prev <- names(which.max(lapply(
                        X = prevCluster,
                        FUN = function(cluster) {
                            length(clusterPaths[[cluster]])
                        }
                    )))
                    currCluster <-
                        unique(clusterInfo[as.character(currTips), ])
                    names(currCluster) <- currCluster
                    # Choose the most ancient cluster which first receive the
                    # mutation
                    curr <- names(which.min(lapply(
                        X = currCluster,
                        FUN = function(cluster) {
                            length(clusterPaths[[cluster]])
                        }
                    )))
                    trans <- paste(prev, curr, sep = "-")
                    if (trans %in% names(transMut)) {
                        transMut[[trans]] <- c(transMut[[trans]], mutation)
                    } else {
                        transMut[[trans]] <- mutation
                    }
                }
            }
        }
        transMut <- lapply(transMut, unique)
        res <- lapply(strsplit(names(transMut), '-'), function(i) {
            n <- paste(i, collapse = '-')
            data.frame(
                "mutation" = transMut[[n]],
                "from" = rep(i[1], length(transMut[[n]])),
                "to" = rep(i[2], length(transMut[[n]]))
            )
        })
        res <- do.call(rbind, res)
    }
    return(res)
}

as.data.frame.multiFixationSites <- function(x, ...) {
    cat("Under development\n")
}

#' @importFrom ape as.phylo
#' @export
ape::as.phylo

#' @export
as.phylo.lineagePath <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.fixationSites <- function(x, ...) {
    paths <- attr(x, "paths")
    res <- as.phylo.lineagePath(paths)
    return(res)
}

#' @export
as.list.sitewiseClusters <- function(x, ...) {
    groupName <- names(x)
    attributes(x) <- NULL
    res <- lapply(x, function(tips) {
        attributes(tips) <- NULL
        return(tips)
    })
    names(res) <- groupName
    return(res)
}
