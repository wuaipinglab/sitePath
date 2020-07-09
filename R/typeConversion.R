#' @rdname as.data.frame
#' @title Convert results to Data Frame
#' @description Convert return of functions in \code{sitePath} package to a
#'   \code{\link{data.frame}} so can be better worked with.
#' @param x A \code{\link{fixationSites}} object.
#' @param row.names NULL or a character vector giving the row names for the data
#'   frame. Missing values are not allowed.
#' @param optional Unimplemented.
#' @param tipname Logical: if the return give fixation mutations between tip
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
    grp <- fixationPath.fixationSites(x = x, minEffectiveSize = 0)
    grp <- as.list.fixationPath(grp)
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
        res <- .mutationTable(x, tree, grp)
        res <- res[, c("mutation", "from", "to")]
    }
    return(res)
}

.mutationTable <- function(fixations, tree, grp) {
    # The cluster each tip belongs to
    clusterInfo <- character()
    # Extract the node path for each tip cluster
    clusterPaths <- list()
    rootNode <- getMRCA(tree, tree[["tip.label"]])
    # Iterate the cluster and grow the two above
    for (cluster in names(grp)) {
        tips <- grp[[cluster]]
        # The cluster named by the tips
        clsNames <- rep(cluster, length(tips))
        names(clsNames) <- as.character(tips)
        clusterInfo <- c(clusterInfo, clsNames)
        # The node path towards the cluster
        ancestral <- getMRCA(tree, tips)
        if (is.null(ancestral)) {
            np <- nodepath(tree, rootNode, tips)
            clusterPaths[[cluster]] <- np[seq_len(length(np) - 1)]
        } else {
            clusterPaths[[cluster]] <- nodepath(tree, rootNode, ancestral)
        }
    }
    # Info for the transition mutation
    prevCls <- character()
    currCls <- character()
    mutName <- character()
    transNode <- integer()
    for (sp in fixations) {
        site <- attr(sp, "site")
        for (mp in sp) {
            for (i in seq_along(mp)[-1]) {
                prevTips <- mp[[i - 1]]
                currTips <- mp[[i]]
                mutation <- paste0(attr(prevTips, "AA"),
                                   site,
                                   attr(currTips, "AA"))
                prev <- unique(clusterInfo[as.character(prevTips)])
                names(prev) <- prev
                # Choose the most recent cluster to stay un-mutated
                prev <- names(which.max(lapply(
                    X = prev,
                    FUN = function(cluster) {
                        length(clusterPaths[[cluster]])
                    }
                )))
                curr <- unique(clusterInfo[as.character(currTips)])
                names(curr) <- curr
                # Choose the most ancient cluster which first receive the
                # mutation
                curr <- names(which.min(lapply(
                    X = curr,
                    FUN = function(cluster) {
                        length(clusterPaths[[cluster]])
                    }
                )))
                # Find the transition node
                currPath <- clusterPaths[[curr]]
                trans <- currPath[length(currPath)]
                # Add the new transition mutation
                prevCls <- c(prevCls, prev)
                currCls <- c(currCls, curr)
                mutName <- c(mutName, mutation)
                transNode <- c(transNode, trans)
            }
        }
    }
    # The mutation between adjacent clusters
    res <- data.frame(
        "mutation" = mutName,
        "from" = prevCls,
        "to" = currCls,
        "node" = transNode
    )
    return(res)
}

#' @export
as.data.frame.SNPsites <- function(x,
                                   row.names = NULL,
                                   optional = FALSE,
                                   ...) {
    res <- attr(x, "allSNP")
    return(res)
}

#' @importFrom tidytree as.treedata
#' @export
tidytree::as.treedata

#' @export
as.treedata.fixationSites <- function(tree, ...) {
    x <- tree
    tree <- as.phylo.fixationSites(x)
    grp <- fixationPath.fixationSites(x, minEffectiveSize = 0)
    grp <- as.list.fixationPath(grp)
    mutTable <- .mutationTable(x, tree, grp)
    transMut <- lapply(X = split(mutTable, mutTable[, "node"]),
                       FUN = "[[",
                       i = "mutation")
    d <- as_tibble(t(vapply(
        X = names(transMut),
        FUN = function(trans) {
            snp <- transMut[[trans]]
            res <- character()
            snpNum <- length(snp)
            for (i in seq_len(snpNum)) {
                res <- paste0(res, snp[i])
                if (i < snpNum) {
                    if (i %% 4 == 0) {
                        res <- paste0(res, ",\n")
                    } else {
                        res <- paste0(res, ", ")
                    }
                }
            }
            res <- c(trans, res)
            names(res) <- c("node", "SNPs")
            return(res)
        },
        FUN.VALUE = character(2)
    )))
    d[["node"]] <- as.integer(d[["node"]])
    d <- full_join(as_tibble(tree), d, by = "node")
    tree <- as.treedata(d)
    tree <- groupOTU(tree, grp, group_name = "Groups")
    return(tree)
}

#' @export
as.treedata.fixationPath <- function(tree, ...) {
    res <- attr(tree, "SNPtracing")
    return(res)
}

#' @importFrom ape as.phylo
#' @export
ape::as.phylo

#' @export
as.phylo.phyMSAmatched <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.lineagePath <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.sitePath <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.fixationSites <- function(x, ...) {
    paths <- attr(x, "paths")
    res <- attr(paths, "tree")
    return(res)
}

#' @export
as.list.fixationPath <- function(x, ...) {
    groupName <- names(x)
    attributes(x) <- NULL
    res <- lapply(x, function(tips) {
        attributes(tips) <- NULL
        return(tips)
    })
    names(res) <- groupName
    return(res)
}
