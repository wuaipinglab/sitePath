#' @export
as.phylo.phyMSAmatched <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.sitePath <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.sitesMinEntropy <- function(x, ...) {
    paths <- attr(x, "paths")
    res <- attr(paths, "tree")
    return(res)
}

#' @export
as.phylo.fixationSites <- function(x, ...) {
    paths <- attr(x, "paths")
    res <- attr(paths, "tree")
    return(res)
}

#' @export
as.treedata.fixationSites <- function(tree, ...) {
    x <- tree
    tree <- as.phylo.fixationSites(x)
    grp <- groupTips.fixationSites(tree = x, tipnames = FALSE)
    mutTable <- .mutationTable(x, tree, grp)
    transMut <- lapply(X = split(mutTable, mutTable[, "node"]),
                       FUN = "[[",
                       i = "mutation")
    tree <- .annotateSNPonTree(tree, transMut)
    tree <- groupOTU(tree, grp, group_name = "Groups")
    return(tree)
}

#' @export
as.treedata.fixationPath <- function(tree, ...) {
    res <- attr(tree, "SNPtracing")
    return(res)
}

#' @rdname as.data.frame
#' @title Convert results to Data Frame
#' @description Convert return of functions in \code{sitePath} package to a
#'   \code{\link{data.frame}} so can be better worked with.
#' @param x A \code{\link{fixationSites}} object.
#' @param row.names NULL or a character vector giving the row names for the data
#'   frame. Missing values are not allowed.
#' @param optional Unimplemented.
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
                                        ...) {
    res <- .mutationTable(
        fixations = x,
        tree = as.phylo.fixationSites(x),
        grp = groupTips.fixationSites(x, tipnames = FALSE)
    )
    res <- res[, c("mutation", "from", "to")]
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
            nodeNames <- names(mp)
            for (i in seq_along(mp)[-1]) {
                prevTips <- mp[[i - 1]]
                currTips <- mp[[i]]
                mutation <- paste0(attr(prevTips, "AA"),
                                   site,
                                   attr(currTips, "AA"))
                trans <- nodeNames[i]
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

#' @export
as.data.frame.parallelSites <- function(x,
                                        row.names = NULL,
                                        optional = FALSE,
                                        ...) {
    tree <- as.phylo.phyMSAmatched(attr(x, "paths"))
    tipNames <- tree[["tip.label"]]
    clustersByPath <- attr(x, "clustersByPath")
    clusterInfo <- character()
    for (gp in clustersByPath) {
        for (tips in gp) {
            # The cluster named by the tips
            clsNames <- rep(attr(tips, "clsName"), length(tips))
            names(clsNames) <- as.character(tipNames[tips])
            clusterInfo <- c(clusterInfo, clsNames)
        }
    }
    # Info for the parallel mutation
    accession <- character()
    clsName <- character()
    mutSite <- integer()
    mutFrom <- character()
    mutTo <- character()
    isFixed <- logical()
    for (sp in x) {
        tips <- extractTips.sitePara(sp)
        for (t in tips) {
            accession <- c(accession, t)
            clsName <- c(clsName, clusterInfo[t])
            mutName <- attr(t, "mutName")
            mutSite <- c(mutSite, rep(mutName[2], length(t)))
            mutFrom <- c(mutFrom, rep(mutName[1], length(t)))
            mutTo <- c(mutTo, rep(mutName[3], length(t)))
            isFixed <- c(isFixed, rep(attr(t, "fixed"), length(t)))
        }
    }
    res <- data.frame(
        "Accession" = accession,
        "group" = clsName,
        "site" = mutSite,
        "mutFrom" = mutFrom,
        "mutTo" = mutTo,
        "fixation" = isFixed
    )
    return(res)
}
