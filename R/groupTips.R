#' @importFrom ape Ntip

#' @rdname groupTips
#' @title The grouping of tree tips
#' @description The tips between divergent nodes or fixation mutations on the
#'   lineages are each gathered as group.
#' @param tree The return from \code{\link{addMSA}}, \code{\link{lineagePath}},
#'   \code{\link{sitesMinEntropy}} or other functions.
#' @param similarity This decides how minor SNPs are to remove. If provided as
#'   fraction between 0 and 1, then the minimum number of SNP will be total tips
#'   times \code{similariy}. If provided as integer greater than 1, the minimum
#'   number will be \code{similariy}. The default \code{similariy} is 0.05 for
#'   \code{lineagePath}.
#' @param simMatrix Deprecated and will not have effect.
#' @param forbidTrivial Does not allow trivial trimming.
#' @param tipnames If return tips as integer or tip names.
#' @param ... Other arguments.
#' @return \code{groupTips} returns grouping of tips.
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' groupTips(tree)
groupTips <- function(tree, ...) {
    UseMethod("groupTips")
}

#' @rdname groupTips
#' @export
groupTips.phyMSAmatched <- function(tree,
                                    similarity = NULL,
                                    simMatrix = NULL,
                                    forbidTrivial = TRUE,
                                    tipnames = TRUE,
                                    ...) {
    paths <- lineagePath.phyMSAmatched(
        tree = tree,
        similarity = simMatrix,
        simMatrix = simMatrix,
        forbidTrivial = forbidTrivial,
        ...
    )
    res <- groupTips.lineagePath(paths, tipnames = tipnames)
    return(res)
}

#' @rdname groupTips
#' @export
groupTips.lineagePath <- function(tree, tipnames = TRUE, ...) {
    paths <- tree
    tree <- attr(paths, "tree")
    # Get the divergent nodes
    divNodes <- divergentNode(paths)
    # The tips and the corresponding ancestral node
    pathNodeTips <- .tipSeqsAlongPathNodes(paths, divNodes)
    # To group the tips by the node right after the divergent point
    res <- list()
    # Iterate through each lineage path
    for (p in paths) {
        # Assume the root node as the first ancestral node
        aNode <- as.character(attr(paths, "rootNode"))
        tips <- integer()
        pathLen <- length(p)
        for (i in seq_len(pathLen)[-pathLen]) {
            currNode <- p[[i]]
            nextNode <- p[[i + 1]]
            # Add the tips of the current node to the group
            tips <- c(tips, pathNodeTips[[as.character(currNode)]])
            # Stop adding the tips to the group and take the group out
            if (nextNode %in% divNodes) {
                res[[aNode]] <- tips
                # The node next to the divergent node is the new ancestral node
                aNode <- as.character(p[[i + 2]])
                tips <- integer()
            }
        }
        # Add the tips of the final node to the group and take the final group
        # out
        res[[aNode]] <- c(tips,
                          pathNodeTips[[as.character(p[[pathLen]])]])
    }
    if (tipnames) {
        res <- lapply(res, function(tips) {
            tree[["tip.label"]][tips]
        })
    }
    return(res)
}

.tipSeqsAlongPathNodes <- function(paths, divNodes) {
    tree <- attr(paths, "tree")
    align <- attr(paths, "align")
    allNodes <- unlist(paths)
    terminalNodes <- vapply(
        X = paths,
        FUN = function(p) {
            p[length(p)]
        },
        FUN.VALUE = integer(1)
    )
    # Get all the nodes that are not at divergent point
    nodes <- setdiff(allNodes, divNodes)
    # Get the sequence of the children tips that are descendant of the nodes.
    # Assign the tip index to the sequences for retrieving the tip name
    nodeAlign <- lapply(nodes, function(n) {
        isTerminal <- FALSE
        if (n %in% terminalNodes) {
            childrenNode <- n
            isTerminal <- TRUE
        } else {
            childrenNode <- tree[["edge"]][which(tree[["edge"]][, 1] == n), 2]
            # Keep the node that is not on the path.
            childrenNode <- setdiff(childrenNode, allNodes)
        }
        res <- .childrenTips(tree, childrenNode)
        attr(res, "align") <- align[res]
        attr(res, "isTerminal") <- isTerminal
        return(res)
    })
    # Assign the node names to the 'nodeAlign' list
    names(nodeAlign) <- nodes
    return(nodeAlign)
}

.childrenTips <- function(tree, node) {
    maxTip <- Ntip(tree)
    children <- integer()
    getChildren <- function(edges, parent) {
        children <<- c(children, parent[which(parent <= maxTip)])
        i <- which(edges[, 1] %in% parent)
        if (length(i) == 0L) {
            return(children)
        } else {
            parent <- edges[i, 2]
            return(getChildren(edges, parent))
        }
    }
    return(getChildren(tree[["edge"]], node))
}

#' @rdname groupTips
#' @export
groupTips.sitesMinEntropy <- function(tree, tipnames = TRUE, ...) {
    .unpackClustersByPath(x = tree, tipnames = tipnames)
}

.unpackClustersByPath <- function(x, tipnames) {
    clustersByPath <- attr(x, "clustersByPath")
    tree <- as.phylo.sitesMinEntropy(x)
    tipLabels <- tree[["tip.label"]]
    if (!tipnames) {
        tipLabels <- seq_along(tipLabels)
    }
    res <- list()
    for (gp in clustersByPath) {
        for (tips in gp) {
            res[[attr(tips, "clsName")]] <- tipLabels[as.integer(tips)]
        }
    }
    return(res)
}

#' @rdname groupTips
#' @export
groupTips.fixationSites <- function(tree, tipnames = TRUE, ...) {
    .unpackClustersByPath(x = tree, tipnames = tipnames)
}

#' @rdname groupTips
#' @export
groupTips.fixationPath <- function(tree, tipnames = TRUE, ...) {
    x <- tree
    groupName <- names(x)
    tree <- attr(x, "tree")
    attributes(x) <- NULL
    if (tipnames) {
        res <- lapply(x, function(tips) {
            attributes(tips) <- NULL
            return(tree[["tip.label"]][tips])
        })
    } else {
        res <- lapply(x, function(tips) {
            attributes(tips) <- NULL
            return(tips)
        })
    }
    names(res) <- groupName
    return(res)
}

#' @rdname similarityMatrix
#' @title Similarity between sequences
#' @description Get similarity between aligned sequences with gap ignored.
#' @param tree The return from \code{\link{addMSA}} function.
#' @return A diagonal matrix of similarity between sequences.
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' simMatrix <- similarityMatrix(tree)
similarityMatrix <- function(tree) {
    sim <- attr(tree, "simMatrix")
    return(sim)
}
