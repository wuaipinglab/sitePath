#' @rdname fixationPath
#' @name fixationPath
#' @title Accumulation of fixed mutation as a tree
#' @description The tips are clustered according to the fixation sites. The
#'   transition of fixation sites will be plotted as a phylogenetic tree. The
#'   length of each branch represents the number of fixation mutation between
#'   two clusters. The name of the tree tips indicate the number of sequences in
#'   the cluster.
#' @param x The return from \code{\link{fixationSites}} function.
#' @param minEffectiveSize The minimum size for a tip cluster.
#' @param ... Further arguments passed to or from other methods.
#' @return An \code{fixationPath} object
#' @importFrom stats na.omit
#' @importFrom tidytree as_tibble full_join as.treedata
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' paths <- lineagePath(tree)
#' mutations <- fixationSites(paths)
#' fixationPath(mutations)
fixationPath.fixationSites <- function(x,
                                       minEffectiveSize = NULL,
                                       ...) {
    grouping <- attr(x, "clustersByPath")
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <-
            mean(lengths(unlist(grouping, recursive = FALSE)))
    }
    # Filter out small sized cluster except the divergent points
    grouping <- lapply(grouping, function(g) {
        g[which(vapply(
            X = g,
            FUN = function(tips) {
                length(tips) > minEffectiveSize || !is.null(attr(tips, "toMerge"))
            },
            FUN.VALUE = logical(1)
        ))]
    })
    # The existing tips in the tree
    tipClusters <- list()
    # The edges and its SNP of the tree
    parentNodes <- integer()
    childrenNodes <- integer()
    edgeSNPs <- list()

    # Keep track of the newly added internal node
    newParentNode <-
        length(unlist(grouping, recursive = FALSE)) + 1L
    parentNode <- newParentNode
    # Keep track of the newly added tip node
    tipNode <- 0L
    # A list to record the fixation sites of the parent nodes
    parentNodesSites <- list()

    for (gpIndex in seq_along(grouping)) {
        # The group to add onto the tree
        gp <- grouping[[gpIndex]]
        # Skip the path if none of its group is qualified
        if (length(gp) == 0) {
            next
        }
        # The initial tips of the group
        currentTips <- gp[[1]]
        # Assume the initial reference site
        refSites <- attr(currentTips, "site")
        # Find where to merge and parent node, update reference site maybe
        for (i in seq_along(tipClusters)) {
            toMerge <- attr(tipClusters[[i]], "toMerge")
            if (!is.null(toMerge) &&
                gpIndex %in% as.integer(names(toMerge))) {
                refSites <- toMerge[[as.character(gpIndex)]]
                parentNode <-
                    parentNodes[which(childrenNodes == i)]
                break
            }
        }
        parentNodesSites[[as.character(parentNode)]] <- refSites

        # Track the tip and internal node
        tipNode <- tipNode + 1L
        startingNode <- tipNode
        tipClusters <- c(tipClusters, list(currentTips))
        # Define the initial edge of the group
        parentNodes <- c(parentNodes, parentNode)
        childrenNodes <- c(childrenNodes, tipNode)
        # SNP of the initial edge is set none
        edgeSNPs <- c(edgeSNPs, list(as.character(na.omit(
            vapply(
                X = names(refSites),
                FUN = function(site) {
                    ref <- refSites[site]
                    snp <- attr(currentTips, "site")[site]
                    if (ref == snp) {
                        return(NA_character_)
                    }
                    return(paste0(ref, site, snp))
                },
                FUN.VALUE = character(1)
            )
        ))))

        # Grow the tree
        for (tipIndex in seq_along(gp)[-1]) {
            tipNode <- tipNode + 1L
            currentTips <- gp[[tipIndex]]
            currentSites <- attr(currentTips, "site")
            # Attach the tip near the most related tips. Assume the reference
            # tips are the most related (least number of SNP)
            mostRelatedTipNode <- startingNode
            leastSNPnum <- sum(refSites != currentSites)
            # Loop through the rest existing tip clusters
            for (otherTipNode in seq_along(tipClusters)[-seq_len(startingNode)]) {
                otherSites <- attr(tipClusters[[otherTipNode]], "site")
                snpNum <- sum(otherSites != currentSites)
                if (snpNum < leastSNPnum) {
                    mostRelatedTipNode <- otherTipNode
                    leastSNPnum <- snpNum
                }
            }
            # Find the direct tree edge to the most related tips
            edgeIndex <-
                which(childrenNodes == mostRelatedTipNode)
            parentNode <- parentNodes[edgeIndex]
            # Tree growing differs according to the edge SNP
            parentSites <-
                parentNodesSites[[as.character(parentNode)]]
            snpSites <- as.character(na.omit(
                vapply(
                    X = names(parentSites),
                    FUN = function(site) {
                        ref <- parentSites[site]
                        snp <- currentSites[site]
                        if (ref == snp) {
                            return(NA_character_)
                        }
                        return(paste0(ref, site, snp))
                    },
                    FUN.VALUE = character(1)
                )
            ))
            edgeSNP <- edgeSNPs[[edgeIndex]]
            sharedWithEdgeSNP <- intersect(snpSites, edgeSNP)
            # A new internal node is needed when no SNP overlap
            if (length(sharedWithEdgeSNP) != 0) {
                newParentNode <- newParentNode + 1L
                # Insert the new internal node to the target edge
                parentNodes[edgeIndex] <- newParentNode
                parentNodes <- c(parentNodes, parentNode)
                childrenNodes <- c(childrenNodes, newParentNode)
                edgeSNPs <- c(edgeSNPs, list(sharedWithEdgeSNP))
                # Update SNP of the directly linked edge
                edgeSNPs[[edgeIndex]] <-
                    setdiff(edgeSNP, sharedWithEdgeSNP)
                # Calculate the site for the new internal node
                siteToChange <- substr(sharedWithEdgeSNP,
                                       2,
                                       nchar(sharedWithEdgeSNP) - 1)
                parentSites[siteToChange] <- substr(
                    sharedWithEdgeSNP,
                    nchar(sharedWithEdgeSNP),
                    nchar(sharedWithEdgeSNP)
                )
                parentNodesSites[[as.character(newParentNode)]] <-
                    parentSites
                # Update the parent node and edge SNP for the current tip node
                parentNode <- newParentNode
                snpSites <- setdiff(snpSites, sharedWithEdgeSNP)
            }
            # Add edge
            parentNodes <- c(parentNodes, parentNode)
            childrenNodes <- c(childrenNodes, tipNode)
            # Add edge SNP
            edgeSNPs <- c(edgeSNPs, list(snpSites))
            # Add the current tips
            tipClusters <- c(tipClusters, list(currentTips))
        }
    }
    names(tipClusters) <- vapply(
        X = tipClusters,
        FUN = attr,
        which = "clsName",
        FUN.VALUE = character(1)
    )
    SNPtracing <- list(
        "edge" = cbind(parentNodes, childrenNodes),
        "edge.length" = lengths(edgeSNPs),
        "Nnode" = length(unique(parentNodes)),
        "tip.label" = names(tipClusters)
    )
    class(SNPtracing) <- "phylo"
    names(edgeSNPs) <- SNPtracing[["edge"]][seq_along(edgeSNPs), 2]
    attr(tipClusters, "SNPtracing") <-
        .annotateSNPonTree(SNPtracing, edgeSNPs)
    class(tipClusters) <- "fixationPath"
    return(tipClusters)
}

.annotateSNPonTree <- function(tree, branchSNPs) {
    d <- as_tibble(t(vapply(
        X = names(branchSNPs),
        FUN = function(n) {
            snp <- branchSNPs[[n]]
            if (length(snp) == 0) {
                res <- NA_character_
            } else {
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
            }
            res <- c(res, n)
            names(res) <- c("SNPs", "node")
            return(res)
        },
        FUN.VALUE = c(character(1), integer(1))
    )))
    d[["node"]] <- as.integer(d[["node"]])
    tree <- as_tibble(tree)
    tree <- full_join(tree, d, by = "node")
    return(as.treedata(tree))
}

#' @export
fixationPath <- function(x, ...) {
    UseMethod("fixationPath")
}

#' @export
print.fixationPath <- function(x, ...) {
    print(names(x))
}

#' @rdname plotFunctions
#' @importFrom ggtree geom_tiplab theme_tree2
#' @importFrom ggrepel geom_label_repel
#' @export
plot.fixationPath <- function(x,
                              y = TRUE,
                              ...) {
    tr <- attr(x, "SNPtracing")
    p <- ggtree(tr) +
        geom_tiplab(hjust = 0.5,
                    align = TRUE,
                    offset = 0.5) +
        theme_tree2()
    if (y) {
        p <- p + geom_label_repel(
            aes(x = branch, label = SNPs),
            fill = 'lightgreen',
            min.segment.length = 0,
            na.rm = TRUE
        )
    }
    return(p)
}
