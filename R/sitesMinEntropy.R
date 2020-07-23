#' @export
sitesMinEntropy.lineagePath <- function(x,
                                        minEffectiveSize = NULL,
                                        searchDepth = 1,
                                        method = c("compare", "insert", "delete"),
                                        ...) {
    paths <- .phyMSAmatch(x)
    tree <- attr(paths, "tree")
    # Set the minimal size of the group during the search
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <-
            length(tree[["tip.label"]]) / length(unique(unlist(paths)))
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    minEffectiveSize <- ceiling(minEffectiveSize)
    # Set the search depth for heuristic search
    if (searchDepth < 1) {
        stop("\"searchDepth\" should be at least 1")
    } else {
        searchDepth <- ceiling(searchDepth)
    }
    # Decide which minimizing strategy
    minimizeEntropy <- switch(
        match.arg(method),
        "compare" = minEntropyByComparing,
        "insert" = minEntropyByInserting,
        "delete" = minEntropyByDeleting
    )
    # Get the divergent nodes
    divNodes <- divergentNode(paths)
    # The tips and matching
    nodeAlign <- .tipSeqsAlongPathNodes(paths = paths,
                                        divNodes = divNodes)
    # Get the MSA numbering
    reference <- attr(paths, "msaNumbering")
    align <- attr(paths, "align")
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    # In case root node does not have any tips (because itself is a divergent
    # node)
    if (!attr(paths, "rootNode") %in% names(nodeAlign)) {
        paths <- lapply(paths, function(p) {
            as.character(setdiff(p[-1], divNodes))
        })
    } else {
        paths <- lapply(paths, function(p) {
            as.character(setdiff(p, divNodes))
        })
    }
    # Turn the site number into index for C++ code
    siteIndices <- reference[loci] - 1
    names(siteIndices) <- as.character(loci)
    # Group the result by path for all loci
    res <- lapply(paths, function(path) {
        # Entropy minimization result for every locus
        lapply(
            X = siteIndices,
            FUN = .runEntropyMinimization,
            path = path,
            nodeAlign = nodeAlign,
            minimizeEntropy = minimizeEntropy,
            minEffectiveSize = minEffectiveSize,
            searchDepth = searchDepth
        )
    })
    # Cluster tips according to fixation sites
    clustersByPath <- .clustersByPath(res)
    clustersByPath <- .mergeClusters(clustersByPath)
    attr(res, "clustersByPath") <-
        .assignClusterNames(clustersByPath)
    attr(res, "paths") <- paths
    class(res) <- "sitesMinEntropy"
    return(res)
}

.runEntropyMinimization <- function(siteIndex,
                                    path,
                                    nodeAlign,
                                    minimizeEntropy,
                                    minEffectiveSize,
                                    searchDepth) {
    # Assign a variable to store the tip names and their info on amino
    # acids. They are the potential fixation segment
    nodeTips <- integer()
    previousAA <- NULL
    currentAA <- NULL
    previousNode <- NULL
    # The input for entropy minimization calculation
    nodeSummaries <- list()
    # Divergent nodes are not included anywhere in the result
    for (node in path) {
        # Get the related descendant tips and related sequences
        nodeTips <- as.integer(names(nodeAlign[[node]]))
        # Frequency of the amino acids at the locus
        aaSummary <- tableAA(nodeAlign[[node]], siteIndex)
        # Associate the amino acid frequency with the tip names
        attr(nodeTips, "aaSummary") <- aaSummary
        # Decide the current fixed amino acid
        if (length(aaSummary) == 1) {
            currentAA <- names(aaSummary)
        } else {
            currentAA <- NULL
        }
        # Attach the node to the previous node if they're both purely fixed
        # and have the same AA fixed.
        if (!is.null(previousAA) &&
            !is.null(currentAA) &&
            previousAA == currentAA) {
            node <- previousNode
            # Combine the tips in the previous node
            nodeTips <- c(nodeSummaries[[node]], nodeTips)
            # Add up the amino acid frequency
            attr(nodeTips, "aaSummary") <-
                attr(nodeSummaries[[node]], "aaSummary") +
                aaSummary
            # Oddly, R uses the name of the first variable when adding two
            # numeric vectors. So there is no need for names (AA) assignment
        }
        # Assign or re-assign the 'nodeTips' with 'aaSummary' to the
        # 'nodeSummaries'
        nodeSummaries[[node]] <- nodeTips
        previousAA <- currentAA
        previousNode <- node
    }
    # Return empty value if the site is purely fixed on the lineage
    seg <- list()
    if (length(nodeSummaries) >= 2) {
        seg <- minimizeEntropy(nodeSummaries,
                               minEffectiveSize,
                               searchDepth)
    }
    return(seg)
}

.clustersByPath <- function(fixations) {
    # Find the clustering for each lineage path
    res <- lapply(fixations, function(sp) {
        # Remove the site purely conserved on the lineage
        group <- list()
        for (site in names(sp)) {
            mp <- sp[[site]]
            if (length(mp) >= 2) {
                toAdd <- lapply(mp, function(tips) {
                    siteChar <- attr(tips, "AA")
                    attributes(tips) <- NULL
                    attr(tips, "site") <- siteChar
                    names(attr(tips, "site")) <- site
                    tips
                })
                group <- c(group, list(toAdd))
            }
        }
        if (length(group) == 0) {
            return(group)
        }
        # Group tips according to fixation points
        res <- group[[1]]
        for (p in group[-1]) {
            for (tips in p) {
                site <- attr(tips, "site")
                # Update grouping for each tips by growing a new list
                newGrouping <- list()
                for (i in seq_along(res)) {
                    gp <- res[[i]]
                    common <- sort(intersect(tips, gp))
                    # No new cluster when the coming tips have no overlap or are
                    # identical to tips in an existing cluster
                    if (length(common) == 0) {
                        newGrouping <- res[seq_len(i)]
                    } else if (identical(sort(gp), sort(tips))) {
                        attr(gp, "site") <- c(attr(gp, "site"), site)
                        if (i + 1 <= length(res)) {
                            trailing <- res[(i + 1):length(res)]
                        } else {
                            trailing <- list()
                        }
                        newGrouping <-
                            c(newGrouping, list(gp), trailing)
                        break
                    } else {
                        # A new cluster formed when there is overlapped between
                        # new coming tips and existing tips in a cluster
                        if (identical(sort(gp), common)) {
                            # The new coming tips includes the current group
                            # The extra tips stay for the next loop
                            tips <- setdiff(tips, gp)
                            # Update the SNP site info for the current group
                            attr(gp, "site") <-
                                c(attr(gp, "site"), site)
                            newGrouping <- c(newGrouping, list(gp))
                        } else if (identical(sort(tips), common)) {
                            # The new coming tips are included in the group
                            # (they are used up at this point)
                            separate <- setdiff(gp, tips)
                            attributes(separate) <- attributes(gp)
                            attr(tips, "site") <-
                                c(attr(gp, "site"), site)
                            if (i + 1 <= length(res)) {
                                trailing <- res[(i + 1):length(res)]
                            } else {
                                trailing <- list()
                            }
                            newGrouping <- c(newGrouping,
                                             list(tips),
                                             list(separate),
                                             trailing)
                            # Go for the next new coming tips
                            break
                        } else {
                            stop("Something's not right")
                        }
                    }
                }
                # The new coming tips are used up and update the grouping
                res <- newGrouping
            }
        }
        return(res)
    })
    res <- res[which(lengths(res) != 0)]
    return(res)
}

#' @export
sitesMinEntropy <- function(x, ...) {
    UseMethod("sitesMinEntropy")
}

#' @export
print.sitesMinEntropy <- function(x, ...) {
    cat("This is a 'sitesMinEntropy' object.", "\n")
}
