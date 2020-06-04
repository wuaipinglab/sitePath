#' @rdname findSites
#' @name findSites
#' @title Finding sites with variation
#' @description Single nucleotide polymorphism (SNP) in the whole package refers
#'   to variation of amino acid. \code{findSNPsite} will try to find SNP in the
#'   multiple sequence alignment. A reference sequence and gap character may be
#'   specified to number the site. This is irrelevant to the intended analysis
#'   but might be helpful to evaluate the performance of \code{fixationSites}.
#' @param tree The return from \code{\link{addMSA}} function
#' @param minSNP Minimum number of amino acid variation to be a SNP
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' SNPsites(tree)
#' @return \code{SNPsite} returns a list of qualified SNP site
#' @export
SNPsites <- function(tree, minSNP = NULL) {
    if (is.null(minSNP)) {
        minSNP <- length(tree[["tip.label"]]) / 10
    }
    alignedSeq <- attr(tree, "align")
    if (is.null(alignedSeq)) {
        stop("No alignment found in \"tree\"")
    }
    seqLen <- unique(nchar(alignedSeq))
    if (length(seqLen) != 1)
        stop("Sequence length not equal")
    reference <- attr(tree, "reference")
    qualified <- integer(0)
    alignedSeq <- strsplit(alignedSeq, "")
    for (i in seq_along(reference)) {
        SNP <- table(vapply(
            X = alignedSeq,
            FUN = "[[",
            FUN.VALUE = character(1),
            reference[i]
        ))
        if (sum(SNP > minSNP) >= 2) {
            qualified <- c(qualified, i)
        }
    }
    return(qualified)
}

#' @export
print.sitePath <- function(x, ...) {
    cat("Site",
        attr(x, "site"),
        "may experience fixation on",
        length(x),
        "path(s):\n\n")
    # A 'sitePath' composes of all the fixation paths for a single site.
    #  So each 'm' represent a single fixation path
    for (m in x) {
        if (length(m) == 2) {
            mutName <-
                paste0(attr(m[[1]], "AA"), attr(x, "site"), attr(m[[2]], "AA"))
            cat(mutName,
                paste0("(", length(m[[1]]), "->", length(m[[2]]), ")"),
                "\n")
        } else {
            mutName <- character(0)
            for (tips in m) {
                aa <- attr(tips, "AA")
                mutName <-
                    c(mutName, paste0(aa, "(", length(tips), ")"))
            }
            cat(paste0(mutName, collapse = " -> "), "\n")
        }
    }
    cat("\nIn the bracket are the number of tips",
        "involved before and after the fixation\n")
}

#' @rdname findSites
#' @name fixationSites
#' @description After finding the \code{\link{lineagePath}} of a phylogenetic
#'   tree, \code{fixationSites} uses the result to find those sites that show
#'   fixation on some, if not all, of the lineages. Parallel evolution is
#'   relatively common in RNA virus. There is chance that some site be fixed in
#'   one lineage but does not show fixation because of different sequence
#'   context.
#' @param paths a \code{lineagePath} object returned from
#'   \code{\link{lineagePath}} function or a \code{phylo} object after
#'   \code{\link{addMSA}}
#' @param minEffectiveSize A vector of two integers to specifiy minimum tree
#'   tips involved before/after mutation. Otherwise the mutation will not be
#'   counted into the return. If more than one number is given, the ancestral
#'   takes the first and descendant takes the second as the minimum. If only
#'   given one number, it's the minimum for both ancestral and descendant.
#' @param searchDepth The function uses heuristic search but the termination of
#'   the search cannot be intrinsically decided. \code{searchDepth} is needed to
#'   tell the search when to stop.
#' @param method The strategy for predicting the fixation. The basic approach is
#'   entropy minimization and can be achieved by adding or removing fixation
#'   point, or by comparing the two.
#' @param ... further arguments passed to or from other methods.
#' @examples
#' fixationSites(lineagePath(tree))
#' @return \code{fixationSites} returns a list of fixation mutations with names
#'   of the tips involved.
#' @importFrom utils tail
#' @importFrom stats na.omit
#' @export
fixationSites.lineagePath <- function(paths,
                                      minEffectiveSize = NULL,
                                      searchDepth = 1,
                                      method = c("compare", "insert", "delete"),
                                      ...) {
    tree <- attr(paths, "tree")
    nTips <- length(tree[["tip.label"]])
    align <- attr(paths, "align")
    # Generate the site mapping from reference
    reference <- attr(paths, "reference")
    # Decide which miniminzing strategy
    minimizeEntropy <- switch(
        match.arg(method),
        "compare" = minEntropyByComparing,
        "insert" = minEntropyByInserting,
        "delete" = minEntropyByDeleting
    )
    # Get the 'minEffectiveSize' for each fixation
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- nTips / length(unique(unlist(paths)))
    } else if (!is.numeric(minEffectiveSize)) {
        stop("\"minEffectiveSize\" only accepts numeric")
    }
    minEffectiveSize <- ceiling(minEffectiveSize)
    # Get the 'searchDepth' for heuristic search
    if (searchDepth < 1) {
        stop("\"searchDepth\" should be at least 1")
    } else {
        searchDepth <- ceiling(searchDepth)
    }
    divNodes <- divergentNode(paths)
    nodeAlign <- .tipSeqsAlongPathNodes(
        paths = paths,
        divNodes = divNodes,
        tree = tree,
        align = align
    )
    res <- .findFixationSite(
        paths = paths,
        tree = tree,
        align = align,
        nodeAlign = nodeAlign,
        divNodes = divNodes,
        reference = reference,
        minimizeEntropy = minimizeEntropy,
        minEffectiveSize = minEffectiveSize,
        searchDepth = searchDepth
    )
    attr(res, "paths") <- paths
    attr(res, "reference") <- reference
    class(res) <- "fixationSites"
    return(res)
}

.childrenTips <- function(tree, node) {
    maxTip <- length(tree$tip.label)
    children <- integer(0)
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
    return(getChildren(tree$edge, node))
}

.tipSeqsAlongPathNodes <- function(paths, divNodes, tree, align) {
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
    # Get the sequence of the children tips that are descendant of
    # the nodes. Assign the tip index to the sequences for
    # retrieving the tip name
    nodeAlign <- lapply(nodes, function(n) {
        isTerminal <- FALSE
        if (n %in% terminalNodes) {
            childrenNode <- n
            isTerminal <- TRUE
        } else {
            childrenNode <- tree$edge[which(tree$edge[, 1] == n), 2]
            # Keep the node that is not on the path.
            childrenNode <- setdiff(childrenNode, allNodes)
        }
        tips <- .childrenTips(tree, childrenNode)
        res <- align[tips]
        attr(res, "isTerminal") <- isTerminal
        names(res) <- tips
        return(res)
    })
    # Assign the node names to the 'nodeAlign' list
    names(nodeAlign) <- nodes
    return(nodeAlign)
}

.findFixationSite <- function(paths,
                              tree,
                              align,
                              nodeAlign,
                              divNodes,
                              reference,
                              minimizeEntropy,
                              minEffectiveSize,
                              searchDepth) {
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    # The variable to store the result from entropy minimization for
    # each path with those purely fixed excluded.
    res <- list()
    # Iterate each path
    rootNode <- attr(paths, "rootNode")
    if (!rootNode %in% names(nodeAlign)) {
        paths <- lapply(paths, function(p)
            p[-1])
    }
    for (path in paths) {
        # Iterate every loci (variant sites)
        for (i in loci) {
            site <- as.character(i)
            # The index to use for cpp code
            s <- reference[i] - 1
            # Assign a variable to store the tip names and their info on
            # amino acids. They are the potential fixation segment
            nodeTips <- integer()
            previousAA <- NULL
            currentAA <- NULL
            previousNode <- NULL
            # The input for entropy minimization calculation
            nodeSummaries <- list()
            # Divergent nodes are not included anywhere in the result
            for (node in as.character(setdiff(path, divNodes))) {
                # Get the related descendant tips and related sequences
                nodeTips <- as.integer(names(nodeAlign[[node]]))
                # Frequency of the amino acids at the site
                aaSummary <- tableAA(nodeAlign[[node]], s)
                # Assoicate the amino acid frequence with the tip names
                attr(nodeTips, "aaSummary") <- aaSummary
                # Decide the current fixed amino acid
                if (length(aaSummary) == 1) {
                    currentAA <- names(aaSummary)
                } else {
                    currentAA <- NULL
                }
                # Attach the node to the preivous node if they're
                # both purely fixed and have the same AA fixed.
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
                    # Oddly, R uses the name of the first variable when
                    # adding two numeric vectors. So there is no need for
                    # names (AA) assignment
                }
                # Assign or re-assign the nodeTips with 'aaSummary'
                # to the 'nodeSummaries'
                nodeSummaries[[node]] <- nodeTips
                previousAA <- currentAA
                previousNode <- node
            }
            # Skip to the next locus if AA is fixed along the whole path
            if (length(nodeSummaries) >= 2) {
                seg <- minimizeEntropy(nodeSummaries,
                                       minEffectiveSize,
                                       searchDepth)
                if (length(seg) >= 2) {
                    targetIndex <- length(res[[site]]) + 1
                    attr(seg, "path") <- path
                    res[[site]][[targetIndex]] <- seg
                    attr(res[[site]], "site") <- i
                    attr(res[[site]], "tree") <- tree
                    class(res[[site]]) <- "sitePath"
                }
            }
        }
    }
    return(res)
}

fixationSites.phylo <- function(paths, ...) {
    align <- attr(paths, "align")
    # Generate the site mapping from reference
    reference <- attr(paths, "reference")
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    res <- fixationSitesSearch(nodepath(paths), align, loci)
    res <- res[which(lengths(res) != 1)]
    return(res)
}

treemerBySite <- function(x, ...) {
    align <- attr(x, "align")
    # Generate the site mapping from reference
    reference <- attr(x, "reference")
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    res <- runTreemerBySite(nodepath(x), align, loci)
    return(res)
}

#' @export
fixationSites <- function(paths, ...)
    UseMethod("fixationSites")

#' @export
print.fixationSites <- function(x, ...) {
    cat("Result for", length(attr(x, "paths")), "paths:\n\n")
    if (length(x) == 0) {
        cat("No multi-fixation found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(attr(x, "reference"), "refSeqName")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}
