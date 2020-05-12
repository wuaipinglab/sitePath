#' @rdname addMSA
#' @name addMSA
#' @title Prepare data for sitePath analysis
#' @description sitePath requires both tree and sequence alignment to do the
#'   analysis. \code{addMSA} wraps \code{read.alignment} function in
#'   \code{seqinr} package and helps match names in tree and sequence alignment.
#'   Either provide the file path to an alignment file and its format or an
#'   alignment object from the return of \code{read.alignment} function. If both
#'   the file path and alignment object are given, the function will use the
#'   sequence in the alignment file.
#' @param tree a \code{phylo} object. This commonly can be from tree parsing
#'   function in \code{ape} or \code{ggtree}. All the \code{tip.label} should be
#'   found in the sequence alignment.
#' @param msaPath The file path to the multiple sequence alignment file
#' @param msaFormat The format of the multiple sequence alignment file
#' @param alignment an \code{alignment} object. This commonly can be from
#'   sequence parsing function in the \code{seqinr} package. Sequence names in
#'   the alignment should include all \code{tip.label} in the tree
#' @return \code{addMSA} returns a \code{phylo} object with matched multiple
#'   sequence alignment
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
#' @importFrom seqinr read.alignment
#' @importFrom methods is
#' @export
addMSA <- function(tree,
                   msaPath = "",
                   msaFormat = "",
                   alignment = NULL) {
    # Get/test the tree object
    if (is(tree, "treedata")) {
        tree <- tree@phylo
    }
    if (!inherits(tree, "phylo")) {
        stop("\"tree\" is not class phylo")
    }
    # Read alignment from the file if applicable
    if (file.exists(msaPath)) {
        alignment <- read.alignment(msaPath, msaFormat)
    } else if (is.null(alignment)) {
        stop("Alignment file \"", msaPath, "\" does not exist")
    }
    # Test the alignment object
    if (!is(alignment, "alignment")) {
        stop("\"alignment\" is not class alignment")
    }
    # Map the names between tree and alignment
    m <- match(tree$tip.label, alignment$nam)
    if (any(is.na(m))) {
        stop("Tree tips and alignment names are not matched")
    } else {
        align <- toupper(alignment$seq[m])
        if (length(unique(nchar(align))) > 1) {
            stop("Sequence lengths are not the same in alignment")
        }
    }
    attr(tree, "align") <- align
    # Use the numbering of MSA as the default site numbering
    attr(tree, "reference") <-
        .checkReference(tree, align, NULL, "-")
    return(tree)
}

.checkReference <- function(tree, align, reference, gapChar) {
    if (is.null(reference)) {
        reference <- seq_len(nchar(align[1]))
    } else {
        if (!is.character(gapChar) ||
            nchar(gapChar) != 1 || length(gapChar) != 1) {
            stop("\"gapChar\" only accepts one single character")
        }
        reference <-
            getReference(align[which(tree$tip.label == reference)], gapChar)
    }
    return(reference)
}

#' @rdname setSiteNumbering
#' @name setSiteNumbering
#' @title Set site numbering to the reference sequence
#' @description A reference sequence can be used to define a global site
#'   numbering scheme for multiple sequence alignment. The gap in the reference
#'   will be skipped so the site ignored in numbering.
#' @param x The object to set site numbering. It could be a \code{phylo} object
#'   after \code{\link{addMSA}} or a \code{lineagePath} object. The function for
#'   \code{fixaitonSites} and \code{multiFixationSites} will be added in later
#'   version.
#' @param reference Name of reference for site numbering. The name has to be one
#'   of the sequences' name. The default uses the intrinsic alignment numbering
#' @param gapChar The character to indicate gap. The numbering will skip the
#'   gapChar for the reference sequence.
#' @param ... further arguments passed to or from other methods.
#' @return A \code{phylo} object with site numbering mapped to reference
#'   sequence
#' @examples
#' data(zikv_tree)
#' msaPath <- system.file('extdata', 'ZIKV.fasta', package = 'sitePath')
#' tree <- addMSA(zikv_tree, msaPath = msaPath, msaFormat = 'fasta')
#' setSiteNumbering(tree)
#' @export
setSiteNumbering.phylo <- function(x,
                                   reference = NULL,
                                   gapChar = "-",
                                   ...) {
    align <- attr(x, "align")
    siteMapping <- .checkReference(x, align, reference, gapChar)
    attr(siteMapping, "refSeqName") <- reference
    attr(x, "reference") <- siteMapping
    return(x)
}

#' @rdname setSiteNumbering
#' @name setSiteNumbering
#' @export
setSiteNumbering.lineagePath <- function(x,
                                         reference = NULL,
                                         gapChar = "-",
                                         ...) {
    align <- attr(x, "align")
    tree <- attr(x, "tree")
    siteMapping <-
        .checkReference(tree, align, reference, gapChar)
    attr(siteMapping, "refSeqName") <- reference
    attr(x, "reference") <- siteMapping
    return(x)
}

#' @export
setSiteNumbering <- function(x, reference, gapChar, ...)
    UseMethod("setSiteNumbering")

# TODO: Need a class called "reference" for site numbering

##' @importFrom ape as.phylo
##' @export
ape::as.phylo

#' @rdname as.phylo.fixationSites
#' @name as.phylo.fixationSites
#' @title Represent site fixation as a tree
#' @description The tips are clustered according to the fixation sites. The
#'   transition of fixation sites will be plotted as a phylogenetic tree. The
#'   length of each branch represents the number of fixation mutation between
#'   two clusters. The name of the tree tips indicate the number of sequences in
#'   the cluster.
#' @param x The return from \code{\link{fixationSites}} function.
#' @param minEffectiveSize The minimum size for a tip cluster
#' @param ... Further arguments passed to or from other methods.
#' @return A \code{phylo} object
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' paths <- lineagePath(tree)
#' mutations <- fixationSites(paths)
#' as.phylo(mutations)
#' @method as.phylo fixationSites
#' @export
as.phylo.fixationSites <- function(x, minEffectiveSize = NULL, ...) {
    grouping <- .transitionClusters(x)
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <-
            mean(lengths(unlist(grouping, recursive = FALSE)))
    }
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

        # The initial tips of the group
        currentTips <- gp[[1]]
        # Assume the initial reference site
        refSites <- attr(currentTips, "site")
        # Find where to merge and parent node, update reference site maybe
        for (i in seq_along(tipClusters)) {
            toMerge <- tipClusters[[i]]
            toMergeIndex <- attr(toMerge, "toMerge")
            if (!is.null(toMergeIndex) && toMergeIndex == gpIndex) {
                refSites <- attr(toMerge, "toMergeRefSites")
                parentNode <- parentNodes[which(childrenNodes == i)]
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
            edgeIndex <- which(childrenNodes == mostRelatedTipNode)
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
    res <- list(
        "edge" = cbind(parentNodes, childrenNodes),
        "edge.length" = lengths(edgeSNPs),
        "Nnode" = length(unique(parentNodes)),
        "tip.label" = as.character(lengths(unlist(
            grouping, recursive = FALSE
        )))
    )
    attr(res, "tipClusters") <- tipClusters
    attr(res, "edgeSNPs") <-
        lapply(seq_along(edgeSNPs), function(i) {
            snp <- edgeSNPs[[i]]
            attr(snp, "edge") <- res[["edge"]][i, ]
            return(snp)
        })
    class(res) <- "phylo"
    return(res)
}

.transitionClusters <- function(fixations) {
    paths <- attr(fixations, "paths")
    tree <- attr(paths, "tree")
    groupByPath <- lapply(paths, function(p) {
        terminalTips <- .childrenTips(tree, p[length(p)])
        # Group fixation results by path rather than site
        group <- list()
        for (sp in fixations) {
            site <- attr(sp, "site")
            for (mp in sp) {
                tips <- mp[[length(mp)]]
                # Filter if not belong to the path
                if (all(terminalTips %in% tips)) {
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
                    if (length(common) == 0) {
                        newGrouping <- res[1:i]
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
                            newGrouping <-
                                c(newGrouping,
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
    # Find the divergent point and remove the overlapped part
    grouping <- list(groupByPath[[1]])

    for (gpIndex in seq_along(groupByPath)[-1]) {
        gp <- groupByPath[[gpIndex]]
        toMergeIndex <- NULL
        divergedIndex <- 0L
        # Loop through to find the most related group
        for (i in seq_along(grouping)) {
            allTips <- unlist(grouping[[i]])
            for (j in seq_along(gp)[-1]) {
                if (all(!gp[[j]] %in% allTips)) {
                    m <- i
                    d <- j
                    break
                }
            }
            if (d > divergedIndex) {
                toMergeIndex <- m
                divergedIndex <- d
            }
        }
        # Find the tips before diverged
        sharedTips <- gp[[divergedIndex - 1]]
        refSites <- attr(sharedTips, "site")
        # The non-shared part
        divergedTips <- setdiff(sharedTips, allTips)
        attr(divergedTips, "site") <- refSites
        # Drop the overlapped part
        grouping[[gpIndex]] <-
            c(list(divergedTips), gp[divergedIndex:length(gp)])
        # Find the most related group
        toMerge <- grouping[[toMergeIndex]]
        # To determine where to add the new (truncated) group
        for (i in seq_along(toMerge)) {
            gpTips <- unlist(gp)
            if (all(!toMerge[[i]] %in% gpTips)) {
                # Find the tips before diverged
                sharedTips <- toMerge[[i - 1]]
                sites <- attr(sharedTips, "site")
                # The non-shared part
                divergedTips <- setdiff(sharedTips, gpTips)
                attr(divergedTips, "site") <- sites
                # The shared part
                sharedTips <- setdiff(sharedTips, divergedTips)
                attr(sharedTips, "site") <- sites
                attr(sharedTips, "toMerge") <- gpIndex
                attr(sharedTips, "toMergeRefSites") <- refSites
                # Reform
                if (i == 2) {
                    preTips <- list()
                } else {
                    preTips <- toMerge[1:(i - 2)]
                }
                grouping[[toMergeIndex]] <- c(preTips,
                                              list(sharedTips),
                                              list(divergedTips),
                                              toMerge[i:length(toMerge)])
                break
            }
        }
    }
    return(grouping)
}

#' @name zikv_align
#' @title Multiple sequence alignment of Zika virus polyprotein
#' @description The raw protein sequences were downloaded from ViPR database
#'   (\url{https://www.viprbrc.org/}) and aliged using MAFFT. with default
#'   settings.
#' @format a \code{alignment} object
#' @usage data(zikv_align)
#' @docType data
"zikv_align"

#' @name zikv_tree
#' @title Phylogenetic tree of Zika virus polyprotein
#' @description Tree was built from \code{\link{zikv_align}} using RAxML with
#'   default settings.The tip ANK57896 was used as outgroup to root the tree.
#' @format a \code{phylo} object
#' @usage data(zikv_tree)
#' @docType data
"zikv_tree"

#' @name h3n2_align
#' @title Multiple sequence alignment of H3N2's HA protein
#' @description The raw protein sequences were downloaded from NCBI database.
#' @format a \code{alignment} object
#' @usage data(h3n2_align)
#' @docType data
"h3n2_align"

#' @name h3n2_tree
#' @title Phylogenetic tree of H3N2's HA protein
#' @description Tree was built from \code{\link{h3n2_align}} using RAxML with
#'   default settings.
#' @format a \code{phylo} object
#' @usage data(h3n2_tree)
#' @docType data
"h3n2_tree"

#' @name zikv_align_reduced
#' @title Truncated data for runnable example
#' @description This is a truncated version of \code{\link{zikv_align}}
#' @format a \code{alignment} object
#' @usage data(zikv_align_reduced)
#' @docType data
"zikv_align_reduced"

#' @name zikv_tree_reduced
#' @title Truncated data for runnable example
#' @description This is a truncated version of \code{\link{zikv_tree}}
#' @format a \code{phylo} object
#' @usage data(zikv_tree_reduced)
#' @docType data
"zikv_tree_reduced"

#' @name h3n2_align_reduced
#' @title Truncated data for runnable example
#' @description This is a truncated version of \code{\link{h3n2_align}}
#' @format a \code{alignment} object
#' @usage data(h3n2_align_reduced)
#' @docType data
"h3n2_align_reduced"

#' @name h3n2_tree_reduced
#' @title Truncated data for runnable example
#' @description This is a truncated version of \code{\link{h3n2_tree}}
#' @format a \code{phylo} object
#' @usage data(h3n2_tree_reduced)
#' @docType data
"h3n2_tree_reduced"

#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

AA_COLORS <- c(
    His = "#8282D2",
    Arg = "#9370DB",
    Lys = "#145AFF",
    Ile = "#55AE3A",
    Phe = "#3232AA",
    Leu = "#0F820F",
    Trp = "#B45AB4",
    Ala = "#C8C8C8",
    Met = "#FFD700",
    Pro = "#DC9682",
    Val = "#2F4F2F",
    Asn = "#00DCDC",
    Cys = "#E6E600",
    Gly = "#666666",
    Ser = "#FF6347",
    Tyr = "#ADD8E6",
    Gln = "#0099CC",
    Thr = "#FA9600",
    Glu = "#8C1717",
    Asp = "#E60A0A",
    gap = "#000000",
    unknown = "#d3d3d3",
    Ile_or_Leu = "#d3d3d3",
    Asp_or_Asn = "#d3d3d3",
    Glu_or_Gln = "#d3d3d3"
)

AA_FULL_NAMES <- c(
    h = "His",
    r = "Arg",
    k = "Lys",
    i = "Ile",
    f = "Phe",
    l = "Leu",
    w = "Trp",
    a = "Ala",
    m = "Met",
    p = "Pro",
    v = "Val",
    n = "Asn",
    c = "Cys",
    g = "Gly",
    s = "Ser",
    y = "Tyr",
    q = "Gln",
    t = "Thr",
    e = "Glu",
    d = "Asp",
    `-` = "gap",
    x = "unknown",
    j = "Ile_or_Leu",
    b = "Asp_or_Asn",
    z = "Glu_or_Gln"
)

AA_SHORT_NAMES <- c(
    His = "H",
    Arg = "R",
    Lys = "K",
    Ile = "I",
    Phe = "F",
    Leu = "L",
    Trp = "W",
    Ala = "A",
    Met = "M",
    Pro = "P",
    Val = "V",
    Asn = "N",
    Cys = "C",
    Gly = "G",
    Ser = "S",
    Tyr = "Y",
    Gln = "Q",
    Thr = "T",
    Glu = "E",
    Asp = "D",
    gap = "-",
    unknown = "X",
    Ile_or_Leu = "J",
    Asp_or_Asn = "B",
    Glu_or_Gln = "Z"
)
