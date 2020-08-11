#' @rdname fixationSites
#' @name fixationSites
#' @title Fixation sites prediction
#' @description After finding the \code{\link{lineagePath}} of a phylogenetic
#'   tree, \code{fixationSites} uses the result to find those sites that show
#'   fixation on some, if not all, of the lineages. Parallel evolution is
#'   relatively common in RNA virus. There is chance that some site be fixed in
#'   one lineage but does not show fixation because of different sequence
#'   context.
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
    # Collect the fixation and group them by site
    res <- list()
    for (segs in x) {
        for (siteName in names(segs)) {
            seg <- segs[[siteName]]
            # There has to be at least one fixation on the lineage
            if (length(seg) >= 2) {
                res[[siteName]] <- c(res[[siteName]],
                                           list(seg))
                # Assign or re-assign the site and class attributes
                attr(res[[siteName]], "site") <-
                    as.integer(siteName)
                class(res[[siteName]]) <- "sitePath"
            }
        }
    }
    # Merge fixation sites on different paths if applicable
    res <- .combineFixations(res, tree, align)
    # Set 'paths' and 'clustersByPath' attributes
    attr(res, "paths") <- paths
    attr(res, "clustersByPath") <- attr(x, "clustersByPath")
    class(res) <- "fixationSites"
    return(res)
}

.combineFixations <- function(fixations, tree, align) {
    # 'res' is going to be the return of this function. Each entry in
    # the list is the 'sitePath' for a site. Each site ('sitePath')
    # consists of 'mutPath' that is named by the starting node name.
    # The fixed AA and number of non-dominant AA is also stored.
    res <- list()
    for (segs in fixations) {
        i <- attr(segs, "site")
        site <- as.character(i)
        res[[site]][[1]] <- lapply(segs[[1]], function(tips) {
            attr(tips, "tipsAA") <- substr(x = align[tips],
                                           start = i,
                                           stop = i)
            return(tips)
        })
        attr(res[[site]], "site") <- i
        attr(res[[site]], "tree") <- tree
        class(res[[site]]) <- "sitePath"
        for (seg in segs[-1]) {
            # Assume a new fixation path is to add.
            targetIndex <- length(res[[site]]) + 1
            # The index to extract the terminal tips of the fixation.
            endIndex <- length(seg)
            finalAA <- attr(seg[[endIndex]], "AA")
            # The following is to decide if any fixation path can be
            # combined.
            existPath <- res[[site]]
            # The fixation before the temrinal tips should be identical
            # and the final fixed amino acid should be the same.
            toCombine <- vapply(
                X = existPath,
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
                    ... = endIndex
                )))
                allTips <-
                    .childrenTips(tree, getMRCA(tree, toCombine))
                if (all(allTips %in% toCombine)) {
                    seg[[endIndex]] <- toCombine
                    attr(seg[[endIndex]], "AA") <- finalAA
                    res[[site]] <- res[[site]][-existIndex]
                    targetIndex <- length(res[[site]]) + 1
                } else {
                    # Add new fixation for the site if no existing
                    # mutation path can be combined with
                    targetIndex <- length(existPath) + 1
                }
            }
            seg <- lapply(seg, function(tips) {
                attr(tips, "tipsAA") <- substr(x = align[tips],
                                               start = i,
                                               stop = i)
                return(tips)
            })
            res[[site]][[targetIndex]] <- seg
            attr(res[[site]], "site") <- i
            attr(res[[site]], "tree") <- tree
            class(res[[site]]) <- "sitePath"
        }
    }
    return(res)
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

#' @export
fixationSites <- function(paths, ...)
    UseMethod("fixationSites")

#' @export
print.fixationSites <- function(x, ...) {
    cat("This is a 'fixationSites' object.\n\nResult for",
        length(attr(x, "paths")),
        "paths:\n\n")
    if (length(x) == 0) {
        cat("No multi-fixation found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(x, "reference")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}

#' @rdname plotFunctions
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tidytree groupOTU
#' @importFrom ggplot2 scale_color_manual guides guide_legend
#' @importFrom ggrepel geom_label_repel
#' @export
#' @examples
#' data(zikv_tree)
#' data(zikv_align)
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' paths <- lineagePath(tree)
#' plot(paths)
#' fixations <- fixationSites(paths)
#' plot(fixations)
#' x <- fixationPath(fixations)
#' plot(x)
plot.fixationSites <- function(x,
                               y = TRUE,
                               ...) {
    tree <- as.treedata.fixationSites(x)
    grp <- levels(attr(tree@phylo, "Groups"))
    grp <- grp[which(grp != "0")]
    groupColors <-
        colorRampPalette(brewer.pal(9, "Set1"))(length(grp))
    names(groupColors) <- grp
    groupColors["0"] <- "black"

    p <- ggtree(tree, aes(color = Groups)) +
        scale_color_manual(values = groupColors) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        theme(legend.position = "left")
    if (y) {
        p <- p + geom_label_repel(
            aes(x = branch, label = SNPs),
            fill = "lightgreen",
            color = "black",
            min.segment.length = 0,
            na.rm = TRUE
        )
    }
    return(p)
}
