#' @name plot.lineagePath
#' @title Visualize phylogenetic lineages
#' @description
#' Visualize \code{\link{lineagePath}} object. A tree diagram will be plotted
#' and paths are black solid line while the trimmed nodes and tips will use
#' grey dashed line.
#' @param x A \code{\link{lineagePath}} object
#' @param y
#' Whether plot the nodes from the \code{extendedSearch} in
#' \code{\link{fixationSites}}
#' @param ... Arguments in \code{plot.phylo} functions.
#' @return
#' The function only makes plot and returns no value
#' (It behaviors like the generic \code{\link{plot}} function).
#' @examples
#' data('zikv_tree')
#' data('zikv_align')
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' plot(lineagePath(tree, 0.996))
#' @importFrom ape plot.phylo
#' @importFrom graphics title
#' @export
plot.lineagePath <- function(x, y = TRUE, ...) {
    tree <- attr(x, "tree")
    tree <- ladderize(tree, right = FALSE)
    if (y) {
        x <- .extendPaths(x, tree)
    }
    nEdges <- length(tree$edge.length)
    color <- rep("#d3d3d3", nEdges)
    lty <- rep(2, nEdges)
    width <- rep(1, nEdges)
    targetEdges <- which(tree$edge[, 2] %in% unique(unlist(x)))
    color[targetEdges] <- "#000000"
    lty[targetEdges] <- 1
    width[targetEdges] <- 2
    # TODO: Emphaszie the nodes along the lineagePath
    plot.phylo(tree, edge.color = color, edge.lty = lty, edge.width = width, show.tip.label = FALSE, 
        ...)
}

AA_COLORS <- c(His = "#8282D2", Arg = "#9370DB", Lys = "#145AFF", Ile = "#55AE3A", 
    Phe = "#3232AA", Leu = "#0F820F", Trp = "#B45AB4", Ala = "#C8C8C8", Met = "#FFD700", 
    Pro = "#DC9682", Val = "#2F4F2F", Asn = "#00DCDC", Cys = "#E6E600", Gly = "#666666", 
    Ser = "#FF6347", Tyr = "#ADD8E6", Gln = "#0099CC", Thr = "#FA9600", Glu = "#8C1717", 
    Asp = "#E60A0A", gap = "#000000", unknown = "#d3d3d3", Ile_or_Leu = "#d3d3d3", 
    Asp_or_Asn = "#d3d3d3", Glu_or_Gln = "#d3d3d3")

AA_FULL_NAMES = c(h = "His", r = "Arg", k = "Lys", i = "Ile", f = "Phe", l = "Leu", 
    w = "Trp", a = "Ala", m = "Met", p = "Pro", v = "Val", n = "Asn", c = "Cys", 
    g = "Gly", s = "Ser", y = "Tyr", q = "Gln", t = "Thr", e = "Glu", d = "Asp", 
    `-` = "gap", x = "unknown", j = "Ile_or_Leu", b = "Asp_or_Asn", z = "Glu_or_Gln")

AA_SHORT_NAMES = c(His = "H", Arg = "R", Lys = "K", Ile = "I", Phe = "F", Leu = "L", 
    Trp = "W", Ala = "A", Met = "M", Pro = "P", Val = "V", Asn = "N", Cys = "C", 
    Gly = "G", Ser = "S", Tyr = "Y", Gln = "Q", Thr = "T", Glu = "E", Asp = "D", 
    gap = "-", unknown = "X", Ile_or_Leu = "J", Asp_or_Asn = "B", Glu_or_Gln = "Z")

#' @rdname plotSingleSite
#' @name plotSingleSite
#' @title Color the tree by a single site
#' @description
#' For \code{lineagePath}, the tree will be colored according to the amino
#' acid of the site. The color scheme tries to assign distinguishable
#' color for each amino acid.
#' @param x
#' A \code{fixationSites} object from \code{\link{fixationSites}} or
#' the return from \code{\link{addMSA}} function.
#' @param site
#' One of the mutations in the \code{fixationSites} object. It should
#' be from the \code{\link{names}} of the object. Or an integer to
#' indicate a site could be provide. The numbering is consistent with
#' the reference defined at \code{\link{fixationSites}}.
#' @param showPath
#' If plot the lineage result from lineagePath.
#' @param reference
#' Name of reference for site numbering. The name has to be one of the
#' sequences' name. The default uses the intrinsic alignment numbering.
#' @param gapChar
#' The character to indicate gap. The numbering will skip the gapChar
#' for the reference sequence.
#' @param ...
#' Arguments in \code{plot.phylo} functions and other arguments.
#' @return
#' The function only makes plot and returns no value
#' (It behaviors like the generic \code{\link{plot}} function).
#' @examples
#' data('zikv_tree')
#' data('zikv_align')
#' tree <- addMSA(zikv_tree, alignment = zikv_align)
#' paths <- lineagePath(tree, 0.996)
#' plotSingleSite(paths, 139)
#' @importFrom ape ladderize
#' @importFrom ape getMRCA
#' @importFrom graphics plot
#' @importFrom graphics legend
#' @export
plotSingleSite.lineagePath <- function(x, site, showPath = FALSE, reference = NULL, 
    gapChar = "-", ...) {
    site <- .checkSite(site)
    tree <- attr(x, "tree")
    tree <- ladderize(tree, right = FALSE)
    align <- attr(x, "align")
    align <- strsplit(tolower(align), "")
    refSeqName <- reference
    reference <- .checkReference(tree, align, reference, gapChar)
    tryCatch(site <- match.arg(as.character(site), seq_along(reference)), error = function(e) {
        stop(paste("\"site\":", site, "is not within the length of reference", refSeqName))
    })
    siteComp <- vapply(align, FUN = "[[", FUN.VALUE = character(1), reference[site])
    nEdges <- length(tree$edge.length)
    color <- rep("#000000", nEdges)
    rootNode <- getMRCA(tree, tree$tip.label)
    group <- list()
    for (i in seq_along(siteComp)) {
        group[[siteComp[[i]]]] <- c(group[[siteComp[[i]]]], i)
    }
    AAnames <- AA_FULL_NAMES[names(group)]
    names(group) <- AA_COLORS[AAnames]
    for (g in names(group)) {
        tip2colorEdge(color, g, tree$edge, group[[g]], rootNode)
    }
    width <- rep(1, nEdges)
    if (showPath) {
        targetEdges <- which(tree$edge[, 2] %in% unique(unlist(x)))
        color[targetEdges] <- "#000000"
        width[targetEdges] <- 2
    }
    plot.phylo(tree, show.tip.label = FALSE, edge.color = color, edge.width = width, 
        main = site, ...)
    legend("topleft", title = "Amino acid", legend = unique(AAnames), fill = AA_COLORS[unique(AAnames)], 
        box.lty = 0)
}

#' @rdname plotSingleSite
#' @description
#' For \code{fixationSites}, it will color the ancestral tips in red,
#' descendant tips in blue and excluded tips in grey.
#' @examples
#' \dontrun{
#' fixations <- fixationSites(paths)
#' plotSingleSite(fixations, 139)
#' }
#' @export
plotSingleSite.fixationSites <- function(x, site, ...) {
    site <- .checkSite(site)
    paths <- attr(x, "paths")
    tree <- attr(paths, "tree")
    tree <- ladderize(tree, right = FALSE)
    rootNode <- getMRCA(tree, tree$tip.label)
    tryCatch(site <- match.arg(as.character(site), choices = names(x)), error = function(e) {
        stop(paste("\"site\":", site, "is not a mutation of fixation"))
    })
    sitePaths <- x[[site]]
    plotName <- character(0)
    nEdges <- length(tree$edge.length)
    color <- rep("#d3d3d3", nEdges)
    lty <- rep(2, nEdges)
    width <- rep(1, nEdges)
    AAnames <- character(0)
    for (sp in sitePaths) {
        aa <- toupper(names(sp))
        plotName <- c(plotName, paste0(aa[1], site, aa[2]))
        aa <- AA_FULL_NAMES[tolower(aa)]
        for (n in c(2, 1)) {
            targetEdges <- tip2Edge(tree$edge, sp[[n]], rootNode)
            color[targetEdges] <- AA_COLORS[aa[n]]
            lty[targetEdges] <- 1
            width[targetEdges] <- 2
        }
        AAnames <- c(AAnames, aa)
    }
    plot.phylo(tree, show.tip.label = FALSE, edge.color = color, edge.lty = lty, 
        edge.width = width, main = paste(unique(plotName), collapse = ", "), ...)
    legend("topleft", title = "Amino acid", legend = AA_SHORT_NAMES[unique(AAnames)], 
        fill = AA_COLORS[unique(AAnames)], box.lty = 0)
}

#' @rdname plotSingleSite
#' @description
#' For \code{multiFixationSites}, it will color the tips which have
#' their site fixed. The color will use the same amino acid color
#' scheme as \code{plotSingleSite.lineagePath}
#' @examples
#' \dontrun{
#' multiFixations <- multiFixationSites(paths)
#' plotSingleSite(multiFixations, 1542)
#' }
#' @export
plotSingleSite.multiFixationSites <- function(x, site, ...) {
    site <- .checkSite(site)
    paths <- attr(x, "paths")
    tree <- attr(paths, "tree")
    tree <- ladderize(tree, right = FALSE)
    rootNode <- getMRCA(tree, tree$tip.label)
    tryCatch(site <- match.arg(as.character(site), choices = names(x)), error = function(e) {
        stop(paste("\"site\":", site, "is not a mutation of fixation"))
    })
    sitePaths <- x[[site]]
    plotName <- character(0)
    nEdges <- length(tree$edge.length)
    color <- rep("#d3d3d3", nEdges)
    lty <- rep(2, nEdges)
    width <- rep(0.5, nEdges)
    AAnames <- character(0)
    for (sp in sitePaths) {
        aaName <- character(0)
        for (tips in rev(sp)) {
            aa <- AA_FULL_NAMES[tolower(attr(tips, "AA"))]
            aaName <- c(aa, aaName)
            targetEdges <- tip2Edge(tree$edge, tips, rootNode)
            color[targetEdges] <- AA_COLORS[aa]
            lty[targetEdges] <- 1
            width[targetEdges] <- 2
        }
        AAnames <- c(AAnames, aaName)
        plotName <- c(plotName, paste0(AA_SHORT_NAMES[aaName], collapse = " -> "))
    }
    plot.phylo(tree, show.tip.label = FALSE, edge.color = color, edge.lty = lty, 
        edge.width = width, ...)
    sepChar <- "\n"
    if (sum(nchar(plotName) <= 18)) {
        sepChar <- ",\t"
    }
    title(main = site, sub = paste(plotName, collapse = sepChar))
    legend("topleft", title = "Amino acid", legend = AA_SHORT_NAMES[unique(AAnames)], 
        fill = AA_COLORS[unique(AAnames)], box.lty = 0)
}

#' @export
plotSingleSite <- function(x, site, ...) UseMethod("plotSingleSite")

.checkSite <- function(site) {
    if (!is.numeric(site) || any(site <= 0) || as.integer(site) != site) {
        stop("Please enter positive integer value for \"site\"")
    }
    if (length(site) != 1) {
        warning(paste("\"site\" has more than one element, only the first", site[1], 
            " will be used."))
    }
    return(site[1])
}
