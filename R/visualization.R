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
    unknown = "#000000",
    Ile_or_Leu = "#000000",
    Asp_or_Asn = "#000000",
    Glu_or_Gln = "#000000"
)

AA_FULL_NAMES = c(
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
    j = 'Ile_or_Leu',
    b = 'Asp_or_Asn',
    z = 'Glu_or_Gln'
)

#' @name plot.fixationSites
#' @title Color the tree by a single site
#' @description
#' Visualize \code{fixationSites} object. It will plot the tree
#' and color the ancestral tips in red, descendant tips in blue and
#' excluded tips in grey.But the plot will color the tree according to
#' the amino acid instead if \code{site} argument is provided.
#' @param x
#' A \code{fixationSites} object from \code{\link{fixationSites}}
#' @param y
#' One of the mutations in the \code{fixationSites} object.
#' It should be from the \code{\link{names}} of the object.
#' Or an integer indicating a site could be provide. The numbering
#' is consistent with the reference defined at 
#' \code{\link{fixationSites}}.
#' @param ... Arguments in \code{plot.phylo} functions.
#' @return 
#' The function only makes plot and returns no value
#' (It behaviors like the generic \code{\link{plot}} function).
#' @importFrom ape ladderize
#' @importFrom ape getMRCA
#' @examples
#' data("zikv_tree")
#' data("zikv_align")
#' tree <- addMSA(zikv_tree, seqs = zikv_align)
#' fixations <- fixationSites(sitePath(tree, 0.996))
#' plot(fixations, names(fixations)[[1]])
#' @importFrom graphics plot
#' @importFrom graphics legend
#' @export
plot.fixationSites <- function(x, y, ...) {
    if (length(y) != 1) {
        stop("\"y\" is not a single integer or character")
    }
    tree <- attr(x, "tree")
    tree <- ladderize(tree, right = FALSE)
    rootNode <- getMRCA(tree, tree$tip.label)
    if (is.character(y)) {
        tryCatch(
            y <- match.arg(y, choices = names(x)),
            error = function(e) {
                stop("\"y\" is character but not a mutation in fixation")
            }
        )
        color <- rep("#d3d3d3", length(tree$edge.length))
        color <-
            tip2colorEdge(color,
                          "#3F51B5",
                          tree$edge,
                          match(x[[y]][[2]], tree$tip.label),
                          rootNode)
        color <-
            tip2colorEdge(color,
                          "#ff0000",
                          tree$edge,
                          match(x[[y]][[1]], tree$tip.label),
                          rootNode)
        plot(
            tree,
            show.tip.label = FALSE,
            edge.col = color,
            main = y
        )
        legend(
            "topleft",
            title = "Lineages",
            c("ancestral", "descendant", "excluded"),
            fill = c("#ff0000", "#3F51B5", "#d3d3d3"),
            box.lty = 0
        )
    } else if (is.numeric(y)) {
        align <- attr(x, "align")
        align <- strsplit(tolower(align), "")
        reference <- attr(x, "reference")
        tryCatch(
            y <- match.arg(as.character(y), seq_along(reference)),
            error = function(e) {
                stop(paste(
                    "\"y\" is integer but not within the length of reference",
                    attr(x, "refSeqName")
                ))
            }
        )
        siteComp <- vapply(align, FUN = "[[", FUN.VALUE = "", reference[y])
        color <- rep("#FFFF00", length(tree$edge.length))
        group <- list()
        for (i in seq_along(siteComp)) {
            group[[siteComp[[i]]]] <- c(group[[siteComp[[i]]]], i)
        }
        AAnames <- AA_FULL_NAMES[names(group)]
        names(group) <- AA_COLORS[AAnames]
        for (g in names(group)) {
            tip2colorEdge(color, g, tree$edge, group[[g]], rootNode)
        }
        plot(
            tree,
            show.tip.label = FALSE,
            edge.col = color,
            main = y
        )
        legend(
            "topleft",
            title = "Amino acid",
            unique(AAnames),
            fill = AA_COLORS[unique(AAnames)],
            box.lty = 0
        )
    } else {
        stop("\"y\" is neither numeric nor integer type")
    }
}

#' @name plot.sitePath
#' @title Visualize phylogenetic lineages
#' @description
#' Visualize \code{\link{sitePath}} object. A tree diagram will be plotted
#' and paths are black solid line while the trimmed nodes and tips will use
#' grey dashed line.
#' @param x A \code{\link{sitePath}} object
#' @param y
#' Whether plot the nodes from the \code{extendedSearch} in
#' \code{\link{fixationSites}}
#' @param ... Arguments in \code{plot.phylo} functions.
#' @return 
#' The function only makes plot and returns no value
#' (It behaviors like the generic \code{\link{plot}} function).
#' @examples
#' data("zikv_tree")
#' data("zikv_align")
#' tree <- addMSA(zikv_tree, seqs = zikv_align)
#' plot(sitePath(tree, 0.996))
#' @export
plot.sitePath <- function(x, y = TRUE, ...) {
    tree <- attr(x, "tree")
    tree <- ladderize(tree, right = FALSE)
    if (y) {
        x <- extendPaths(x, tree)
    }
    nEdges <- length(tree$edge.length)
    color <- rep("#d3d3d3", nEdges)
    lty <- rep(2, nEdges)
    width <- rep(1, nEdges)
    targetEdges <- which(tree$edge[, 2] %in% unique(unlist(x)))
    color[targetEdges] <- "#000000"
    lty[targetEdges] <- 1
    width[targetEdges] <- 2
    plot(
        tree,
        edge.col = color,
        edge.lty = lty,
        edge.width = width,
        show.tip.label = FALSE,
        ...
    )
}
