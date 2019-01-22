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

#' @rdname Visualization
#' @name plot.fixationSites
#' @title Color the tree by a single site
#' @description
#' The result of \code{\link{fixationSites}} can be visualized using
#' the funciton. It will plot the tree and color the ancestral tips
#' in red, descendant tips in blue and excluded tips in grey.
#' But the plot will color the tree according to the amino acid instead
#' if \code{site} argument is provided.
#' @param x
#' a \code{fixationSites} object from \code{\link{fixationSites}}
#' @param y
#' One of the mutations in the \code{\link{fixationSites}} object.
#' It should be from the \code{\link{names}} of the object.
#' Or an integer indicating a site could be provide. The numbering
#' is consistent with the reference defined at \code{\link{fixationSites}}.
#' @param ... further arguments passed to or from other methods.
#' @importFrom ape ladderize
#' @importFrom ape getMRCA
#' @examples
#' data("zikv_fixations")
#' plot(zikv_fixations, names(zikv_fixations)[[1]])
#' @return plot and color the tree
#' @importFrom graphics plot
#' @importFrom graphics legend
#' @export
plot.fixationSites <- function(x, y, ...) {
    if (length(y) != 1) {
        stop("site is not a single integer or character")
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
            y <- match.arg(as.character(y), 1:length(reference)),
            error = function(e) {
                stop(paste(
                    "site is integer but not within the length of reference",
                    attr(x, "refSeqName")
                ))
            }
        )
        siteComp <- sapply(align, "[[", reference[y])
        color <- rep("#FFFF00", length(tree$edge.length))
        group <- list()
        for (i in 1:length(siteComp)) {
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
        stop("site is neither numeric nor integer type")
    }
}

# plot.fixationSites <- function(x, y, ...) {
#     tree <- attr(x, "tree")
#     if (length(y) != 1) {
#         stop("site is not a single integer or character")
#     }
#     if (is.character(y)) {
#         tryCatch(
#             y <- match.arg(y, choices = names(x)),
#             error = function(e) {
#                 stop("\"y\" is character but not a mutation in fixationSites")
#             }
#         )
#         if (!suppressWarnings(require("ggtree")))
#             group <- c(x[[y]],
#                        list(excluded = setdiff(tree$tip.label, unlist(x[[y]]))))
#         groupCols <- c("#d3d3d3", "#ff0000", "#3F51B5")
#         names(groupCols) <- c("excluded", "ancestral", "descendant")
#         return(plotColoredTree(tree, group, groupCols, y, "Lineage"))
#     } else if (is.numeric(y)) {
#         align <- attr(x, "align")
#         align <- strsplit(tolower(align), "")
#         reference <- attr(x, "reference")
#         tryCatch(
#             y <- match.arg(as.character(y), 1:length(reference)),
#             error = function(e) {
#                 stop(paste(
#                     "site is integer but not within the length of reference",
#                     attr(x, "refSeqName")
#                 ))
#             }
#         )
#         y <- reference[y]
#         siteComp <- sapply(align, "[[", y)
#         names(siteComp) <- tree$tip.label
#         group <- list()
#         for (s in names(siteComp)) {
#             group[[siteComp[[s]]]] <- c(group[[siteComp[[s]]]], s)
#         }
#         names(group) <- AA_FULL_NAMES[names(group)]
#         groupCols <- AA_COLORS[names(group)]
#         return(plotColoredTree(tree, group, groupCols, y, "Amino acid"))
#     } else {
#         stop("site is neither numeric nor integer type")
#     }
# }
#
# plotColoredTree <-
#     function(tree, group, groupCols, title, legendTitle) {
#         p <- ggtree(groupOTU(tree, group), aes(color = group)) +
#             scale_color_manual(values = c('white', groupCols)) +
#             ggtitle(title) +
#             theme(legend.position = "left") +
#             guides(color = guide_legend(override.aes = list(size = 3),
#                                         title = legendTitle))
#         return(p)
#     }
