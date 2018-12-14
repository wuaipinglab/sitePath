AA_COLORS <- c(
  His = "#8282D2", Arg = "#9370DB", Lys = "#145AFF",
  Ile = "#55AE3A", Phe = "#3232AA", Leu = "#0F820F",
  Trp = "#B45AB4", Ala = "#C8C8C8", Met = "#FFD700",
  Pro = "#DC9682", Val = "#2F4F2F", Asn = "#00DCDC",
  Cys = "#E6E600", Gly = "#666666", Ser = "#FF6347",
  Tyr = "#ADD8E6", Gln = "#0099CC", Thr = "#FA9600",
  Glu = "#8C1717", Asp = "#E60A0A", gap = "#000000",
  unknown = "#000000", Ile_or_Leu = "#000000",
  Asp_or_Asn = "#000000", Glu_or_Gln = "#000000"
)

AA_FULL_NAMES = c(
  h = "His", r = "Arg", k = "Lys",
  i = "Ile", f = "Phe", l = "Leu",
  w = "Trp", a = "Ala", m = "Met",
  p = "Pro", v = "Val", n = "Asn",
  c = "Cys", g = "Gly", s = "Ser",
  y = "Tyr", q = "Gln", t = "Thr",
  e = "Glu", d = "Asp", `-` = "gap",
  x = "unknown", j = 'Ile_or_Leu',
  b = 'Asp_or_Asn', z = 'Glu_or_Gln'
)

#' @rdname Visualization
#' @name plot.fixationSites
#' @title Color the tree by a single site
#' @description 
#' The result of \code{\link{fixationSites}} can be visualized using the funciton.
#' It will plot the tree and color the tips before fixation in red, 
#' tips after fixation in blue and excluded tips in grey. 
#' But the plot will color the tree according to the amino acid instead
#' if \code{site} argument is provided.
#' @param fixationSites a \code{fixationSites} object from \code{\link{fixationSites}}
#' @param site 
#' One of the mutations in the \code{\link{fixationSites}} object. It should be from 
#' the \code{\link{names}} of the object.Or an integer indicating a site could be provide.
#' The numbering is consistent with the reference defined at \code{\link{fixationSites}}.
#' @return a ggtree plot
#' @import ggplot2
#' @import ggtree
#' @export
plot.fixationSites <- function(fixationSites, site) {
  tree <- attr(fixationSites, "tree")
  if (length(site) != 1) {stop("site is not a single integer or character")}
  if (is.character(site)) {
    tryCatch(
      site <- match.arg(site, choices = names(fixationSites)),
      error = function(e) {
        stop("site is character but not a mutation in fixationSites")
      }
    )
    group <- c(
      fixationSites[[site]],
      list(excluded = setdiff(tree$tip.label, unlist(fixationSites[[site]])))
    )
    groupCols <- c("#d3d3d3", "#ff0000", "#3F51B5")
    names(groupCols) <- c("excluded", "ancestral", "descendant")
    return(plotColoredTree(tree, group, groupCols, site, "Lineage"))
  } else if (is.numeric(site)) {
    align <- attr(fixationSites, "align")
    align <- strsplit(tolower(align), "")
    reference <- attr(fixationSites, "reference")
    tryCatch(
      site <- match.arg(as.character(site), 1:length(reference)),
      error = function(e) {
        stop(paste(
          "site is integer but not within the length of reference", 
          attr(fixationSites, "refSeqName")
        ))
      }
    )
    site <- reference[site]
    siteComp <- sapply(align, "[[", site)
    names(siteComp) <- tree$tip.label
    group <- list()
    for (s in names(siteComp)) {
      group[[siteComp[[s]]]] <- c(group[[siteComp[[s]]]], s)
    }
    names(group) <- AA_FULL_NAMES[names(group)]
    groupCols <- AA_COLORS[names(group)]
    return(plotColoredTree(tree, group, groupCols, site, "Amino acid"))
  } else {
    stop("site is neither numeric nor integer type")
  }
}

plotColoredTree <- function(tree, group, groupCols, title, legendTitle) {
  p <- ggtree(groupOTU(tree, group), aes(color = group)) +
    scale_color_manual(values = c('white', groupCols)) +
    ggtitle(title) +
    theme(legend.position = "left") +
    guides(
      color = guide_legend(
        override.aes = list(size = 3),
        title = legendTitle
      )
    )
  return(p)
}
