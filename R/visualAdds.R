#' @title Mutation map
#' @param mutations the return of \code{\link{phyloMutations}} function
#' @description plot the linked clades and mutations
#' @return NULL
#' @examples 
#' \dontrun{
#' treeFile <- system.file("test.tree", package = "sitePath")
#' alignFile <- system.file("test.aligned.fasta", package = "sitePath")
#' matched <- readTreeAlign(treeFile, "newick", alignFile, "fasta")
#' mutations <- phyloMutations(matched, 0.95, "AA")
#' plot(mutations)
#' }
#' @export

plot.phyloMutations <- function(mutations) {
  for (name in names(mutations)) {
    edges <- strsplit(names(mutations[[name]]), "~")
    g <- graph::graphNEL(nodes = unique(unlist(edges)), edgemode = "directed")
    for (edge in edges) {g <- graph::addEdge(edge[1], edge[2], g)}
    Rgraphviz::plot(
      g, main = name,
      edgeAttrs = list(
        label = sapply(mutations[[name]], function(m) {
          paste(m, collapse = " ")
        })
      ),
      attrs = list(
        node = list(fillcolor = "lightblue"),
        edge = list(color = "cyan")
      )
    )
  }
}

#' @title Internal node plot
#' @param mutations the return of \code{\link{phyloMutations}} function
#' @description plot tree with its internal numbering
#' @return NULL
#' @export

nodePlot <- function(mutations) {
  p <- ggtree::ggtree(attr(mutations, "tree")) + 
    ggtree::geom_text2(ggplot2::aes(subset=!isTip, label=node), hjust=-.1)
  return(p)
}
