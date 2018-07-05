#' @param treeFile Path to the file stored tree string
#' @param treeFormat Specify tree format
#' @param alignFile Path to the file stored sequence alignment
#' @param alignFormat Specify alignment format
#' @return a treeAlignMatch object
#' @export

readTreeAlign <- function(
  treeFile, treeFormat = c("nexus", "newick", "beast"),
  alignFile, alignFormat = c("fasta", "clustal", "phylip", "mase", "msf")
) {
  tree <- switch (
    treeFormat,
    nexus = read.nexus(treeFile),
    newick = read.tree(treeFile),
    beast = ggtree::read.beast(treeFile)
  )
  align <- read.alignment(alignFile, alignFormat)
  if(!identical(sort(tree$tip.label), sort(align$nam))) {
    stop("tree tips and alignment names are not matched")
  }
  return(structure(
    list(tree = tree, align = align), 
    class = "treeAlignMatch"
  ))
}

#' @param tree a phylo object
#' @param align an alignment object
#' @param similarity similarity threshold for tree trimming
#' @return grouping of tips
#' @export

groupTips <- function(tree, align, similarity) {
  return(trimTree(
    lapply(nodepath(tree), rev),
    strsplit(unlist(align$seq), "")[match(tree$tip.label, align$nam)], 
    similarity
  ))
}

#' @param tree a phylo object
#' @param align an alignment object
#' @param similarity similarity threshold for tree trimming
#' @param seqType specify sequence type ("DNA", "AA")
#' @param model substitution model for reconstruction ancestral sequence
#' @param siteMode specify the mutational mode of return
#' @return Linked clades and the mutations
#' @export

ancestralMutations <- function(
  tree, align, similarity, seqType = c("DNA", "AA"),
  model = NULL, siteMode = c(1, 2)
) {
  fit <- pml(tree, as.phyDat(align, seqType))
  if (is.null(model)) {
    model <- switch (seqType, DNA = "GTR", AA = "JTT")
  }
  fit <- optim.pml(fit, model = model)
  anc.ml <- ancestral.pml(fit, type = "ml", return = "phyDat")
  alignAR <- phyDat2alignment(anc.ml)
  return(mutationPath(
    lapply(nodepath(tree), rev),
    strsplit(alignAR$seq, ""),
    similarity, siteMode
  ))
}

#' @param tree a phylo object
#' @param mutations the return of \code{ancestralMutations} function
#' @return clade and corresponding tips
#' @export

mutations2tips <- function(tree, mutations) {
  res <- lapply(mutations, function(m) {
    edges <- lapply(strsplit(names(m), "~"), function(e) {
      tips <- lapply(as.integer(e), function(n) {
        return(Descendants(tree, n, "tips"))
      })
      names(tips) <- e
      return(tips)
    })
    names(edges) <- names(m)
    return(edges)
  })
  return(res)
}

#' @param mutations the return of \code{ancestralMutations} function
#' @description plot the linked clades and mutations
#' @return NULL

mutations2graphviz <- function(mutations) {
  for (name in names(mutations)) {
    edges <- strsplit(names(mutations[[name]]), "~")
    g <- graph::graphNEL(nodes = unique(unlist(edges)), edgemode = "directed")
    for (edge in edges) {g <- graph::addEdge(edge[1], edge[2], g)}
    Rgraphviz::plot(
      g, edgeAttrs = list(
        label = sapply(mutations[[name]], function(m) {
          paste(m, collapse = " ")
        })
      ),
      attrs = list(
        node=list(fillcolor="lightgreen"),
        edge=list(color="cyan")
        # graph=list(rankdir="LR")
      )
    )
  }
}
