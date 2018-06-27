library(phangorn)
library(seqinr)
library(Rgraphviz)
library(ggtree)
library(Rcpp)
sourceCpp("treemer.cpp")

groupTips <- function(tree, align, similarity = 0.95) {
  if (!is(tree, "phylo")) {
    stop("object 'tree' is not of class 'phylo'")
  } else if (!is(align, "alignment")) {
    stop("object 'align' is not of class 'alignment")
  }
  tipPath <- lapply(nodepath(tree), rev)
  seqs <- strsplit(
    sapply(align$seq, function(seq) {seq}), ""
  )[match(tree$tip.label, align$nam)]
  if (anyNA(seqs)) {
    stop("tree tips and alignment names are not matched")
  }
  return(trimTree(tipPath, seqs, similarity))
}

ancestralMutations <- function(
  tree, align, seqType, model,
  similarity = 0.95
) {
  if (!is(tree, "phylo")) {
    stop("object 'tree' is not of class 'phylo'")
  } else if (!is(align, "alignment")) {
    stop("object 'align' is not of class 'alignment")
  } else if(!identical(sort(tree$tip.label), sort(align$nam))) {
    stop("tree tips and alignment names are not matched")
  }
  tipPath <- lapply(nodepath(tree), rev)
  seqsAR <- as.phyDat(align, seqType)
  fit <- pml(tree, seqsAR)
  fit <- optim.pml(fit, model = model)
  anc.ml <- ancestral.pml(fit, type = "ml", return = "phyDat")
  alignAR <- phyDat2alignment(anc.ml)
  seqsAR <- strsplit(alignAR$seq, "")
  mutations <- mutationPath(tipPath, seqsAR, similarity)
  return(mutations)
}

filterMutation <- function(mutations) {
  
}

mutations2graphviz <- function(mutations) {
  edges <- strsplit(names(mutations), "~")
  g <- graphNEL(nodes = unique(unlist(edges)), edgemode = "directed")
  for (edge in edges) {g <- addEdge(edge[1], edge[2], g)}
  plot(
    g, edgeAttrs = list(
      label = sapply(mutations, function(m) {paste(m, collapse = ",")})
    ),
    attrs = list(
      node=list(fillcolor="lightgreen"),
      edge=list(color="cyan")
      # graph=list(rankdir="LR")
    )
  )
}
