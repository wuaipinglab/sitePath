#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

#' @title Group Tips
#' @description
#' Tree tips are grouped by their sequence similarity and members in a group
#' are constrained to share a same ancestral node.
#' Similarity between two tips is derived from their multiple sequence alignment.
#' The site doesn't count into total length if both are gap.
#' So similarity is calculated as number of matched divided by revised total length
#' @param treeAlignMatch a \code{\link{treeAlignMatch}} object
#' @param similarity similarity threshold for tree trimming
#' @seealso \code{\link{phyloMutations}}
#' @return grouping of tips
#' @export
groupTips <- function(treeAlignMatch, similarity) {
  align <- phyDat2alignment(treeAlignMatch$align)
  grouping <- trimTree(
    lapply(nodepath(treeAlignMatch$tree), rev),
    strsplit(align$seq, ""), 
    similarity
  )
  res <- lapply(grouping, function(g) {
    return(treeAlignMatch$tree$tip.label[g])
  })
  attr(res, "tree") <- treeAlignMatch$tree
  return(res)
}

#' @title Predicted mutations on nodes
#' @description 
#' Under development
#' @param treeAlignMatch a \code{\link{treeAlignMatch}} object
#' @param similarity similarity threshold for tree trimming
#' @return an ancestralMutation object
#' @export
ancestralMutations <- function(treeAlignMatch, similarity, ...) {
  fit <- pml(treeAlignMatch$tree, treeAlignMatch$align, ...)
  fit <- optim.pml(fit, ...)
  anc.ml <- ancestral.pml(fit, return = "phyDat", ...)
  if (length(anc.ml) < max(treeAlignMatch$tree$edge[,1])) {
    stop("The tree needs to be rooted")
  }
  alignAR <- phyDat2alignment(anc.ml)
  res <- mutationList(
    lapply(nodepath(treeAlignMatch$tree), rev),
    strsplit(alignAR$seq, ""), similarity
  )
  attr(res, "tree") <- treeAlignMatch$tree
  attr(res, "class") <- "ancestralMutations"
  return(res)
}

#' @export
print.ancestralMutations <- function(mutations) {
  cat(sum(lengths(mutations)), "mutations on", length(mutations), "nodes.")
}

#' @export


#' @title Predicted mutations on main path
#' @description
#' Sequence of each anchestral node is predict by \code{phangorn}.
#' The tips are grouping by similarity in the same way as \code{\link{groupTips}} does.
#' Common ancestral nodes of each group become the new terminals.
#' The tip stays as terminal if it's the only member in the group.
#' With new terminals, skeleton of the new tree describes the main evolution path.
#' Mutation events that happened on the main path are collected and
#' are categorized as (alternate, coevolve, fixed). In each category is a list.
#' The names are the node number where mutation detected, linked by "~" while
#' the corresponding mutations are under each name.
#' @param treeAlignMatch a \code{\link{treeAlignMatch}} object
#' @param similarity similarity threshold for tree trimming
#' @param model substitution model for reconstruction ancestral sequence
#' @param siteMode specify the mutational mode of return
#' @seealso \code{\link{groupTips}}, \code{\link{pml}}
#' @return a phyloMutations object
#' @export
phyloMutations <- function(
  treeAlignMatch, similarity, siteMode = c(1, 2, 3), ...
) {
  fit <- pml(treeAlignMatch$tree, treeAlignMatch$align, ...)
  fit <- optim.pml(fit, ...)
  anc.ml <- ancestral.pml(fit, type = "ml", return = "phyDat")
  if (length(anc.ml) < max(treeAlignMatch$tree$edge[,1])) {
    stop("The tree needs to be rooted")
  }
  alignAR <- phyDat2alignment(anc.ml)
  res <- mutationPath(
    lapply(nodepath(treeAlignMatch$tree), rev),
    strsplit(alignAR$seq, ""), similarity, siteMode
  )
  attr(res, "tree") <- treeAlignMatch$tree
  attr(res, "class") <- "phyloMutations"
  return(res)
}

#' @export
print.phyloMutations <- function(mutations) {
  cat("Mutation mode(s):", paste(names(mutations), collapse = ", "), "\n")
  print(attr(mutations, "tree"))
}
