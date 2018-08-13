#' @title Tree and Alignment
#' @description
#' Read tree and its matched alignment file. 
#' The return type is the input for other analysis in \code{sitePath} package
#' @param treeFile Path to the file stored phylogenetic tree
#' @param treeFormat Specify tree format ("nexus", "newick", "beast")
#' @param alignFile Path to the file stored sequence alignment
#' @param alignFormat Specify alignment format ("fasta", "clustal", "interleaved", "sequential")
#' @return a \code{\link{treeAlignMatch}} object
#' @export
readTreeAlign <- function(
  treeFile, treeFormat = c("nexus", "newick", "beast"),
  alignFile, alignFormat = c("fasta", "clustal", "interleaved", "sequential")
) {
  tree <- switch(
    treeFormat,
    nexus = ape::read.nexus(file = treeFile),
    newick = ape::read.tree(file = treeFile),
    beast = ggtree::read.beast(file = treeFile)@phylo
  )
  align <- ape::read.dna(
    file = alignFile,
    format = alignFormat,
    as.character = TRUE
  )
  align <- ape::as.alignment(align)
  return(treeAlignMatch(tree, align))
}

#' @title Tree and Alignment
#' @description
#' Assembly a phylo and an alignment object into a treeAlignMatch object.
#' The \code{tip.label} of phylo and \code{nam} of alignment should be matched
#' one by one
#' @param tree an object of \code{phylo}
#' @param align an object of \code{alignment}
#' @return a treeAlignMatch object 
#' @export
treeAlignMatch <- function(tree, align) {
  if (!is(tree, "phylo")) {
    stop("tree is not a phylo object")
  } else if (!is(align, "alignment")) {
    stop("align is not a alignment object")
  }
  align <- toupper(align$seq)[match(tree$tip.label, align$nam)]
  if(any(is.na(align))) {
    stop("tree tips and alignment names are not matched")
  } else if (length(unique(nchar(align))) > 1) {
    stop("Sequence lengths are not the same in alignment")
  }
  return(structure(
    list(tree = tree, align = align), 
    class = "treeAlignMatch"
  ))
}

#' @export
print.treeAlignMatch <- function(matched) {
  cat(unique(nchar(matched$align)), "long alignment\n")
  print(matched$tree)
}
