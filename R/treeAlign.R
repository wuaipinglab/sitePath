#' @title Tree and Alignment
#' @description
#' Read tree and its matched alignment file. 
#' The return type is the input for other analysis in \code{sitePath} package
#' @param treeFile Path to the file stored phylogenetic tree
#' @param treeFormat Specify tree format ("nexus", "newick", "beast")
#' @param alignFile Path to the file stored sequence alignment
#' @param alignFormat Specify alignment format ("fasta", "clustal", "phylip", "mase", "msf")
#' @param seqType Sequence type ("DNA", "AA")
#' @return a \code{\link{treeAlignMatch}} object
#' @export
readTreeAlign <- function(
  treeFile, treeFormat = c("nexus", "newick", "beast"),
  alignFile, alignFormat = c("fasta", "clustal", "phylip", "mase", "msf"),
  seqType = c("DNA", "AA")
) {
  tree <- switch (
    treeFormat,
    nexus = ape::read.nexus(file = treeFile),
    newick = ape::read.tree(file = treeFile),
    beast = ggtree::read.beast(file = treeFile)@phylo
  )
  align <- phangorn::read.phyDat(
    file = alignFile,
    format = alignFormat,
    type = seqType
  )
  if(any(is.na(match(tree$tip.label, attr(align, "names"))))) {
    align <- seqinr::read.alignment(alignFile, format = alignFormat)
    align <- phangorn::as.phyDat(align, type = seqType)
  }
  return(treeAlignMatch(tree, align))
}

#' @title Tree and Alignment
#' @description
#' Assembly a phylo and a phyDat object into a treeAlignMatch object.
#' The \code{tip.label} of phylo and \code{names} of phyDat should be matched
#' one by one
#' @param tree an object of \code{phylo}
#' @param align an object of \code{phyDat}
#' @return a treeAlignMatch object 
#' @export
treeAlignMatch <- function(tree, align) {
  if (!is(tree, "phylo")) {
    stop("tree is not a phylo object")
  } else if (!is(align, "phyDat")) {
    stop("align is not a phyDat object")
  } else if(any(is.na(match(tree$tip.label, attr(align, "names"))))) {
    stop("tree tips and alignment names are not matched")
  }
  align = phangorn::phyDat2alignment(subset(align, tree$tip.label))
  return(structure(
    list(tree = tree, align = align), 
    class = "treeAlignMatch"
  ))
}

#' @export
print.treeAlignMatch <- function(matched) {
  print(matched$align)
  print(matched$tree)
}
