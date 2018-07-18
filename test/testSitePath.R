tree <- ggtree::read.beast(
  system.file("m.trees", package = "sitePath")
)@phylo

align <- read.phyDat(
  system.file("m.aligned.fasta", package = "sitePath"),
  format = "fasta", type = "DNA"
)

matched <- treeAlignMatch(tree, align)
grouping <- groupTips(matched, 0.95)
mutations <- phyloMutations(matched, similarity = 0.95)
tipName <- toTips(mutations)
plot(mutations)
nodePlot(mutations)
