library(sitePath)

matched <- readTreeAlign(
  system.file("test.tree", package = "sitePath"), "newick",
  system.file("test.aligned.fasta", package = "sitePath"), "fasta"
)
grouping <- groupTips(matched, 0.95)
mutations <- ancestralMutations(matched, 0.95, "AA")
nodeName <- mutations2Nodes(mutations)
mutations2graphviz(mutations)
