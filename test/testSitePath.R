matched <- readTreeAlign(
  system.file("m.trees", package = "sitePath"), "beast",
  system.file("m.aligned.fasta", package = "sitePath"), "fasta",
  "AA"
)
grouping <- groupTips(matched, 0.95)
mutations <- phyloMutations(matched, 0.95)
tipName <- toTips(mutations)
plot(mutations)
nodePlot(mutations)
