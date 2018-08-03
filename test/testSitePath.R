tree <- ggtree::read.beast(
  system.file("m.trees", package = "sitePath")
)@phylo

align <- phangorn::read.phyDat(
  system.file("m.aligned.fasta", package = "sitePath"),
  format = "fasta", type = "DNA"
)

tree <- read.tree(
  system.file("test.tree", package = "sitePath")
)
align <- seqinr::read.alignment(
  system.file("test.aligned.fasta", package = "sitePath"),
  format = "fasta"
)
align <- as.phyDat(align, "AA")

matched <- treeAlignMatch(tree, align)
grouping <- groupTips(matched, 0.9)
mutations <- ancestralMutations(matched, 0.9)
# plot(mutations)

ggtree(groupClade(tree, node = 226), aes(linetype = group))
