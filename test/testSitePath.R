# tree <- read.tree(
#   system.file("test.tree", package = "sitePath")
# )
# align <- seqinr::read.alignment(
#   system.file("test.aligned.fasta", package = "sitePath"),
#   format = "fasta"
# )
# align <- as.phyDat(align, "AA")

tree <- ggtree::read.beast(
  system.file("m.trees", package = "sitePath")
)@phylo

align <- phangorn::read.phyDat(
  system.file("m.aligned.fasta", package = "sitePath"),
  format = "fasta", type = "DNA"
)

matched <- treeAlignMatch(tree, align)
grouping <- groupTips(matched, 0.98)
divPoints <- getSitePath(matched, 0.98)

library(ggplot2)
library(ggtree)

group <- NULL
ggtree(groupClade(tree, node = 148), aes(linetype = group))

isTip <- NULL
node <- NULL
ggtree(matched$tree) +
  geom_text2(aes(subset = !isTip, label = node), hjust = -.3)

ggtree(matched$tree) +
  geom_text2(aes(subset = 185, label = node), hjust = -.3)
