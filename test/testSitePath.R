# tree <- read.tree(
#   system.file("test.tree", package = "sitePath")
# )
# align <- seqinr::read.alignment(
#   system.file("test.aligned.fasta", package = "sitePath"),
#   format = "fasta"
# )
# align <- as.phyDat(align, "AA")

library(ggplot2)
library(ggtree)

tree <- ggtree::read.beast(
  system.file("m.trees", package = "sitePath")
)@phylo

align <- seqinr::read.alignment(
  system.file("m.aligned.fasta", package = "sitePath"),
  format = "fasta", forceToLower = FALSE
)

matched <- treeAlignMatch(tree, align)
grouping <- groupTips(matched, 0.98)
sitePath <- getSitePath(matched, 0.98)
findSites(sitePath, "isolated", 1)


group <- NULL
ggtree(groupClade(tree, node = 148), aes(linetype = group))

isTip <- NULL
node <- NULL
ggtree(tree) +
  geom_text2(aes(subset = !isTip, label = node), hjust = -.3)

ggtree(matched$tree) +
  geom_text2(aes(subset = 185, label = node), hjust = -.3)
