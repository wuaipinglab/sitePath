source("treemer.R")

tree <- read.tree(
  "/home/chengyang/PycharmProjects/FLUSITE/tests/dummy/RAxML_bestTree.mixed"
)
align <- read.alignment(
  "/home/chengyang/PycharmProjects/FLUSITE/tests/dummy/Alignment.fasta",
  "fasta"
)
grouping <- groupTips(tree, align)
mutations <- ancestralMutations(tree, align, "AA", "JTT", 1)
mutations2graphviz(mutations)
