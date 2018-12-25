context("Treemer")

test_that("Topology-dependent trimming", {
  tree <- ape::read.tree(
    system.file("ZIKV.newick", package = "sitePath")
  )
  tree <- ape::root(tree, "ANK57896")
  align <- seqinr::read.alignment(
    system.file("ZIKV.fasta", package = "sitePath"),
    format = "fasta"
  )
  simMatrix <- similarityMatrix(tree, align)
  minSim <- min(simMatrix)
  step <- round(minSim - 1, 3) / 50
  for (s in seq(1, minSim, step)) {
    grouping <- groupTips(tree, align, s, simMatrix, forbidTrivial = FALSE)
    expect_equal(
      sort(unlist(grouping, use.names = FALSE)),
      sort(tree$tip.label)
    )
    for (g in names(grouping)) {
      an <- as.integer(g)
      descendant <- grouping[[g]]
      if (length(descendant) == 1) {
        expect_equal(descendant, tree$tip.label[an])
      } else {
        expect_equal(an, ape::getMRCA(tree, descendant))
      }
      expect_equal(
        sort(descendant), 
        sort(ggtree::get.offspring.tip(tree, an))
      )
    }
    paths <- sitePath(tree, align, s, simMatrix, forbidTrivial = FALSE)
    for (p in paths) {
      expect_equal(ape::nodepath(tree, from = p[1], to = p[length(p)]), p)
    }
    if (length(paths) == 0) {
      break
    }
  }
})
