context("Treemer")

test_that("Topology-dependent trimming", {
    data("zikv_align")
    data("zikv_tree")
    tree <- addMSA(zikv_tree, seqs = zikv_align)
    simMatrix <- similarityMatrix(tree)
    minSim <- min(simMatrix)
    step <- round(minSim - 1, 3) / 50
    for (s in seq(1, minSim, step)) {
        grouping <-
            groupTips(
                tree,
                similarity = s,
                simMatrix = simMatrix,
                forbidTrivial = FALSE,
                tipnames = FALSE
            )
        expect_equal(sort(unlist(grouping, use.names = FALSE)),
                     1:length(tree$tip.label))
        for (g in names(grouping)) {
            an <- as.integer(g)
            descendant <- grouping[[g]]
            if (length(descendant) == 1) {
                expect_equal(descendant, an)
            } else {
                expect_equal(an, ape::getMRCA(tree, descendant))
            }
            expect_equal(sort(descendant),
                         sort(sitePath:::ChildrenTips(tree, an)))
        }
        paths <-
            sitePath(
                tree,
                similarity = s,
                simMatrix = simMatrix,
                forbidTrivial = FALSE
            )
        for (p in paths) {
            expect_equal(ape::nodepath(tree, from = p[1], to = p[length(p)]), p)
        }
        if (length(paths) == 0) {
            break
        }
    }
})
