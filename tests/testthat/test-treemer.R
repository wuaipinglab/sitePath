context("test-treemer")

test_that("Topology-dependent trimming", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
    nTips <- length(tree[["tip.label"]])
    simMatrix <- similarityMatrix(tree)
    expect_identical(colnames(simMatrix), tree$tip.label)
    expect_identical(row.names(simMatrix), tree$tip.label)
    minSim <- min(simMatrix)
    for (s in seq(1, minSim, length.out = 10)) {
        grouping <- groupTips(
            tree = tree,
            similarity = s,
            simMatrix = simMatrix,
            forbidTrivial = FALSE,
            tipnames = FALSE
        )
        expect_equal(sort(unlist(grouping, use.names = FALSE)), 1:nTips)
        for (g in names(grouping)) {
            an <- as.integer(g)
            index <- which(row.names(simMatrix) %in% grouping[[g]])
            expect_true(all(simMatrix[index, index] > s))
            descendant <- grouping[[g]]
            if (length(descendant) == 1) {
                expect_equal(descendant, an)
            } else {
                expect_equal(an, ape::getMRCA(tree, descendant))
            }
            expect_equal(sort(descendant),
                         sort(sitePath:::.childrenTips(tree, an)))
        }
    }
})
