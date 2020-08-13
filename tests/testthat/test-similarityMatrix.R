context("test-similarityMatrix")

test_that("Calculate similarity matrix", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tree <- addMSA(tree = zikv_tree_reduced,
                   alignment = zikv_align_reduced)
    tipNames <- as.phylo(tree)[["tip.label"]]
    simMatrix <- similarityMatrix(tree)
    expect_identical(colnames(simMatrix), tipNames)
    expect_identical(row.names(simMatrix), tipNames)
})
