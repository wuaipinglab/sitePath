context("test-typeConversion")

test_that("as.phylo", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tree <- addMSA(tree = zikv_tree_reduced,
                   alignment = zikv_align_reduced)
    paths <- lineagePath(tree)
    mutations <- fixationSites(paths)
    if (ape::is.binary(zikv_tree_reduced)) {
        expect_equal(zikv_tree_reduced, as.phylo(tree))
        expect_equal(zikv_tree_reduced, as.phylo(paths))
        expect_equal(zikv_tree_reduced, as.phylo(mutations))
    } else {
        resolved <- ape::multi2di(zikv_tree_reduced, random = FALSE)
        expect_equal(resolved, as.phylo(tree))
        expect_equal(resolved, as.phylo(paths))
        expect_equal(resolved, as.phylo(mutations))
    }
})
