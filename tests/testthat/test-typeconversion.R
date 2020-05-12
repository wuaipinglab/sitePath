context("test-typeconversion")

test_that("Conversion to phylo works", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tree <-
        addMSA(h3n2_tree_reduced, alignment = h3n2_align_reduced)
    nTips <- length(tree$tip.label)
    paths <- lineagePath(tree)
    mutations <- fixationSites(paths)
    x <- as.phylo(mutations)
    checkOutput <- capture.output(ape::checkValidPhylo(x))
    expect_false(any(grepl("FATAL", checkOutput)))
})
