test_that("The output is valid phylo", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tree <- addMSA(tree = h3n2_tree_reduced,
                   alignment = h3n2_align_reduced)
    nTips <- length(as.phylo(tree)$tip.label)
    paths <- lineagePath(tree)
    mutations <- fixationSites(paths)
    expect_false(any(duplicated(unlist(
        attr(mutations, "clustersByPath")
    ))))
    x <- fixationPath(mutations)
    tr <- attr(x, "SNPtracing")@phylo
    checkOutput <- capture.output(ape::checkValidPhylo(tr))
    expect_false(any(grepl("FATAL", checkOutput)))
})
