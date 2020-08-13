context("test-groupTips")

test_that("The grouped tips include all tree tips", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tree <- addMSA(tree = zikv_tree_reduced,
                   alignment = zikv_align_reduced)
    tipNames <- as.phylo(tree)[["tip.label"]]
    for (s in seq(0.1, 0.05, length.out = 5)) {
        # The grouping from the divergent nodes
        paths <- lineagePath(tree, similarity = s)
        grouped <- groupTips(paths)
        allTips <- unlist(grouped)
        names(allTips) <- NULL
        expect_identical(sort(allTips), sort(tipNames))
        # The grouping from the fixation mutation
        minEntropy <- sitesMinEntropy(paths)
        grouped <- groupTips(minEntropy)
        allTips <- unlist(grouped)
        names(allTips) <- NULL
        expect_identical(sort(allTips), sort(tipNames))
        # The grouping from the SNP tracing
        snpTracing <- fixationPath(minEntropy, minEffectiveSize = 0)
        grouped <- groupTips(snpTracing)
        allTips <- unlist(grouped)
        names(allTips) <- NULL
        expect_identical(sort(allTips), sort(tipNames))
    }
})
