context("test-multiplefixations")

test_that("Constrains in multiFixationSites work", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tree <-
        addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
    nTips <- length(tree$tip.label)
    paths <- lineagePath(tree, similarity = 0.1)
    expect_warning(object = mutations <-
                       multiFixationSites(paths, samplingTimes = 10))
    for (sp in mutations) {
        site <- attr(sp, "site")
        for (mp in sp) {
            expect_lte(sum(lengths(mp)), nTips)
            for (i in seq_along(mp)) {
                nodeTips <- mp[[i]]
                expect_equal(sort(nodeTips), sort(unique(nodeTips)))
            }
        }
    }
})
