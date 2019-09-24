context("test-multiplefixations")

test_that("Constrains in multiFixationSites work", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tree <-
        addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
    paths <- lineagePath(tree)
    expect_warning(object = mutations <-
                       multiFixationSites(paths, samplingTimes = 10))
    for (sp in mutations) {
        site <- attr(sp, "site")
        for (mp in sp) {
            aa <- lapply(mp, function(tips) {
                sum <- zikv_align_reduced$seq[which(is.element(
                    el = zikv_align_reduced$nam,
                    set = tree$tip.label[tips]
                ))]
                sum <- table(sapply(sum, substring, site, site))
            })
            for (i in seq_along(mp)) {
                nodeTips <- mp[[i]]
                aaD <- toupper(names(which.max(aa[[i]])))
                tipsAA <- attr(nodeTips, "tipsAA")
                tipsAA <- sitePath:::tableAA(tipsAA, 0)
                aaD2 <- toupper(names(which.max(aa[[i]])))
                expect_equal(aaD, aaD2)
            }
        }
    }
})
