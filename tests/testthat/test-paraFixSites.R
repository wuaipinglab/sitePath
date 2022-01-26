test_that("The function works", {
    data(zikv_tree)
    data(zikv_align)

    zikv_paths <- addMSA(zikv_tree, alignment = zikv_align)
    zikv_entropy <- sitesMinEntropy(zikv_paths)

    expect_identical(paraFixSites(zikv_tree, alignment = zikv_align),
                     paraFixSites(zikv_paths))
    expect_identical(paraFixSites(zikv_paths), paraFixSites(zikv_entropy))
    expect_error(paraFixSites(zikv_entropy, category = "union"), NA)
    expect_error(paraFixSites(zikv_entropy, category = "parallelOnly"), NA)
    expect_error(paraFixSites(zikv_entropy, category = "fixationOnly"), NA)
})
