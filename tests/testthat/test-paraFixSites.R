test_that("The function works", {
    data(zikv_tree)
    data(zikv_align)

    zikv_paths <- addMSA(zikv_tree, alignment = zikv_align)
    zikv_entropy <- sitesMinEntropy(zikv_paths)

    expect_error(paraFixSites(zikv_entropy), NA)
    expect_error(paraFixSites(zikv_entropy, category = "union"), NA)
    expect_error(paraFixSites(zikv_entropy, category = "parallelOnly"), NA)
    expect_error(paraFixSites(zikv_entropy, category = "fixationOnly"), NA)
})
