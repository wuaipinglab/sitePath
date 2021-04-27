test_that("The function works", {
    data(zikv_tree)
    data(zikv_align)

    zikv_tr <- addMSA(zikv_tree, alignment = zikv_align)

    zikv_snp <- SNPsites(zikv_tr)

    sites <- allSitesName(zikv_snp)
    expect_type(sites, "character")
    expect_equal(length(sites), length(zikv_snp))

    zikv_p <- lineagePath(zikv_tr)

    zikv_entropy <- sitesMinEntropy(zikv_p)

    zikv_fixed <- fixationSites(zikv_entropy)
    sites <- allSitesName(zikv_fixed)
    expect_type(sites, "character")
    expect_equal(length(sites), length(zikv_fixed))

    zikv_para <- parallelSites(zikv_entropy, minSNP = 1)
    sites <- allSitesName(zikv_para)
    expect_type(sites, "character")
    expect_equal(length(sites), length(zikv_para))
})
