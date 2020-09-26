test_that("The function works", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tr <- addMSA(tree = zikv_tree_reduced,
                 alignment = zikv_align_reduced)
    snp <- SNPsites(tr)
    expect_error(as.data.frame(snp), NA)
    p <- lineagePath(tr)
    fixedSites <- fixationSites(p)
    expect_error(as.data.frame(fixedSites), NA)
})
