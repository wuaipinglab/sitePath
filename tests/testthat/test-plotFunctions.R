test_that("The plot functions work", {
    data(zikv_align)
    data(zikv_tree)
    tr <- addMSA(zikv_tree,
                 alignment = zikv_align)

    p <- lineagePath(tr)
    expect_error(plot(p), NA)

    fixedSites <- fixationSites(p)
    expect_error(plot(fixedSites), NA)

    # Test the constrain for 'sitePath'
    sp <- extractSite(fixedSites, 139)
    nSP <- length(sp)
    expect_error(plot(sp, select = nSP + 1))

    fp <- fixationPath(fixedSites)
    expect_error(plot(fp), NA)
})
