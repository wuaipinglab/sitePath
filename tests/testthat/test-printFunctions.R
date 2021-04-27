test_that("The print functions work", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)

    tr <- addMSA(zikv_tree_reduced,
                 alignment = zikv_align_reduced)
    expect_error(capture.output(print(tr)), NA)

    snp <- SNPsites(tr)
    expect_error(capture.output(print(snp)), NA)

    p <- lineagePath(tr)
    expect_error(capture.output(print(p)), NA)

    snp2 <- SNPsites(p)
    expect_error(capture.output(print(snp2)), NA)

    minEntropy <- sitesMinEntropy(p)
    expect_error(capture.output(print(minEntropy)), NA)

    fixedSites <- fixationSites(minEntropy)
    expect_error(capture.output(print(fixedSites)), NA)

    sp <- extractSite(fixedSites, 139)
    expect_error(capture.output(print(sp)), NA)

    fp <- fixationPath(fixedSites)
    expect_error(capture.output(print(fp)), NA)

    paraSites <- parallelSites(minEntropy, minSNP = 1)
    expect_error(capture.output(print(paraSites)), NA)

    para <- extractSite(paraSites, 106)
    expect_error(capture.output(print(para)), NA)
})
