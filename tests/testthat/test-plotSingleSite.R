context("test-plotSingleSite")

test_that("The function plotSingleSite works for AA", {
    data(zikv_align)
    data(zikv_tree)
    tr <- addMSA(tree = zikv_tree,
                 alignment = zikv_align,
                 seqType = "AA")
    p <- lineagePath(tr)
    expect_error(plotSingleSite(p, 2, showPath = TRUE), NA)
    expect_error(plotSingleSite(p, 3, showPath = FALSE), NA)
    muts <- fixationSites(p)
    expect_error(plotSingleSite(muts, 139), NA)
})

test_that("The function plotSingleSite works for DNA", {
    data(sars2_align)
    data(sars2_tree)
    tr <- addMSA(tree = sars2_tree,
                 alignment = sars2_align,
                 seqType = "DNA")
    p <- lineagePath(tr)
    expect_error(plotSingleSite(p, 2, showPath = TRUE), NA)
    expect_error(plotSingleSite(p, 3, showPath = FALSE), NA)
    muts <- fixationSites(p)
    expect_error(plotSingleSite(muts, 8517), NA)
})

test_that("Warnings in plotSingleSite works", {
    data(zikv_align)
    data(zikv_tree)
    tr <- addMSA(zikv_tree, alignment = zikv_align)
    p <- lineagePath(tr)
    expect_error(plotSingleSite(p, 0))
    muts <- fixationSites(p)
    expect_error(extractSite(muts, 0))
    expect_error(plotSingleSite(muts, 0))
    sp <- extractSite(muts, 139)
    nSP <- length(sp)
    expect_error(extractTips(sp, nSP + 1))
    expect_error(plot(sp, nSP + 1))
})
