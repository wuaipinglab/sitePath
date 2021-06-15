test_that("The plot functions work", {
    data(h3n2_align_reduced)
    data(h3n2_tree_reduced)

    paths <- addMSA(h3n2_tree_reduced,
                    alignment = h3n2_align_reduced)
    expect_error(plot(paths), NA)
    expect_error(plotMutSites(paths), NA)

    minEntropy <- sitesMinEntropy(paths)

    fixedSites <- fixationSites(minEntropy)
    expect_error(plot(fixedSites), NA)
    expect_error(plotMutSites(fixedSites), NA)

    # Test the constrain for 'sitePath'
    for (site in allSitesName(fixedSites)) {
        sp <- extractSite(fixedSites, site)
        expect_error(plot(sp, select = length(sp) + 1))
    }

    paraSites <- parallelSites(minEntropy)
    expect_error(plot(paraSites), NA)
    expect_error(plotMutSites(paraSites), NA)

    categories <- c("intersect", "union", "parallelOnly", "fixationOnly")
    for (mutMode in categories) {
        pf <- paraFixSites(minEntropy)
        expect_error(plotMutSites(pf), NA)
    }

    fp <- fixationPath(fixedSites)
    expect_error(plot(fp), NA)
})
