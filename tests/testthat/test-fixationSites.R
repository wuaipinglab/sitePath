context("test-fixationSites")

test_that("Constrains in fixationSites work", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tree <- addMSA(tree = h3n2_tree_reduced,
                   alignment = h3n2_align_reduced)
    tipNames <- as.phylo(tree)$tip.label
    nTips <- length(tipNames)
    p <- lineagePath(tree)
    fixSites <- fixationSites(p)
    minEntropy <- sitesMinEntropy(p)
    # Test the two function give the same result
    expect_identical(fixSites, fixationSites(minEntropy))
    expect_false(any(duplicated(unlist(
        attr(fixSites, "clustersByPath")
    ))))
    # Test the number of tips before and after each fixation mutation is enough
    minEffectiveSize <- nTips / length(unique(unlist(p)))
    for (sp in fixSites) {
        site <- attr(sp, "site")
        for (mp in sp) {
            # The number of tips should be less than total tree tips
            expect_lte(sum(lengths(mp)), nTips)
            # AA/NT summary of the tips on the 'mutPath'
            aa <- lapply(mp, function(tips) {
                matchIndex <-
                    which(h3n2_align_reduced[["nam"]] %in% tipNames[tips])
                sum <- h3n2_align_reduced[["seq"]][matchIndex]
                sum <- table(sapply(sum, substring, site, site))
            })
            for (i in seq_along(mp)) {
                # The tips in each group before or after fixation mutation
                nodeTips <- mp[[i]]
                # Tips should be unique and more than 'minEffectiveSize'
                expect_equal(sort(nodeTips), sort(unique(nodeTips)))
                expect_gte(length(nodeTips), minEffectiveSize)
                # The dominant AA/NT should be the fixed one
                aaD <- toupper(names(which.max(aa[[i]])))
                fixedAA <- attr(nodeTips, "AA")
                attributes(fixedAA) <- NULL
                expect_equal(fixedAA, aaD)
            }
        }
    }
})
