context("test-fixationSites")

test_that("Constrains in fixationSites work", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tree <- addMSA(tree = h3n2_tree_reduced,
                   alignment = h3n2_align_reduced)
    tipNames <- as.phylo(tree)$tip.label
    nTips <- length(tipNames)
    paths <- lineagePath(tree)
    mutations <- fixationSites(paths)
    expect_false(any(duplicated(unlist(
        attr(mutations, "clustersByPath")
    ))))
    minEffectiveSize <- nTips / length(unique(unlist(paths)))
    for (sp in mutations) {
        site <- attr(sp, "site")
        for (mp in sp) {
            expect_lte(sum(lengths(mp)), nTips)
            aa <- lapply(mp, function(tips) {
                matchIndex <-
                    which(h3n2_align_reduced$nam %in% tipNames[tips])
                sum <- h3n2_align_reduced$seq[matchIndex]
                sum <- table(sapply(sum, substring, site, site))
            })
            for (i in seq_along(mp)) {
                nodeTips <- mp[[i]]
                expect_equal(sort(nodeTips), sort(unique(nodeTips)))
                expect_gte(length(nodeTips), minEffectiveSize)
                aaD <- toupper(names(which.max(aa[[i]])))
                fixedAA <- attr(nodeTips, "AA")
                attributes(fixedAA) <- NULL
                expect_equal(fixedAA, aaD)
            }
        }
    }
})
