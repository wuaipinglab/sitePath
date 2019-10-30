context("test-SNPsites")

test_that("Restrication applied in SNPsites", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tree <-
        addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)

    nTips <- length(tree$tip.label)
    minT <- floor(nTips / 10)
    maxT <- ceiling(nTips / 2)
    for (n in seq(minT, maxT, 20)) {
        for (snp in SNPsites(tree, minSNP = n)) {
            aa <- table(sapply(zikv_align_reduced$seq, substring, snp, snp))
            expect_gte(sum(aa > n), 2)
        }
    }
})

context("test-fixationSites")

test_that("Constrains in fixationSites work", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tree <-
        addMSA(h3n2_tree_reduced, alignment = h3n2_align_reduced)
    nTips <- length(tree$tip.label)
    paths <- lineagePath(tree)
    mutations <- fixationSites(paths)
    minEffectiveSize <-
        length(tree$tip.label) / length(unique(unlist(paths)))
    for (sp in mutations) {
        site <- attr(sp, "site")
        for (mp in sp) {
            expect_lte(sum(lengths(mp)), nTips)
            aa <- lapply(mp, function(tips) {
                matchIndex <-
                    which(h3n2_align_reduced$nam %in% tree$tip.label[tips])
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
