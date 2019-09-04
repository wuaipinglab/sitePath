data(zikv_align_reduced)
data(zikv_tree_reduced)
tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)

context("test-SNPsites")

test_that("Restrication applied in SNPsites", {
    nTips <- length(tree$tip.label)
    minT <- floor(nTips/10)
    maxT <- ceiling(nTips/2)
    for (n in seq(minT, maxT, 20)) {
        for (snp in SNPsites(tree, minSNP = n)) {
            aa <- table(sapply(zikv_align_reduced$seq, substring, snp, snp))
            expect_gte(sum(aa > n), 2)
        }
    }
})

data(h3n2_tree_reduced)
data(h3n2_align_reduced)

context("test-fixationSites")

test_that("Constrains in fixationSites work", {
    tree <- addMSA(h3n2_tree_reduced, alignment = h3n2_align_reduced)
    paths <- lineagePath(tree)
    mutations <- fixationSites(paths)
    minEffectiveSize <- length(tree$tip.label)/length(unique(unlist(paths)))
    for (sp in mutations) {
        site <- attr(sp, "site")
        for (m in sp) {
            aa <- lapply(m, function(g) {
                sum <- h3n2_align_reduced$seq[which(h3n2_align_reduced$nam %in% tree$tip.label[g])]
                sum <- table(sapply(sum, substring, site, site))
            })
            for (i in seq_along(m)) {
                nodeTips <- m[[i]]
                expect_gte(length(nodeTips), minEffectiveSize)
                aaD <- toupper(names(which.max(aa[[i]])))
                fixedAA <- attr(nodeTips, "AA")
                attributes(fixedAA) <- NULL
                expect_equal(fixedAA, aaD)
            }
        }
    }
})
