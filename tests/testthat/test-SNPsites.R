context("test-SNPsites")

test_that("Restrication applied in SNPsites", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tree <- addMSA(tree = zikv_tree_reduced,
                   alignment = zikv_align_reduced)

    nTips <- length(as.phylo(tree)$tip.label)
    minT <- floor(nTips / 10)
    maxT <- ceiling(nTips / 2)
    for (n in seq(minT, maxT, 20)) {
        for (snp in SNPsites(tree, minSNP = n)) {
            aa <- table(sapply(zikv_align_reduced$seq, substring, snp, snp))
            expect_gte(sum(aa > n), 2)
        }
    }
})
