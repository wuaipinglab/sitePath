data(zikv_align)
data(zikv_tree)
tree <- addMSA(zikv_tree, seqs = zikv_align)

context("test-SNPsites")

test_that("Restrication applied in SNPsites", {
    nTips <- length(tree$tip.label)
    minT <- floor(nTips / 10)
    maxT <- ceiling(nTips / 2)
    for (n in seq(minT, maxT, 5)) {
        for (snp in SNPsites(tree, minSNP = n)) {
            aa <- table(sapply(zikv_align$seq, substring, snp, snp))
            expect_gte(sum(aa > n), 2)
        }
    }
})

context("test-fixationSites")

test_that("Restrication applied in fixationSites", {
    paths <- sitePath(tree)
    for (ta in seq(0, 0.49, 0.05)) {
        for (td in seq(0, 0.49, 0.05)) {
            mutations <- fixationSites(paths, tolerance = c(ta, td))
            for (m in names(mutations)) {
                nLen <- nchar(m)
                site <-
                    substring(m, c(1, 2, nLen), c(1, nLen - 1, nLen))
                aa <- lapply(mutations[[m]], function(g) {
                    sum <- zikv_align$seq[which(zikv_align$nam %in% g)]
                    sum <-
                        table(sapply(sum, substring, site[2], site[2]))
                })
                ancAA <- aa[[1]]
                ancD <- which.max(ancAA)
                expect_equal(toupper(names(ancD)), site[1])
                expect_lte(sum(ancAA[-ancD]), floor(ta * sum(ancAA)))
                descAA <- aa[[2]]
                descD <- which.max(descAA)
                expect_equal(toupper(names(descD)), site[3])
                expect_lte(sum(descAA[-descD]), floor(td * sum(descAA)))
            }
        }
    }
})
