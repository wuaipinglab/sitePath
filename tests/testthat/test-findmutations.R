data(zikv_align)
data(zikv_tree)
tree <- addMSA(zikv_tree, alignment = zikv_align)

context("test-SNPsites")

test_that("Restrication applied in SNPsites", {
    nTips <- length(tree$tip.label)
    minT <- floor(nTips / 10)
    maxT <- ceiling(nTips / 2)
    for (n in seq(minT, maxT, 10)) {
        for (snp in SNPsites(tree, minSNP = n)) {
            aa <- table(sapply(zikv_align$seq, substring, snp, snp))
            expect_gte(sum(aa > n), 2)
        }
    }
})

context("test-fixationSites")

test_that("Constrains in fixationSites work", {
    paths <- lineagePath(tree)
    for (ta in seq(0, 0.49, 0.2)) {
        for (td in seq(0, 0.49, 0.2)) {
            mutations <- fixationSites(paths, tolerance = c(ta, td))
            for (sp in mutations) {
                site <- attr(sp, "site")
                for (m in sp) {
                    aa <- lapply(m, function(g) {
                        sum <- 
                            zikv_align$seq[which(zikv_align$nam %in% tree$tip.label[g])]
                        sum <-
                            table(sapply(sum, substring, site, site))
                    })
                    ancAA <- aa[[1]]
                    ancD <- which.max(ancAA)
                    expect_equal(toupper(names(ancD)), names(m)[1])
                    expect_lte(sum(ancAA[-ancD]), floor(ta * sum(ancAA)))
                    descAA <- aa[[2]]
                    descD <- which.max(descAA)
                    expect_equal(toupper(names(descD)), names(m)[2])
                    expect_lte(sum(descAA[-descD]), floor(td * sum(descAA)))
                }
            }
        }
    }
})
