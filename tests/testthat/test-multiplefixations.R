context("test-multiplefixations")

data(h3n2_tree_reduced)
data(h3n2_align_reduced)

test_that("Constrains in multiFixationSites work", {
    tree <- addMSA(h3n2_tree_reduced, alignment = h3n2_align_reduced)
    paths <- lineagePath(tree)
    mutations <- multiFixationSites(paths)
    minEffectiveSize <-
        length(tree$tip.label) / length(unique(unlist(paths)))
    for (sp in mutations) {
        site <- attr(sp, "site")
        for (m in sp) {
            aa <- lapply(m, function(g) {
                sum <-
                    h3n2_align_reduced$seq[which(h3n2_align_reduced$nam %in% tree$tip.label[g])]
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
