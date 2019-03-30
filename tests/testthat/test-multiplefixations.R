context("test-multiplefixations")

data(h3n2_tree)
data(h3n2_align)

h3n2_tree <- ape::drop.tip(h3n2_tree, 1:1000)
h3n2_align$seq <- lapply(h3n2_align$seq, FUN = function(s) {
    substring(s, 20, 550)
})

test_that("Constrains in multiFixationSites work", {
    tree <- addMSA(h3n2_tree, alignment = h3n2_align)
    paths <- lineagePath(tree)
    mutations <- multiFixationSites(paths)
    minEffectiveSize <- length(tree$tip.label) / (length(paths) * 10)
    for (sp in mutations) {
        site <- attr(sp, "site")
        for (m in sp) {
            aa <- lapply(m, function(g) {
                sum <- 
                    h3n2_align$seq[which(h3n2_align$nam %in% tree$tip.label[g])]
                sum <-
                    table(sapply(sum, substring, site, site))
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
