context("test-parallelSites")

test_minEffectiveSize <- function(tips, minSNP) {
    tipNum <- length(unlist(tips))
    isFixed <- vapply(
        X = tips,
        FUN = attr,
        FUN.VALUE = logical(1),
        which = "fixed"
    )
    # This check could be improved because the 'minSNP' constrain
    # applies on the two lineage separately
    expect_true(tipNum >= minSNP * 2  ||
                    sum(isFixed) >= 2 ||
                    tipNum >= minSNP &&
                    any(isFixed))
}

test_that("Constrains in parallelSites works", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tree <- addMSA(tree = h3n2_tree_reduced,
                   alignment = h3n2_align_reduced)
    paths <- lineagePath(tree)
    x <- sitesMinEntropy(paths)
    for (s in c(1, 5, 15)) {
        mutations <- parallelSites(x, minSNP = s, mutMode = "all")
        for (sp in mutations) {
            mutTips <- extractTips(sp)
            test_minEffectiveSize(mutTips, s)
        }
        mutations <-
            parallelSites(x, minSNP = s, mutMode = "exact")
        for (sp in mutations) {
            mutTips <- extractTips(sp)
            mutModeTips <- list()
            for (tips in mutTips) {
                mutName <- attr(tips, "mutName")[4]
                mutModeTips[[mutName]] <- c(mutModeTips[[mutName]],
                                            list(tips))
            }
            for (tips in mutModeTips) {
                test_minEffectiveSize(tips, s)
            }
        }
        mutations <-
            parallelSites(x, minSNP = s, mutMode = "pre")
        for (sp in mutations) {
            mutTips <- extractTips(sp)
            for (tips in mutTips) {
                mutName <- attr(tips, "mutName")[1]
                mutModeTips[[mutName]] <- c(mutModeTips[[mutName]],
                                            list(tips))
            }
            for (tips in mutModeTips) {
                test_minEffectiveSize(tips, s)
            }
        }
        mutations <-
            parallelSites(x, minSNP = s, mutMode = "post")
        for (sp in mutations) {
            mutTips <- extractTips(sp)
            for (tips in mutTips) {
                mutName <- attr(tips, "mutName")[3]
                mutModeTips[[mutName]] <- c(mutModeTips[[mutName]],
                                            list(tips))
            }
            for (tips in mutModeTips) {
                test_minEffectiveSize(tips, s)
            }
        }
    }
})
