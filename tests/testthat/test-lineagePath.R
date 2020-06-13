context("test-lineagPath")

test_that("SNP-dependent trimming", {
    data(zikv_align)
    data(zikv_tree)
    tree <- addMSA(zikv_tree, alignment = zikv_align)
    nTips <- length(tree[["tip.label"]])
    align <- attr(tree, "align")
    reference <- attr(tree, "reference")
    similarity <- nTips * 0.1
    paths <- lineagePath(
        tree = tree,
        similarity = similarity,
        forbidTrivial = FALSE
    )
    # Exclude the invariant sites
    loci <- which(vapply(
        X = seq_along(reference),
        FUN = function(s) {
            s <- reference[s]
            length(unique(substr(align, s, s))) > 1
        },
        FUN.VALUE = logical(1)
    ))
    majorSNPsites <- list()
    for (site in loci) {
        for (p in paths) {
            seqs <- align[sitePath:::.childrenTips(tree, p[length(p)])]
            siteSummary <- sitePath:::tableAA(seqs, site - 1)
            siteChar <- names(siteSummary)[which(siteSummary == length(seqs))]
            if (length(siteSummary) != 0) {
                s <- as.character(site)
                if (s %in% names(majorSNPsites)) {
                    majorSNPsites[[s]] <- c(majorSNPsites[[s]], siteChar)
                } else {
                    majorSNPsites[[s]] <- siteChar
                }
            }
        }
    }
    for (site in names(majorSNPsites)) {
        siteSummary <- sitePath:::tableAA(align, as.integer(site) - 1)
        expect_true(any(siteSummary > similarity))
    }
})
