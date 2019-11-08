library(sitePath)

test_that("setSiteNumbering works", {
    data(zikv_align)
    data(zikv_tree)
    tree <- addMSA(zikv_tree, alignment = zikv_align)
    tipNames <- zikv_tree[["tip.label"]]
    align <- attr(tree, "align")
    for (i in seq_along(tipNames)) {
        tip <- tipNames[[i]]
        refSeq <- align[[i]]
        expect_identical(refSeq,
                         toupper(zikv_align[["seq"]][[which(zikv_align[["nam"]] == tip)]]))
        paths <- lineagePath(tree)
        paths <- setSiteNumbering(paths, tip, "C")
        reference <- attr(paths, "reference")
        refSeqName <- attr(reference, "refSeqName")
        expect_equal(tip, refSeqName)
        ungappedRefSeq <- gsub("C", "", refSeq)
        aaEqual <- all(sapply(seq_along(reference), function(n) {
            # Corresponding site index for the original sequence
            m <- reference[[n]]
            substr(ungappedRefSeq, n, n) == substr(refSeq, m, m)
        }))
        expect_true(aaEqual)

        tree <- setSiteNumbering(tree, tip, "G")
        reference <- attr(tree, "reference")
        refSeqName <- attr(reference, "refSeqName")
        expect_equal(tip, refSeqName)
        ungappedRefSeq <- gsub("G", "", refSeq)
        aaEqual <- all(sapply(seq_along(reference), function(n) {
            # Corresponding site index for the original sequence
            m <- reference[[n]]
            substr(ungappedRefSeq, n, n) == substr(refSeq, m, m)
        }))
        expect_true(aaEqual)
    }
})
