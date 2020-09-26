test_that("The return value contains correct extra info", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tr <- addMSA(tree = h3n2_tree_reduced,
                 alignment = h3n2_align_reduced)
    p <- lineagePath(tr)
    minEntropy <- sitesMinEntropy(p)
    # The result is grouped by path
    expect_true(length(p) == length(minEntropy))
    for (segs in minEntropy) {
        for (seg in segs) {
            expect_true(all(sapply(names(seg), function(node) {
                tips <- seg[[node]]
                # The major amino acid/nucleotide should be same as the fixed
                # and the node names should be the same
                dominantAA <-
                    names(which.max(attr(tips, "aaSummary")))
                dominantAA == attr(tips, "AA") &&
                node == attr(tips, "node")
            })))
        }
    }
})
