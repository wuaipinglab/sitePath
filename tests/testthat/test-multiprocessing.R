test_that("Multiprocessing works", {
    data(h3n2_tree)
    data(h3n2_align)

    tr <- addMSA(tree = h3n2_tree, alignment = h3n2_align)
    p <- lineagePath(tr)
    minEntropy <- sitesMinEntropy(p)

    options(cl.cores = "auto")
    tr_mp <- addMSA(tree = h3n2_tree, alignment = h3n2_align)
    p_mp <- lineagePath(tr)
    minEntropy_mp <- sitesMinEntropy(p_mp)
    options(cl.cores = NULL)

    expect_identical(tr, tr_mp)
    expect_identical(minEntropy, minEntropy_mp)
})
