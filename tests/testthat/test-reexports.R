data(zikv_align_reduced)
data(zikv_tree_reduced)
tr <- addMSA(tree = zikv_tree_reduced,
               alignment = zikv_align_reduced)
p <- lineagePath(tr)
fixedSites <- fixationSites(p)

test_that("The reexported 'as.phylo' function works", {
    if (ape::is.binary(zikv_tree_reduced)) {
        expect_equal(zikv_tree_reduced, as.phylo(tr))
        expect_equal(zikv_tree_reduced, as.phylo(p))
        expect_equal(zikv_tree_reduced, as.phylo(fixedSites))
    } else {
        resolved <- ape::multi2di(zikv_tree_reduced, random = FALSE)
        expect_equal(resolved, as.phylo(tr))
        expect_equal(resolved, as.phylo(p))
        expect_equal(resolved, as.phylo(fixedSites))
    }
})

test_that("The reexported 'as.treedata' function works", {
    expect_error(as.treedata(fixedSites), NA)
    fp <- fixationPath(fixedSites)
    expect_error(as.treedata(fp), NA)
})
