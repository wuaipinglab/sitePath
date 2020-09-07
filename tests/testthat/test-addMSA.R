test_that("Read from file works", {
    tree_file <- system.file("extdata",
                             "ZIKV.newick",
                             package = "sitePath")
    alignment_file <- system.file("extdata",
                                  "ZIKV.fasta",
                                  package = "sitePath")
    tr <- read.tree(tree_file)
    expect_error(addMSA(tr, alignment_file), NA)
    expect_error(addMSA(tr, msaPath = ""))

    data(zikv_align)
    expect_identical(addMSA(tr, alignment_file),
                     addMSA(tr, alignment_file, alignment = zikv_align))
})

test_that("Read from object works", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    expect_error(addMSA(tree = zikv_tree_reduced,
                        alignment = zikv_align_reduced),
                 NA)

})
