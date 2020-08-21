test_minSNPandSiteSkipping <- function(minEntropy, unambiguous) {
    # Use different 'minSNP' constrains and 'mutMode'
    mutModes <- c("pre", "all", "post", "exact")
    for (minSNP in c(1, 5, 15)) {
        for (i in seq_along(mutModes)) {
            mutations <- parallelSites(minEntropy,
                                       minSNP = minSNP,
                                       mutMode = mutModes[i])
            for (site in names(mutations)) {
                sp <- mutations[[site]]
                mutTips <- extractTips(sp)
                # Group the tips according to the 'mutMode'
                mutModeTips <- list()
                for (tips in mutTips) {
                    mutName <- attr(tips, "mutName")[i]
                    mutModeTips[[mutName]] <-
                        c(mutModeTips[[mutName]],
                          list(tips))
                }
                for (tips in mutModeTips) {
                    mutChars <- unlist(lapply(tips, function(i) {
                        attr(i, "mutName")[c(1, 3)]
                    }), use.names = FALSE)
                    # The gap and ambiguous character should have been ignored
                    expect_true(all(mutChars %in% unambiguous),
                                label = paste(site, mutModes[i], minSNP))
                    tipNum <- length(unlist(tips))
                    isFixed <- vapply(
                        X = tips,
                        FUN = attr,
                        FUN.VALUE = logical(1),
                        which = "fixed"
                    )
                    # This check could be improved because the 'minSNP'
                    # constrain applies on the two lineage separately. Here it
                    # was tested for all the result of a site
                    expect_true(
                        tipNum >= minSNP * 2  ||
                            sum(isFixed) >= 2 ||
                            tipNum >= minSNP &&
                            any(isFixed)
                    )
                }
            }
        }
    }
}

test_that("The function works for amino acid", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tr <- addMSA(tree = h3n2_tree_reduced,
                 alignment = h3n2_align_reduced)
    gapChar <- "S"
    p <- lineagePath(tr)
    p <- setSiteNumbering(p, gapChar = gapChar)
    # Test the input of 'minSNP' and 'mutMode'
    expect_error(parallelSites(p, minSNP = "a"))
    expect_error(parallelSites(p, mutMode = "asdfa"))
    minEntropy <- sitesMinEntropy(p)
    unambiguous <- setdiff(sitePath:::AA_UNAMBIGUOUS, gapChar)
    test_minSNPandSiteSkipping(minEntropy, unambiguous)
})

test_that("The function works for nucleotide", {
    data(sars2_tree)
    data(sars2_align)
    tr <- addMSA(tree = sars2_tree,
                 alignment = sars2_align,
                 seqType = "DNA")
    gapChar <- "G"
    p <- lineagePath(tr)
    p <- setSiteNumbering(p, gapChar = gapChar)
    # Test the input of 'minSNP' and 'mutMode'
    expect_error(parallelSites(p, minSNP = "a"))
    expect_error(parallelSites(p, mutMode = "asdfa"))
    minEntropy <- sitesMinEntropy(p)
    unambiguous <- setdiff(sitePath:::NT_UNAMBIGUOUS, gapChar)
    test_minSNPandSiteSkipping(minEntropy, unambiguous)
})
