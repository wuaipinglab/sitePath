test_that("Using SNP to get phylogenetic lineages", {
    data(zikv_align)
    data(zikv_tree)
    tr <- addMSA(zikv_tree, alignment = zikv_align)
    expect_true(is(tr, "phyMSAmatched"))
    tipNames <- zikv_tree[["tip.label"]]
    nTips <- length(tipNames)
    # The tree is bifurcated and should be identical to the original tree
    expect_identical(as.phylo(tr), zikv_tree)
    align <- attr(tr, "align")
    reference <- attr(tr, "msaNumbering")
    # Test the input of 'similarity'
    expect_error(lineagePath(tree = tr, similarity = "0.96"))
    expect_error(lineagePath(tree = tr, similarity = -0.1))
    expect_error(lineagePath(
        tree = tr,
        similarity = 1,
        forbidTrivial = FALSE
    ))
    # Use 0.1 as the input for 'similarity'
    similarity <- nTips * 0.1
    paths <- lineagePath(tree = tr,
                         similarity = similarity,
                         forbidTrivial = FALSE)
    # Exclude the invariant sites
    loci <- attr(tr, "loci")
    majorSNPsites <- list()
    for (site in loci) {
        for (p in paths) {
            # The descendant tips of the terminal node of the path
            seqs <-
                align[sitePath:::.childrenTips(as.phylo(tr), p[length(p)])]
            siteSummary <- sitePath:::tableAA(seqs, site - 1)
            # Get the SNP possibly exclusive to the descendant tips
            siteChar <-
                names(siteSummary)[which(siteSummary == length(seqs))]
            # Add the SNP to 'majorSNPsites'
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
    # The SNP should meet the 'similarity' threshold
    for (site in names(majorSNPsites)) {
        siteSummary <- sitePath:::tableAA(align, as.integer(site) - 1)
        expect_true(any(siteSummary > similarity))
    }
})

test_that("The sneakPeek function works", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tr <- addMSA(zikv_tree_reduced,
                 alignment = zikv_align_reduced)
    rangeOfResults <- sneakPeek(tr)
    expect_error(sneakPeek(tr), NA)
    for (i in seq_len(nrow(rangeOfResults))) {
        similarity <- rangeOfResults[i, "similarity"]
        pathNum <- rangeOfResults[i, "pathNum"]
        p <- lineagePath(rangeOfResults, similarity = similarity)
        expect_true(is(p, "lineagePath"))
        expect_true(length(p) == pathNum,
                    label = paste("Simialrity:", similarity))
    }
})
