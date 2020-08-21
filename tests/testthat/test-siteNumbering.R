test_referenceNumbering <- function(x, tip, refSeq, gapChar) {
    # The 'gapChar' can only be a single character
    expect_error(setSiteNumbering(x, tip, "--"))
    expect_error(setSiteNumbering(x, tip, c("-", "-")))
    # Use an arbitrary character as the gap character
    x <- setSiteNumbering(x, tip, gapChar)
    reference <- attr(x, "msaNumbering")
    refSeqName <- attr(x, "reference")
    # The correct reference sequence name can be retrieved
    expect_equal(tip, refSeqName)
    # Get the numbering of the reference sequence
    refNumbering <- seq_along(reference)
    # Create an ungapped reference sequence and compared with the gapped
    # reference sequence using their own numbering schemes
    ungappedRefSeq <- gsub(gapChar, "", refSeq)
    # The amino acid/nucleotide are the same for each mapped site
    aaEqual <- sapply(refNumbering, function(n) {
        # Corresponding site index for the original sequence
        m <- reference[[n]]
        substr(ungappedRefSeq, n, n) == substr(refSeq, m, m)
    })
    expect_true(all(aaEqual))
}

test_siteSkipping <- function(x, tip, gapChar, minSkipSize) {
    # The skip size cannot be zero or negative
    expect_error(setSiteNumbering(x, tip, gapChar, 0))
    expect_error(setSiteNumbering(x, tip, gapChar, -2))
    # Use an arbitrary character as the gap character
    x <- setSiteNumbering(x, tip, gapChar, minSkipSize)
    # Get all the unambiguous character
    if (attr(x, "seqType") == "AA") {
        unambiguous <- setdiff(sitePath:::AA_UNAMBIGUOUS, gapChar)
    } else {
        unambiguous <- setdiff(sitePath:::NT_UNAMBIGUOUS, gapChar)
    }
    # The min number of tips having unambiguous characters for a site
    tipNum <- length(attr(x, "tree")[["tip.label"]])
    if (minSkipSize < 1) {
        minUnambiguous <- tipNum - tipNum * minSkipSize
    } else {
        minUnambiguous <- tipNum - minSkipSize
    }
    # The aligned sequences, reference numbering and loci
    align <- attr(x, "align")
    reference <- attr(x, "msaNumbering")
    loci <- attr(x, "loci")
    # Test the loci meet the skip size threshold
    lociOK <- sapply(loci, function(i) {
        siteSummary <- sitePath:::tableAA(align, reference[i] - 1)
        siteChars <- names(siteSummary)
        # All the unambiguous characters in the site summary
        unambiguousChars <-
            siteChars[which(siteChars %in% unambiguous)]
        sum(siteSummary[unambiguousChars]) > minUnambiguous
    })
    expect_true(all(lociOK))
    # Test the non-loci doesn't meet the skip size threshold
    nonLociOK <- sapply(
        X = setdiff(seq_along(reference), loci),
        FUN = function(i) {
            siteSummary <- sitePath:::tableAA(align, reference[i] - 1)
            siteChars <- names(siteSummary)
            # All the unambiguous characters in the site summary
            unambiguousChars <-
                siteChars[which(siteChars %in% unambiguous)]
            # The site is skipped also when it's completely conserved
            if (length(unambiguousChars) == 1) {
                return(TRUE)
            }
            sum(siteSummary[unambiguousChars]) <= minUnambiguous
        }
    )
    expect_true(all(nonLociOK))
}

test_that("The function works for amino acid", {
    data(zikv_align)
    data(zikv_tree)
    tr <- addMSA(zikv_tree, alignment = zikv_align)
    tipNames <- zikv_tree[["tip.label"]]
    align <- attr(tr, "align")
    # Find the index of sequence which has gap character
    for (i in grep('-', align)) {
        # Select a tip as the arbitrary reference
        tip <- tipNames[[i]]
        refSeq <-
            toupper(zikv_align[["seq"]][[which(zikv_align[["nam"]] == tip)]])
        # The reference sequence from the raw data and processed data should be
        # the same
        expect_identical(refSeq, align[[i]])
        # Test the numbering on 'phyMSAmatched' object
        test_referenceNumbering(tr, tip, refSeq, '-')
        test_siteSkipping(tr, tip, "C", 0.8)
        # Test the numbering on the 'lineagePath' object
        p <- lineagePath(tr)
        test_referenceNumbering(p, tip, refSeq, '-')
        test_siteSkipping(p, tip, "G", 0.9)
    }
})

test_that("The function works for nucleotide", {
    data(sars2_align)
    data(sars2_tree)
    tr <- addMSA(sars2_tree,
                 alignment = sars2_align,
                 seqType = "DNA")
    tipNames <- sars2_tree[["tip.label"]]
    align <- attr(tr, "align")
    # Find the index of sequence which has gap character
    for (i in grep('-', align)) {
        # Select a tip as the arbitrary reference
        tip <- tipNames[[i]]
        refSeq <-
            toupper(sars2_align[["seq"]][[which(sars2_align[["nam"]] == tip)]])
        # The reference sequence from the raw data and processed data should be
        # the same
        expect_identical(refSeq, align[[i]])
        # Test the numbering on 'phyMSAmatched' object
        test_referenceNumbering(tr, tip, refSeq, '-')
        test_siteSkipping(tr, tip, "C", 0.8)
        # Test the numbering on the 'lineagePath' object
        p <- lineagePath(tr)
        test_referenceNumbering(p, tip, refSeq, '-')
        test_siteSkipping(p, tip, "G", 0.9)
    }
})
