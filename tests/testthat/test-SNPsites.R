test_that("Works for amino acid", {
    data(zikv_align_reduced)
    data(zikv_tree_reduced)
    tr <- addMSA(tree = zikv_tree_reduced,
                 alignment = zikv_align_reduced)
    tipNames <- zikv_tree_reduced[["tip.label"]]
    nTips <- length(tipNames)
    # Test the input of 'minSNP'
    expect_error(SNPsites(tr, minSNP = nTips))
    expect_error(SNPsites(tr, minSNP = "notRight"))
    # Set an arbitrary 'gapChar'
    gapChar <- "G"
    tr <- setSiteNumbering(tr, gapChar = gapChar)
    unambiguous <- setdiff(sitePath:::AA_UNAMBIGUOUS, gapChar)
    align <- attr(tr, "align")
    # Set up a range of 'minSNP' values
    minT <- floor(nTips / 10)
    maxT <- ceiling(nTips / 2)
    for (minSNP in seq(minT, maxT, 20)) {
        snp <- SNPsites(tr, minSNP = minSNP)
        allSNP <- attr(snp, "allSNP")
        allSNP <- split(x = allSNP,
                        f = allSNP[, c("Pos", "SNP")],
                        drop = TRUE)
        for (siteSNP in allSNP) {
            siteName <- unique(siteSNP[["Pos"]])
            snpAA <- unique(siteSNP[["SNP"]])
            # The amino acid should not be gap or ambiguous
            expect_true(all(snpAA %in% unambiguous))
            # At least two kinds of amino acid over the 'minSNP' threshold
            snpNum <- sapply(snpAA, function(AA) {
                sum(siteSNP[["SNP"]] == AA)
            })
            expect_true(all(snpNum >= minSNP))
            # Summarize the amino acid of the SNP site and the number should be
            # consistent with the result
            siteCharSummary <- table(sapply(
                X = align,
                FUN = substring,
                first = siteName,
                last = siteName
            ))
            expect_true(all(siteCharSummary[snpAA] == snpNum))
        }
    }
})

test_that("Works for nucleotide", {
    data(sars2_align)
    data(sars2_tree)
    tr <- addMSA(tree = sars2_tree,
                 alignment = sars2_align)
    tipNames <- sars2_tree[["tip.label"]]
    nTips <- length(tipNames)
    # Test the input of 'minSNP'
    expect_error(SNPsites(tr, minSNP = nTips))
    expect_error(SNPsites(tr, minSNP = "notRight"))
    # Set an arbitrary 'gapChar'
    gapChar <- "G"
    tr <- setSiteNumbering(tr, gapChar = gapChar)
    unambiguous <- setdiff(sitePath:::NT_UNAMBIGUOUS, gapChar)
    align <- attr(tr, "align")
    # Set up a range of 'minSNP' values
    minT <- floor(nTips / 10)
    maxT <- ceiling(nTips / 2)
    for (minSNP in seq(minT, maxT, 20)) {
        snp <- SNPsites(tr, minSNP = minSNP)
        allSNP <- attr(snp, "allSNP")
        allSNP <- split(x = allSNP,
                        f = allSNP[, c("Pos", "SNP")],
                        drop = TRUE)
        for (siteSNP in allSNP) {
            siteName <- unique(siteSNP[["Pos"]])
            snpAA <- unique(siteSNP[["SNP"]])
            # The amino acid should not be gap or ambiguous
            expect_true(all(snpAA %in% unambiguous))
            # At least two kinds of amino acid over the 'minSNP' threshold
            snpNum <- sapply(snpAA, function(AA) {
                sum(siteSNP[["SNP"]] == AA)
            })
            expect_true(all(snpNum >= minSNP))
            # Summarize the amino acid of the SNP site and the number should be
            # consistent with the result
            siteCharSummary <- table(sapply(
                X = align,
                FUN = substring,
                first = siteName,
                last = siteName
            ))
            expect_true(all(siteCharSummary[snpAA] == snpNum))
        }
    }
})
