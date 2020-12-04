test_ambiguousIgnored <- function(fixSites,
                                  tree,
                                  alignment,
                                  gapChar,
                                  minEffectiveSize) {
    paths <- attr(fixSites, "paths")
    tipNames <- tree[["tip.label"]]
    nTips <- length(tipNames)
    p <- attr(fixSites, "paths")
    if (is.null(minEffectiveSize)) {
        minEffectiveSize <- attr(paths, "minSize")
    }
    # Get the divergent nodes
    divNodes <- sitePath:::divergentNode(paths)
    # The tips and matching
    pathNodeTips <-
        sitePath:::.tipSeqsAlongPathNodes(paths, divNodes)
    # In case root node does not have any tips
    excludedNodes <- divNodes
    rootNode <- attr(paths, "rootNode")
    if (!rootNode %in% names(pathNodeTips)) {
        excludedNodes <- c(rootNode, excludedNodes)
    }
    pathsWithSeqs <- lapply(paths, function(path) {
        path <- as.character(setdiff(path, excludedNodes))
        names(path) <- path
        pathNodeAlign <- pathNodeTips[path]
        attr(pathNodeAlign, "pathTipNum") <-
            sum(lengths(pathNodeAlign))
        return(pathNodeAlign)
    })
    pathTipNums <-
        vapply(pathsWithSeqs, attr, integer(1), "pathTipNum")
    minEffectiveSize <-
        ceiling(minEffectiveSize * (min(pathTipNums) / max(pathTipNums)))
    # Get all the unambiguous character
    if (attr(p, "seqType") == "AA") {
        unambiguous <- setdiff(sitePath:::AA_UNAMBIGUOUS, gapChar)
    } else {
        unambiguous <- setdiff(sitePath:::NT_UNAMBIGUOUS, gapChar)
    }
    for (sp in fixSites) {
        for (mp in sp) {
            # The number of tips should be less than total tree tips
            expect_lte(sum(lengths(mp)), nTips)
            # AA/NT summary of the tips on the 'mutPath'
            aa <- lapply(mp, function(tips) {
                site <- attr(sp, "site")
                matchIndex <-
                    which(alignment[["nam"]] %in% tipNames[tips])
                sum <- alignment[["seq"]][matchIndex]
                sum <- table(sapply(sum, substring, site, site))
            })
            unambiguousNum <- 0
            for (i in seq_along(mp)) {
                # The tips in each group before or after fixation mutation
                nodeTips <- mp[[i]]
                # Tips should be unique and more than 'minEffectiveSize'
                expect_equal(sort(nodeTips), sort(unique(nodeTips)))
                # The dominant AA/NT should be the fixed one
                aaD <- toupper(names(which.max(aa[[i]])))
                fixedAA <- attr(nodeTips, "AA")
                attributes(fixedAA) <- NULL
                if (fixedAA %in% unambiguous) {
                    unambiguousNum <- unambiguousNum + 1
                }
                if (is.null(attr(nodeTips, "toMerge"))) {
                    expect_gte(length(nodeTips), minEffectiveSize)
                    expect_equal(fixedAA, aaD)
                }
            }
            expect_gte(unambiguousNum, 2)
        }
    }
}

test_that("The function works for amino acid", {
    data(h3n2_tree_reduced)
    data(h3n2_align_reduced)
    tr <- addMSA(tree = h3n2_tree_reduced,
                 alignment = h3n2_align_reduced)
    p <- lineagePath(tr)
    # Set an arbitrary gap character
    gapChar <- "R"
    p <- setSiteNumbering(p, gapChar = gapChar)
    # Test the input of 'minEffectiveSize' and 'searchDepth'
    expect_error(fixationSites(paths = p, minEffectiveSize = -1))
    expect_error(fixationSites(paths = p, minEffectiveSize = "3"))
    expect_error(fixationSites(paths = p, searchDepth = -1))
    # Here comes the real deal
    fixSites <- fixationSites(p)
    minEntropy <- sitesMinEntropy(p)
    # Test the two function give the same result
    expect_identical(fixSites, fixationSites(minEntropy))
    expect_false(any(duplicated(unlist(
        attr(fixSites, "clustersByPath")
    ))))
    # Test the number of tips before and after each fixation mutation is enough
    test_ambiguousIgnored(fixSites,
                          h3n2_tree_reduced,
                          h3n2_align_reduced,
                          gapChar,
                          NULL)
})

test_that("The function works for nucleotide", {
    data(sars2_align)
    data(sars2_tree)
    tr <- addMSA(sars2_tree,
                 alignment = sars2_align,
                 seqType = "DNA")
    # Set an arbitrary gap character
    gapChar <- "G"
    tr <- setSiteNumbering(tr, gapChar = gapChar)
    p <- lineagePath(tr)
    fixSites <- fixationSites(p)
    minEntropy <- sitesMinEntropy(p)
    # Test the two function give the same result
    expect_identical(fixSites, fixationSites(minEntropy))
    # Test the number of tips before and after each fixation mutation is enough
    test_ambiguousIgnored(fixSites,
                          sars2_tree,
                          sars2_align,
                          gapChar,
                          NULL)
})
