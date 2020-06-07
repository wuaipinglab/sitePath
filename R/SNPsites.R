#' @rdname SNPsites
#' @title Finding sites with variation
#' @description Single nucleotide polymorphism (SNP) in the whole package refers
#'   to variation of amino acid. \code{SNPsite} will try to find SNP in the
#'   multiple sequence alignment. A reference sequence and gap character may be
#'   specified to number the site.
#' @param tree The return from \code{\link{addMSA}} function.
#' @param minSNP Minimum number of amino acid variation to be a SNP.
#' @return A \code{SNPsites} object
#' @importFrom stats complete.cases
#' @export
#' @examples
#' data(zikv_tree_reduced)
#' data(zikv_align_reduced)
#' tree <- addMSA(zikv_tree_reduced, alignment = zikv_align_reduced)
#' SNPsites(tree)
SNPsites <- function(tree, minSNP = NULL) {
    align <- attr(tree, "align")
    if (is.null(align)) {
        stop("No alignment found in \"tree\"")
    } else if (length(unique(nchar(align))) != 1) {
        stop("Sequence length not equal")
    }
    if (is.null(minSNP)) {
        minSNP <- length(tree[["tip.label"]]) / 10
    }
    sites <- attr(tree, "reference")
    majorSNP <- vapply(
        X = sites,
        FUN = function(s) {
            aaSummary <- tableAA(align, s - 1)
            gapIndex <- which(names(aaSummary) == '-')
            if (length(gapIndex) != 0) {
                aaSummary <- aaSummary[-gapIndex]
            }
            names(aaSummary)[which.max(aaSummary)]
        },
        FUN.VALUE = character(1)
    )
    snpColors <- vapply(
        X = AA_FULL_NAMES,
        FUN = function(i) {
            AA_COLORS[[i]]
        },
        FUN.VALUE = character(1)
    )
    names(snpColors) <- toupper(names(snpColors))
    allSNP <- lapply(names(align), function(ac) {
        seq <- align[[ac]]
        res <- as.data.frame(t(vapply(
            X = seq_along(majorSNP),
            FUN = function(site) {
                s <- as.integer(site)
                snp <- substr(seq, s, s)
                if (snp != majorSNP[[site]] && snp != '-') {
                    return(c(site, snp))
                } else {
                    return(c(site, NA_character_))
                }
            },
            FUN.VALUE = c(integer(1), character(1))
        )))
        res <- res[complete.cases(res), ]
        colnames(res) <- c("Pos", "SNP")
        res[["Pos"]] <- as.integer(res[["Pos"]])
        res[["Accession"]] <- rep(ac, nrow(res))
        return(res)
    })
    allSNP <- do.call(rbind, allSNP)
    allSNP <- allSNP[, c(3, 1, 2)]
    snpSummary <- as.data.frame(table(allSNP[["Pos"]],
                                      allSNP[["SNP"]]))
    allSNP <- merge(
        x = allSNP,
        y = snpSummary,
        by.x = c("Pos", "SNP"),
        by.y = c("Var1", "Var2"),
        all.x = TRUE
    )
    allSNP <- allSNP[which(allSNP[, "Freq"] >= minSNP),
                     c(3, 1, 2, 4)]
    rownames(allSNP) <- NULL
    res <- sort(unique(allSNP[["Pos"]]))
    attr(res, "allSNP") <- allSNP
    attr(res, "tree") <- tree
    class(res) <- "SNPsites"
    return(res)
}

#' @export
print.SNPsites <- function(x, ...) {
    x <- as.integer(x)
    print(x)
}
