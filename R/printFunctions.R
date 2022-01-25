#' @export
print.alignment <- function(x, ...) {
    cat("MSA with", length(x[["seq"]]), "sequences.\n\n")
    cat("Sequence names:\n")
    cat("  ", paste(head(x[["nam"]]), collapse = ", "), ", ...\n\n")
    cat("Sequence length:", nchar(x[["seq"]][[1]]), "\n")
}

#' @export
print.phyMSAmatched <- function(x, ...) {
    cat("This is a 'phyMSAmatched' object.\n")
}

#' @export
print.lineagePath <- function(x, ...) {
    cat(
        "This is a 'lineagePath' object.\n\n",
        length(x),
        " lineage paths using ",
        attr(x, "minSize"),
        " as \"major SNP\" threshold \n",
        sep = ""
    )
}

#' @export
print.SNPsites <- function(x, ...) {
    cat("This is a 'SNPsites' object.", "\n\n")
    x <- allSitesName(x)
    cat(paste(x, collapse = " "))
}

#' @export
print.sitesMinEntropy <- function(x, ...) {
    cat("This is a 'sitesMinEntropy' object.", "\n\n")
    loci <- allSitesName.sitesMinEntropy(x)
    cat("There are", length(loci), "sites with enough variation.\n")
}

#' @export
print.fixationSites <- function(x, ...) {
    cat("This is a 'fixationSites' object.\n\nResult for",
        length(attr(x, "paths")),
        "paths:\n\n")
    if (length(x) == 0) {
        cat("No multi-fixation found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(x, "reference")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}

#' @export
print.sitePath <- function(x, ...) {
    cat("Site",
        attr(x, "site"),
        "may experience fixation on",
        length(x),
        "path(s):\n\n")
    # A 'sitePath' consists of all the fixation paths for a single site. So each
    # 'm' represent a single fixation path
    for (m in x) {
        if (length(m) == 2) {
            mutName <-
                paste0(attr(m[[1]], "AA"), attr(x, "site"), attr(m[[2]], "AA"))
            cat(mutName,
                paste0("(", length(m[[1]]), "->", length(m[[2]]), ")"),
                "\n")
        } else {
            mutName <- character(0)
            for (tips in m) {
                aa <- attr(tips, "AA")
                mutName <-
                    c(mutName, paste0(aa, "(", length(tips), ")"))
            }
            cat(paste0(mutName, collapse = " -> "), "\n")
        }
    }
    cat("\nIn the bracket are the number of tips",
        "involved before and after the fixation\n")
}

#' @export
print.parallelSites <- function(x, ...) {
    cat("This is a 'parallelSites' object.\n\nResult for",
        length(attr(x, "paths")),
        "paths:\n\n")
    if (length(x) == 0) {
        cat("No parallel site found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(x, "reference")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}

#' @export
print.sitePara <- function(x, ...) {
    cat(
        "This is a 'sitePara' object.\n\nSite",
        attr(x, "site"),
        "may have parallel mutation on",
        length(x),
        "pair of paths:\n\n"
    )
    mutSummary <- table(vapply(
        X = extractTips.sitePara(x),
        FUN = function(mutTips) {
            attr(mutTips, "mutName")[4]
        },
        FUN.VALUE = character(1)
    ))
    mutInfo <- character()
    for (mutName in names(mutSummary)) {
        mutInfo <-
            c(mutInfo, paste0(mutName, "(", mutSummary[[mutName]], ")"))
    }
    cat(
        paste0(mutInfo, collapse = ", "),
        "\n\nIn the bracket are the number of tips",
        "involved in the mutation\n"
    )
}

#' @export
print.paraFixSites <- function(x, ...) {
    cat("This is a 'paraFixSites' object\n\n")
    for (type in c("fixation", "parallel", "paraFix")) {
        sites <- allSitesName.paraFixSites(x, type = type)
        if (length(sites)) {
            cat(type,
                " sites:\n",
                paste0(sites, collapse = ", "),
                "\n\n", sep = "")
        } else {
            cat("No", type, "sites found.\n\n")
        }
    }
}

#' @export
print.fixationIndels <- function(x, ...) {
    cat("This is a 'fixationIndels' object.\n\nResult for",
        length(attr(x, "paths")),
        "paths:\n\n")
    if (length(x) == 0) {
        cat("No multi-fixation found\n")
    } else {
        cat(paste(names(x), collapse = " "), "\n")
        refSeqName <- attr(x, "reference")
        if (is.null(refSeqName)) {
            cat("No reference sequence specified.",
                "Using alignment numbering\n")
        } else {
            cat("Reference sequence: ", refSeqName, "\n", sep = "")
        }
    }
}

#' @export
print.indelPath <- function(x, ...) {
    cat(
        "Site(s)",
        attr(x, "indelSites"),
        "may experience fixation on",
        length(x),
        "group(s) of tips:\n\n"
    )
    for (i in x) {
        cat("(", length(i), ") ", sep = "")
    }
    cat("\nIn the bracket are the number of tips",
        "involved in the fixation\n")
}

#' @export
print.fixationPath <- function(x, ...) {
    print(names(x))
}
