#' @title Find sites
#' @description
#' To find site in sitePath
#' @param paths a \code{sitePath} object
#' @param n the index of the target path
#' @param allseq reconstructed sequences
#' @importFrom ape nodepath
#' @importFrom phangorn Descendants phyDat2alignment
#' @return sites and tip names
#' @export
findSites.sitePath <- function(paths, n, allseq) {
  tree <- attr(paths, "tree")
  attr(paths, "tree") <- NULL
  if (!is(allseq, "phyDat")) {
    stop("allseq is not a class of phyDat")
  }
  allseq <- toupper(phyDat2alignment(allseq)$seq)
  fixed <- lapply(paths[[n]], function(node) {
    tips <- unlist(Descendants(x = tree, node = node, type = "tips"))
    ancestors <- strsplit(allseq[nodepath(phy = tree, from = paths[[n]][1], to = node)], "")
    children <- strsplit(allseq[tips], "")
    mutations <- character(0)
    for (i in 1:length(children[[1]])) {
      a <- unique(sapply(ancestors, "[[", i))
      c <- unique(sapply(children, "[[", i))
      if (length(a) == 1 && length(c) == 1 && a != c) {
        mutations <- c(mutations, paste(a, i, c, sep = ""))
      }
    }
    if (length(mutations) == 0) {
      return(character(0))
    } else {
      return(list(tips = tree$tip.label[tips], mutations = mutations))
    }
  })
  names(fixed) <- paths[[n]]
  return(fixed[which(lengths(fixed) != 0)])
}

#' @export
findSites <- function(x, ...) UseMethod("findSites")
