#' @title Format conversion
#' @description
#' The function reformats a \code{\link{phyloMutations}} object.
#' Instead of sorting by linked nodes, the mutations are summarized by sites 
#' @param mutations A \code{\link{phyloMutations}} object
#' @seealso \code{\link{toTips}}
#' @return site and corresponding nodes
#' @export

toNodes <- function(mutations) {
  res <- lapply(mutations, function(sites) {
    linksAsList <- branch2Site(attr(mutations, "tree"), sites)
    return(lapply(linksAsList, unlist))
  })
  return(res)
}

branch2Site <- function(tree, sites) {
  res <- list()
  for (link in names(sites)) {
    muts <- sites[[link]]
    if (length(muts) != 0){
      nodes <- unlist(strsplit(link, "~"))
      for (m in muts) {
        len <- nchar(m)
        mfrom <- substr(m, 1, 1)
        mto <- substr(m, len, len)
        pos <- substr(m, 2, len - 1)
        alleles <- c(mfrom, mto)
        names(alleles) <- nodes
        if (length(res[[pos]]) == 0) {
          res[[pos]] <- list(alleles)
        } else {
          res[[pos]] <- c(res[[pos]], list(alleles))
        }
      }
    }
  }
  return(res)
}

#' @title Format conversion
#' @description
#' Reformat a \code{\link{phyloMutations}} object.
#' The whole tree will be chopped at the branches where mutation sits. 
#' Descendants derived from nodes on the same tree segment are grouped.
#' @seealso \code{\link{toNodes}}
#' @param mutations A \code{\link{phyloMutations}} object
#' @return site and corresponding tips
#' @export

toTips <- function(mutations) {
  tree <- attr(mutations, "tree")
  root <- getRoot(tree)
  endNodes <- sapply(attr(mutations, "evolPath"), function(ep) {
    return(ep[length(ep)])
  })
  pathNodes <- unlist(attr(mutations, "evolPath"))
  res <- lapply(mutations, function(sites) {
    groupings <- lapply(branch2Site(tree, sites), function(linkList) {
      links <- unlist(linkList)
      linkList <- lapply(linkList, function(l) {
        as.integer(names(l))
      })
      nodes <- as.integer(names(links))
      keepTos <- c(nodes[seq(1, length(links) - 1, 2)], endNodes)
      keepFroms <- c(nodes[seq(2, length(links), 2)], endNodes)
      segPath <- list()
      for (keepTo in keepTos) {
        segPath <- c(segPath, list(nodepath(tree, root, keepTo)))
        for (keepFrom in keepFroms) {
          segPath <- c(segPath, list(nodepath(tree, keepFrom, keepTo)))
        }
      }
      allSegPoint <- c(root, nodes, endNodes)
      qualified <- list()
      for (seg in segPath) {
        if (length(intersect(allSegPoint, seg)) <= 2 && !list(sort(seg)) %in% linkList) {
          qualified <- c(qualified, list(seg))
        }
      }
      combined <- list(qualified[[1]])
      for (seg in qualified[2:length(qualified)]) {
        found <- FALSE
        for (i in 1:length(combined)) {
          g <- combined[[i]]
          if (length(intersect(g, seg) != 0)) {
            combined[[i]] <- union(g, seg)
            found <- TRUE
            break
          }
        }
        if (!found) {
          combined <- c(combined, list(seg))
        }
      }
      grouping <- lapply(combined, function(seg) {
        tips <- integer(0)
        for (n in seg) {
          for (child in Children(tree, n)) {
            if (!child %in% pathNodes) {
              if (child <= length(tree$tip.label)) {
                tips <- c(tips, child)
              } else {
                tips <- c(tips, unlist(Descendants(tree, child, type = "tips")))
              }
            }
          }
        }
        return(tree$tip.label[tips])
      })
      return(grouping)
    })
    return(groupings)
  })
  return(res)
}
