#' @export
findSites.sitePath <- function(paths, type = c("isolated", "fused"), n) {
  findFixed <- function(index, path, nodeTips) {
    res <- character(0)
    
    before <- nodeTips[as.character(path[1:(index - 1)])]
    if (all(is.na(names(before)))) return(res)
    before <- unlist(before[!is.na(names(before))])
    before <- strsplit(before, "")
    
    after <- nodeTips[as.character(path[index:length(path)])]
    if (all(is.na(names(after)))) return(res)
    after <- unlist(after[!is.na(names(after))])
    after <- strsplit(after, "")
    
    for (i in 1:length(before[[1]])) {
      b <- unique(sapply(before, "[[", i))
      a <- unique(sapply(after, "[[", i))
      if (length(b) == 1 && length(a) == 1 && a != b) {
        res <- c(res, paste(b, i, a, sep = ""))
      }
    }
    attr(res, "from") <- names(before)
    attr(res, "to") <- names(after)
    return(res)
  }
  path <- paths[[type]][[n]]
  fixed <- lapply(
    1:length(path), findFixed, path,
    attr(paths, "nodeTips")
  )
  names(fixed) <- path
  return(fixed)
}

#' @export
findSites <- function(x, ...) UseMethod("findSites")
