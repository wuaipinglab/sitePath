#' @export
findSites.sitePath <- function(paths, n) {
  fixed <- findFixed(paths[[n]], paths, attr(paths, "nodeTips"))
  return(fixed)
}

#' @export
findSites <- function(x, ...) UseMethod("findSites")

findFixed.isolated <- function(path, paths, nodeTips) {
  findMutation <- function(index, path, nodeTips) {
    before <- nodeTips[as.character(path[1:(index - 1)])]
    if (all(is.na(names(before)))) return(character(0))
    before <- before[!is.na(names(before))]
    names(before) <- NULL
    before <- strsplit(unlist(before), "")
    
    after <- nodeTips[as.character(path[index:length(path)])]
    if (all(is.na(names(after)))) return(character(0))
    after <- after[!is.na(names(after))]
    names(after) <- NULL
    after <- strsplit(unlist(after), "")
    
    mutations <- character(0)
    for (i in 1:length(before[[1]])) {
      b <- unique(sapply(before, "[[", i))
      a <- unique(sapply(after, "[[", i))
      if (length(b) == 1 && length(a) == 1 && a != b) {
        mutations <- c(mutations, paste(b, i, a, sep = ""))
      }
    }
    if (length(mutations) == 0) {
      return(character(0))
    } else {
      return(list(from  = names(before), to = names(after), mutations = mutations))
    }
  }
  res <- lapply(1:length(path), findMutation, path, attr(paths, "nodeTips"))
  names(res) <- path
  return(res[which(lengths(res) != 0)])
}

findFixed.fused <- function(path, paths, nodeTips) {
  findMutation <- function(index, path, nodeTips) {
    before <- nodeTips[as.character(path[1:(index - 1)])]
    if (all(is.na(names(before)))) return(character(0))
    before <- before[!is.na(names(before))]
    names(before) <- NULL
    before <- strsplit(unlist(before), "")
    
    after <- nodeTips[as.character(c(path[index:length(path)], attr(path, "endNodes")))]
    after <- after[!is.na(names(after))]
    names(after) <- NULL
    after <- strsplit(unlist(after), "")
    
    mutations <- character(0)
    for (i in 1:length(before[[1]])) {
      b <- unique(sapply(before, "[[", i))
      a <- unique(sapply(after, "[[", i))
      if (length(b) == 1 && length(a) == 1 && a != b) {
        mutations <- c(mutations, paste(b, i, a, sep = ""))
      }
    }
    if (length(mutations) == 0) {
      return(character(0))
    } else {
      return(list(from  = names(before), to = names(after), mutations = mutations))
    }
  }
  res <- lapply(1:length(path), findMutation, path, attr(paths, "nodeTips"))
  names(res) <- path
  return(res[which(lengths(res) != 0)])
}

findFixed <- function(x, ...) UseMethod("findFixed")
