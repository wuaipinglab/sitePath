as.data.frame.fixationSites <- function(x, ...) {
    cat("Under development\n")
}

as.data.frame.multiFixationSites <- function(x, ...) {
    cat("Under development\n")
}

#' @importFrom ape as.phylo
#' @export
ape::as.phylo

#' @export
as.phylo.lineagePath <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.fixationSites <- function(x, ...) {
    paths <- attr(x, "paths")
    res <- as.phylo.lineagePath(paths)
    return(res)
}

#' @export
as.list.sitewiseClusters <- function(x, ...) {
    groupName <- names(x)
    attributes(x) <- NULL
    res <- lapply(x, function(tips) {
        attributes(tips) <- NULL
        return(tips)
    })
    names(res) <- groupName
    return(res)
}
