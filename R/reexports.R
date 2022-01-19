#' @importFrom ape as.phylo read.tree
#' @importFrom tidytree as.treedata

#' @export
ape::read.tree

#' @export
seqinr::read.alignment

#' @export
ape::as.phylo

#' @export
as.phylo.phyMSAmatched <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.sitePath <- function(x, ...) {
    res <- attr(x, "tree")
    return(res)
}

#' @export
as.phylo.sitesMinEntropy <- function(x, ...) {
    paths <- attr(x, "paths")
    res <- attr(paths, "tree")
    return(res)
}

#' @export
as.phylo.fixationSites <- function(x, ...) {
    paths <- attr(x, "paths")
    res <- attr(paths, "tree")
    return(res)
}

as.phylo.fixationIndels <- function(x, ...) {
    paths <- attr(x, "paths")
    res <- as.phylo(paths)
    return(res)
}

#' @export
tidytree::as.treedata

#' @export
as.treedata.fixationSites <- function(tree, ...) {
    extraArgs <- list(...)
    mutTable <- .mutationTable(tree)
    transMut <- lapply(X = split(mutTable, mutTable[, "node"]),
                       FUN = "[[",
                       i = "mutation")
    .node <- extraArgs[[".node"]]
    if (is.null(.node)) {
        .node <- groupTips.fixationSites(tree)
    }
    tree <- groupOTU(as.phylo.fixationSites(tree), .node)
    tree <- .annotateSNPonTree(tree, transMut)
    return(tree)
}

as.treedata.fixationIndels <- function(tree, ...) {
    for (sites in names(tree)) {
        indelPath <- tree[[sites]]
    }
}

#' @export
as.treedata.fixationPath <- function(tree, ...) {
    res <- attr(tree, "SNPtracing")
    return(res)
}
