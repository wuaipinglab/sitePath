#' @title ReadTreeAlign
#' @description Read tree and its matched alignment file
#' @param treeFile Path to the file stored tree string
#' @param treeFormat Specify tree format
#' @param alignFile Path to the file stored sequence alignment
#' @param alignFormat Specify alignment format
#' @return a treeAlignMatch object
#' @export

readTreeAlign <- function(
  treeFile, treeFormat = c("nexus", "newick", "beast"),
  alignFile, alignFormat = c("fasta", "clustal", "phylip", "mase", "msf")
) {
  tree <- switch (
    treeFormat,
    nexus = read.nexus(treeFile),
    newick = read.tree(treeFile),
    beast = ggtree::read.beast(treeFile)
  )
  align <- read.alignment(alignFile, alignFormat)
  if(!identical(sort(tree$tip.label), sort(align$nam))) {
    stop("tree tips and alignment names are not matched")
  }
  return(structure(
    list(tree = tree, align = align), 
    class = "treeAlignMatch"
  ))
}

#' @title Group Tips
#' @description Group tree tips by their sequence similarity
#' @param treeAlignMatch a \code{\link{treeAlignMatch}} object
#' @param similarity similarity threshold for tree trimming
#' @return grouping of tips
#' @export
#' @examples
#' \dontrun{
#' treeFile <- system.file("test.tree", package = "sitePath")
#' alignFile <- system.file("test.aligned.fasta", package = "sitePath")
#' matched <- readTreeAlign(treeFile, "newick", alignFile, "fasta")
#' grouping <- groupTips(matched, 0.95)
#' }

groupTips <- function(treeAlignMatch, similarity) {
  return(group(
    treeAlignMatch$tree,
    treeAlignMatch$align,
    similarity
  ))
}

group <- function(tree, align, similarity) {
  grouping <- trimTree(
    lapply(nodepath(tree), rev),
    strsplit(unlist(align$seq), "")[match(tree$tip.label, align$nam)], 
    similarity
  )
  res <- lapply(grouping, function(g) {
    return(tree$tip.label[g])
  })
  attr(res, "tree") <- tree
  return(res)
}

#' @title Ancestral Mutations
#' @description 
#' Find the mutations in predicted ancestral sequence and corresponding clades
#' @param treeAlignMatch a \code{\link{treeAlignMatch}} object
#' @param similarity similarity threshold for tree trimming
#' @param seqType specify sequence type ("DNA", "AA")
#' @param model substitution model for reconstruction ancestral sequence
#' @param siteMode specify the mutational mode of return
#' @return Linked clades and the mutations
#' @export
#' @examples
#' \dontrun{
#' treeFile <- system.file("test.tree", package = "sitePath")
#' alignFile <- system.file("test.aligned.fasta", package = "sitePath")
#' matched <- readTreeAlign(treeFile, "newick", alignFile, "fasta")
#' mutations <- ancestralMutations(matched, 0.95, "AA")
#' }

ancestralMutations <- function(
  treeAlignMatch, similarity, seqType = c("DNA", "AA"),
  model = NULL, siteMode = c(1, 2, 3)
) {
  fit <- pml(treeAlignMatch$tree, as.phyDat(treeAlignMatch$align, seqType))
  if (is.null(model)) {
    model <- switch (seqType, DNA = "GTR", AA = "JTT")
  }
  fit <- optim.pml(fit, model = model)
  anc.ml <- ancestral.pml(fit, type = "ml", return = "phyDat")
  if (length(anc.ml) < max(treeAlignMatch$tree$edge[,1])) {
    stop("The tree might be unrooted yet. Use ape::root() function to root.")
  }
  alignAR <- phyDat2alignment(anc.ml)
  res <- mutationPath(
    lapply(nodepath(treeAlignMatch$tree), rev),
    strsplit(alignAR$seq, ""),
    similarity, siteMode
  )
  attr(res, "tree") <- treeAlignMatch$tree
  return(res)
}

#' @title Format conversion
#' @description Reformat the return from \code{\link{ancestralMutations}}
#' @param mutations the return of \code{\link{ancestralMutations}} function
#' @return site and corresponding nodes
#' @export

mutations2Nodes <- function(mutations) {
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
#' @description Reformat the return from \code{\link{ancestralMutations}}
#' @param mutations the return of \code{\link{ancestralMutations}} function
#' @return site and corresponding tips
#' @export

mutations2Tips <- function(mutations) {
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

#' @title Mutation map
#' @param mutations the return of \code{\link{ancestralMutations}} function
#' @description plot the linked clades and mutations
#' @return NULL
#' @examples 
#' \dontrun{
#' treeFile <- system.file("test.tree", package = "sitePath")
#' alignFile <- system.file("test.aligned.fasta", package = "sitePath")
#' matched <- readTreeAlign(treeFile, "newick", alignFile, "fasta")
#' mutations <- ancestralMutations(matched, 0.95, "AA")
#' mutations2graphviz(mutations)
#' }
#' @export

mutations2graphviz <- function(mutations) {
  for (name in names(mutations)) {
    edges <- strsplit(names(mutations[[name]]), "~")
    g <- graph::graphNEL(nodes = unique(unlist(edges)), edgemode = "directed")
    for (edge in edges) {g <- graph::addEdge(edge[1], edge[2], g)}
    Rgraphviz::plot(
      g, main = name,
      edgeAttrs = list(
        label = sapply(mutations[[name]], function(m) {
          paste(m, collapse = " ")
        })
      ),
      attrs = list(
        node = list(fillcolor = "lightblue"),
        edge = list(color = "cyan")
      )
    )
  }
}

#' @title Internal node plot
#' @param mutations the return of \code{\link{ancestralMutations}} function
#' @description plot tree with its internal numbering
#' @return NULL
#' @export

nodePlot <- function(mutations) {
  p <- ggtree::ggtree(attr(mutations, "tree")) + 
    ggtree::geom_text2(ggplot2::aes(subset=!isTip, label=node), hjust=-.1)
  return(p)
}

#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL
