library(Rcpp)
library(phangorn)
library(seqinr)
library(Rgraphviz)

tree <- read.tree(
  "/home/chengyang/PycharmProjects/FLUSITE/tests/dummy/RAxML_bestTree.mixed"
)
tipPath <- lapply(nodepath(tree), rev)

aligns <- read.alignment(
  "/home/chengyang/PycharmProjects/FLUSITE/tests/dummy/Alignment.fasta",
  "fasta"
)
seqs <- strsplit(
  sapply(aligns$seq, function(seq) {seq}), ""
)[match(tree$tip.label, aligns$nam)]
if (anyNA(seqs)) {
  stop("Size not matched for tree tips and alignment")
}

seqsAR <- as.phyDat(aligns, "AA")
fit = pml(tree, seqsAR)
fit = optim.pml(fit, model = "JTT")
anc.ml = ancestral.pml(fit, "ml")
ancestral2phyDat <- function(x) {
  eps <- 1.0e-5
  contr <- attr(x, "contrast")
  ind1 <- which( apply(contr, 1, function(x)sum(x > eps)) == 1L)
  ind2 <- which( contr[ind1, ] > eps, arr.ind = TRUE)
  pos <- ind2[match(seq_len(ncol(contr)),  ind2[,2]),1]
  res <- lapply(x, function(x, pos) pos[max.col(x)], pos)
  attributes(res) <- attributes(x)
  return(res)
}
alignsAR <- phyDat2alignment(ancestral2phyDat(anc.ml))
seqsAR <- strsplit(alignsAR$seq, "")

sourceCpp("treemer.cpp")
grouping <- trimTree(tipPath, seqs, 0.95); length(grouping)

mutations <- mutationPath(tipPath, seqsAR, 0.95)

edgeLabels <- names(mutations)
edges <- strsplit(edgeLabels, "~")

g <- new(
  "graphNEL",
  nodes=unique(unlist(edges)),
  edgemode="directed"
)

for (edge in edges) {g <- addEdge(edge[1], edge[2], g)}

eAttrs <- list()
eAttrs$label <- sapply(mutations, function(m) {paste(m, collapse = "\\\n")})
plot(g, edgeAttrs=eAttrs)
