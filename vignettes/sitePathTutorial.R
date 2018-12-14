## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(ape)
tree <- read.tree(system.file("ZIKV.newick", package = "sitePath"))

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(ggtree))
tree <- read.tree(system.file("ZIKV.newick", package = "sitePath"))

## ------------------------------------------------------------------------
outgroup <- readLines(system.file("ZIKV_outgroup.txt", package = "sitePath"))
print(outgroup)
tree <- ape::root(tree, outgroup)

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(seqinr))
align <- read.alignment(system.file("ZIKV.fasta", package = "sitePath"), "fasta")

## ------------------------------------------------------------------------
suppressPackageStartupMessages(library(phangorn))
align <- read.phyDat(system.file("ZIKV.fasta", package = "sitePath"), format = "fasta", type = "AA")
align <- phyDat2alignment(align)

## ------------------------------------------------------------------------
treeS4 <- read.beast(system.file("extdata/BEAST", "beast_mcc.tree", package="ggtree"))
treeS3 <- treeS4@phylo

## ---- fig.width=6, fig.height=4------------------------------------------
library(sitePath)
preassessment <- sneakPeek(tree, align)

## ------------------------------------------------------------------------
paths <- sitePath(tree, align, 0.996)
paths

## ------------------------------------------------------------------------
grouping <- groupTips(tree, align, 0.996)

## ------------------------------------------------------------------------
fixations <- fixationSites(paths)
fixations

## ------------------------------------------------------------------------
plot(fixations, "S139N")

## ------------------------------------------------------------------------
plot(fixations, names(fixations)[6])

## ---- fig.show="hold"----------------------------------------------------
plot(fixations, 139)
plot(fixations, 763)

## ---- fig.show="hold"----------------------------------------------------
snps <- SNPsites(tree, align)
plot(fixations, snps[4])
plot(fixations, snps[5])

