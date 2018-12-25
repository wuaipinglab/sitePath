## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE-------------------------------------------------------
library(ape)
library(ggtree)

zikv_tree <- system.file("ZIKV.newick", package = "sitePath")
tree <- read.tree(zikv_tree)

## ------------------------------------------------------------------------
tree <- root(tree, "ANK57896")

## ------------------------------------------------------------------------
beast_file <- system.file("beast_mcc.tree", package="sitePath")
treeS4 <- read.beast(beast_file)
treeS3 <- treeS4@phylo

## ----message=FALSE-------------------------------------------------------
library(seqinr)

align <- read.alignment(system.file("ZIKV.fasta", package = "sitePath"), "fasta")

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

## ---- fig.show="hold"----------------------------------------------------
plot(fixations, "S139N")
plot(fixations, names(fixations)[6])

## ---- fig.show="hold"----------------------------------------------------
plot(fixations, 139)
plot(fixations, 763)

## ---- fig.show="hold"----------------------------------------------------
snps <- SNPsites(tree, align)
plot(fixations, snps[4])
plot(fixations, snps[5])

