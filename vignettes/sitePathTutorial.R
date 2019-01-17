## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----import_tree, message=FALSE------------------------------------------
library(ape)
library(ggtree)

zikv_tree <- system.file("ZIKV.newick", package = "sitePath")
tree <- read.tree(zikv_tree)

## ----root_tree-----------------------------------------------------------
tree <- root(tree, "ANK57896")

## ----S4phylo-------------------------------------------------------------
beast_file <- system.file("beast_mcc.tree", package="sitePath")
treeS4 <- read.beast(beast_file)
treeS3 <- treeS4@phylo

## ----read_alignment, message=FALSE---------------------------------------
library(seqinr)

align <- read.alignment(system.file("ZIKV.fasta", package = "sitePath"), "fasta")

## ----sneakPeek_plot, fig.width=6, fig.height=4---------------------------
library(sitePath)

preassessment <- sneakPeek(tree, align)

## ----get_sitePath--------------------------------------------------------
paths <- sitePath(tree, align, 0.996)
paths

## ----find_fixations------------------------------------------------------
fixations <- fixationSites(paths)
fixations

## ----get_tipNames--------------------------------------------------------
fixations$S139N

## ---- plot_fixations, fig.show="hold", fig.width=4-----------------------
par(mar = c(1,1,1,1))
plot(fixations, "S139N")
plot(fixations, names(fixations)[6])

## ----group_tips----------------------------------------------------------
grouping <- groupTips(tree, align, 0.996)

## ----plot_sites, fig.show="hold", fig.width=4----------------------------
par(mar = c(1,1,1,1))
plot(fixations, 139)
plot(fixations, 763)

## ----find_SNP, fig.show="hold", fig.width=4------------------------------
snps <- SNPsites(tree, align)
par(mar = c(1,1,1,1))
plot(fixations, snps[4])
plot(fixations, snps[5])

