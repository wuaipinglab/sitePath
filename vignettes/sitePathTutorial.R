## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----import_tree, message=FALSE------------------------------------------
library(ape)

tree <- read.tree(system.file("extdata", "ZIKV.newick", package = "sitePath"))

## ----root_tree-----------------------------------------------------------
tree <- root(tree, "ANK57896")

## ----add_alignment, message=FALSE----------------------------------------
library(sitePath)

alignment_file <- system.file("extdata", "ZIKV.fasta", package = "sitePath")
tree <- addMSA(tree, alignment_file, "fasta")

## ----sneakPeek_plot, fig.width=6, fig.height=4---------------------------
preassessment <- sneakPeek(tree)

## ----get_sitePath--------------------------------------------------------
paths <- sitePath(tree, 0.996)
paths

## ----find_fixations------------------------------------------------------
fixations <- fixationSites(paths)
fixations

## ----get_tipNames--------------------------------------------------------
fixations$S139N

## ----plot_fixations, fig.show="hold", fig.width=4------------------------
par(mar = c(1,1,1,1))
plot(fixations, "S139N")
plot(fixations, names(fixations)[6])

## ----group_tips----------------------------------------------------------
grouping <- groupTips(tree, 0.996)

## ----plot_sites, fig.show="hold", fig.width=4----------------------------
par(mar = c(1,1,1,1))
plot(fixations, 139)
plot(fixations, 763)

## ----find_SNP, fig.show="hold", fig.width=4------------------------------
snps <- SNPsites(tree)
par(mar = c(1,1,1,1))
plot(fixations, snps[4])
plot(fixations, snps[5])

