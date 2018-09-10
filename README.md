sitePath
========
An R package for weird phylogenetic analysis

Installation
------------
To install the package from github:
```
install.packages("devtools")
devtools::install_github("Takkoona/sitePath")
```
To use plotting function and read annotated nexus file (from beast). [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) is needed.

### Linux
Use the code above to install

### Windows
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) is needed for the installation. (There might be no compatible Rtools for a newly released R, in which case an [older R  version](https://cran.r-project.org/bin/windows/base/old/) is required)

### Mac
Unable to install from github so far

Example
-------
```
tree <- ggtree::read.beast(
  system.file("m.trees", package = "sitePath")
)@phylo

align <- seqinr::read.alignment(
  system.file("m.aligned.fasta", package = "sitePath"),
  format = "fasta", forceToLower = FALSE
)

matched <- treeAlignMatch(tree, align)
grouping <- groupTips(matched, 0.98)
sitePath <- getSitePath(matched, 0.98)
findSites(sitePath, 1)
```
