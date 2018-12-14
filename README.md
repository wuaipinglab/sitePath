sitePath
========
An R package for finding site with lineage-dependent fixation in RNA virus

Installation
------------
To install the package from github:
```r
install.packages("devtools")
devtools::install_github("Takkoona/sitePath")
```
To use plotting function and read annotated nexus file (from beast),  [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) is needed.

### Linux
Use the code above to install

### Windows
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) is needed for the installation. (There might be no compatible Rtools for a newly released R, in which case an [older R  version](https://cran.r-project.org/bin/windows/base/old/) is required)

### Mac
Unable to install from github so far

Example
-------
```r
library(sitePath)
tree <- ape::read.tree(
  system.file("ZIKV.newick", package = "sitePath")
)
outgroup <- readLines(system.file("ZIKV_outgroup.txt", package = "sitePath"))
tree <- ape::root(tree, outgroup)
align <- seqinr::read.alignment(
  system.file("ZIKV.fasta", package = "sitePath"),
  format = "fasta"
)
(paths <- sitePath(tree, align, 0.996))
mutations <- fixationSites(paths)
(predSites <- as.numeric(sapply(
  names(mutations), function(m) {substr(m, 2, nchar(m) - 1)}
)))
```
