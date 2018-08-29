sitePath
========
An R package for weird phylogenetic analysis

Installation
------------
To install the package from github:
```
install.packages("devtools")
devtools::install_github("Takkoona/sitePath", ref = "ancestor")
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
library(sitePath)

tree <- ape::read.tree(
  system.file("zika.tree", package = "sitePath")
)
align <- seqinr::read.alignment(
  system.file("zika.fasta", package = "sitePath"), "fasta"
)

sitePath <- getSitePath(tree, align, 0.996)
mutations <- findFixed(sitePath)
```
