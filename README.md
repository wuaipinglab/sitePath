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
To use plotting function and read annotated nexus file (from beast). [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html) is needed.

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
align <- seqinr::read.alignment(
  system.file("zika_genome.fasta", package = "sitePath"), "fasta"
)

tree <- ape::read.tree(
  system.file("zika_genome.tree", package = "sitePath")
)

(sitePath <- sitePath(tree, align, 0.996))
mutations <- findFixed(sitePath)
```
