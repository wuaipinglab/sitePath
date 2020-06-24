# sitePath: an R package for detection of site fixation in molecular evolution

A more detailed tutorial can be found
[here](https://bioconductor.org/packages/release/bioc/vignettes/sitePath/inst/doc/sitePathTutorial.html).

## 1\. Installation

[R programming language](https://cran.r-project.org/) \>= 3.6.0 is
required to use `sitePath`.

#### 1.1. Install from Bioconductor

The stable release is avaiable on
[Bioconductor](https://bioconductor.org/packages/release/bioc/html/sitePath.html).

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sitePath")
```

### 1.2. Install from GitHub

The package is in the experimental stage but gives the newest feature:

``` r
install.packages("devtools")
devtools::install_github("wuaipinglab/sitePath")
```

## 2\. A QuickStart

### 2.1. Import data

Both the phylogenetic tree and the matching multiple sequence alignment
are required.

``` r
library(ape)
library(sitePath)

# The file names of your tree and multiple sequence alignment
tree_file <- system.file("extdata", "ZIKV.newick", package = "sitePath")
alignment_file <- system.file("extdata", "ZIKV.fasta", package = "sitePath")

tree <- read.tree(tree_file)
tree <- addMSA(tree, alignment_file, "fasta")
```

### 2.2. Detect fixation sites

Use the `lineagePath` function to resolve major lineages (the choice of
threshold really depends). Then use the `fixationSites` function to
detect fixation sites.

``` r
paths <- lineagePath(tree, 0.05)

fixations <- fixationSites(paths)
print(fixations)
```

    ## Result for 3 paths:
    ## 
    ## 109 139 894 988 2074 2086 2634 3045 3144 107 1118 3353 1143 2842 3398 
    ## No reference sequence specified. Using alignment numbering
