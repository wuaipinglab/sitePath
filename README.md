# sitePath: an R package for detection of site fixation in molecular evolution

## Getting help

Post on Bioconductor [support site](https://support.bioconductor.org/)
if having trouble using `sitePath`. Or open an
[issue](https://github.com/wuaipinglab/sitePath/issues/new?assignees=&labels=&template=bug_report.md&title=)
if a bug is found.

## Installation

[R programming language](https://cran.r-project.org/) \>= 3.6.0 is
required to use `sitePath`.

The stable release is available on
[Bioconductor](https://bioconductor.org/packages/release/bioc/html/sitePath.html).

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sitePath")
```

The installation from [GitHub](https://github.com/wuaipinglab/sitePath/)
is in experimental stage but gives the newest feature:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("wuaipinglab/sitePath")
```

## QuickStart

Both the phylogenetic tree and the matching multiple sequence alignment
files are required.

``` r
library(sitePath)

# The file names of your tree and multiple sequence alignment
tree_file <- system.file("extdata", "ZIKV.newick", package = "sitePath")
alignment_file <- system.file("extdata", "ZIKV.fasta", package = "sitePath")

tree <- read.tree(tree_file)
tree <- addMSA(tree, alignment_file, "fasta")
```

Use the `lineagePath` function to resolve major lineages (the choice of
threshold really depends). Then use the `fixationSites` function to
detect fixation sites.

``` r
paths <- lineagePath(tree)

fixations <- fixationSites(paths)
print(fixations)
```

    ## This is a 'fixationSites' object.
    ## 
    ## Result for 3 paths:
    ## 
    ## 109 139 894 988 2074 2086 2634 3045 3144 107 1118 3353 1143 2842 3398 
    ## No reference sequence specified. Using alignment numbering
