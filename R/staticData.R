#' @useDynLib sitePath
#' @importFrom Rcpp sourceCpp
NULL

#' @rdname sars2_align
#' @title Multiple sequence alignment of SARS-CoV-2 genome CDS
#' @description The raw sequences were downloaded from GISAID database
#'   (\url{https://www.gisaid.org/}) and aligned using MAFFT
#'   (\url{https://mafft.cbrc.jp/alignment/software/}) with default settings.
#' @format an \code{alignment} object
#' @usage data(sars2_align)
#' @docType data
"sars2_align"

#' @rdname sars2_tree
#' @title Phylogenetic tree of SARS-CoV-2 genome CDS
#' @description Tree was built from \code{\link{sars2_align}} using RAxML
#'   (\url{http://www.exelixis-lab.org/}) with default settings. The tip
#'   \code{EPI_ISL_402125} was used as the outgroup to root the tree.
#' @format a \code{phylo} object
#' @usage data(sars2_tree)
#' @docType data
"sars2_tree"

#' @rdname zikv_align
#' @title Multiple sequence alignment of Zika virus polyprotein
#' @description The raw protein sequences were downloaded from ViPR database
#'   (\url{https://www.viprbrc.org/}) and aligned using MAFFT
#'   (\url{https://mafft.cbrc.jp/alignment/software/}). with default settings.
#' @format an \code{alignment} object
#' @usage data(zikv_align)
#' @docType data
"zikv_align"

#' @rdname zikv_tree
#' @title Phylogenetic tree of Zika virus polyprotein
#' @description Tree was built from \code{\link{zikv_align}} using RAxML
#'   (\url{http://www.exelixis-lab.org/}) with default settings. The tip
#'   ANK57896 was used as outgroup to root the tree.
#' @format a \code{phylo} object
#' @usage data(zikv_tree)
#' @docType data
"zikv_tree"

#' @rdname zikv_align
#' @description \code{zikv_align_reduced} is a truncated version of
#'   \code{zikv_align}
#' @format an \code{alignment} object
#' @usage data(zikv_align_reduced)
#' @docType data
"zikv_align_reduced"

#' @rdname zikv_tree
#' @description \code{zikv_tree_reduced} is a truncated version of
#'   \code{zikv_tree}
#' @format a \code{phylo} object
#' @usage data(zikv_tree_reduced)
#' @docType data
"zikv_tree_reduced"

#' @rdname h3n2_align
#' @title Multiple sequence alignment of H3N2's HA protein
#' @description The raw protein sequences were downloaded from NCBI database and
#'   aligned using MAFFT (\url{https://mafft.cbrc.jp/alignment/software/}).
#' @format an \code{alignment} object
#' @usage data(h3n2_align)
#' @docType data
"h3n2_align"

#' @rdname h3n2_tree
#' @title Phylogenetic tree of H3N2's HA protein
#' @description Tree was built from \code{\link{h3n2_align}} using RAxML
#'   (\url{http://www.exelixis-lab.org/}) with default settings.
#' @format a \code{phylo} object
#' @usage data(h3n2_tree)
#' @docType data
"h3n2_tree"

#' @rdname h3n2_align
#' @description \code{h3n2_align_reduced} is a truncated version of
#'   \code{h3n2_align}
#' @format an \code{alignment} object
#' @usage data(h3n2_align_reduced)
#' @docType data
"h3n2_align_reduced"

#' @rdname h3n2_tree
#' @description \code{h3n2_tree_reduced} is a truncated version of
#'   \code{h3n2_tree}
#' @format a \code{phylo} object
#' @usage data(h3n2_tree_reduced)
#' @docType data
"h3n2_tree_reduced"

AA_COLORS <- c(
    "His" = "#8282D2",
    "Arg" = "#9370DB",
    "Lys" = "#145AFF",
    "Ile" = "#55AE3A",
    "Phe" = "#3232AA",
    "Leu" = "#0F820F",
    "Trp" = "#B45AB4",
    "Ala" = "#C8C8C8",
    "Met" = "#FFD700",
    "Pro" = "#DC9682",
    "Val" = "#2F4F2F",
    "Asn" = "#00DCDC",
    "Cys" = "#E6E600",
    "Gly" = "#666666",
    "Ser" = "#FF6347",
    "Tyr" = "#ADD8E6",
    "Gln" = "#0099CC",
    "Thr" = "#FA9600",
    "Glu" = "#8C1717",
    "Asp" = "#E60A0A",
    "gap" = "#d3d3d3",
    "unknown" = "#000000",
    "Ile_or_Leu" = "#000000",
    "Asp_or_Asn" = "#000000",
    "Glu_or_Gln" = "#000000"
)

AA_FULL_NAMES <- c(
    "h" = "His",
    "r" = "Arg",
    "k" = "Lys",
    "i" = "Ile",
    "f" = "Phe",
    "l" = "Leu",
    "w" = "Trp",
    "a" = "Ala",
    "m" = "Met",
    "p" = "Pro",
    "v" = "Val",
    "n" = "Asn",
    "c" = "Cys",
    "g" = "Gly",
    "s" = "Ser",
    "y" = "Tyr",
    "q" = "Gln",
    "t" = "Thr",
    "e" = "Glu",
    "d" = "Asp",
    "-" = "gap",
    "x" = "unknown",
    "j" = "Ile_or_Leu",
    "b" = "Asp_or_Asn",
    "z" = "Glu_or_Gln"
)

AA_SHORT_NAMES <- c(
    "His" = "H",
    "Arg" = "R",
    "Lys" = "K",
    "Ile" = "I",
    "Phe" = "F",
    "Leu" = "L",
    "Trp" = "W",
    "Ala" = "A",
    "Met" = "M",
    "Pro" = "P",
    "Val" = "V",
    "Asn" = "N",
    "Cys" = "C",
    "Gly" = "G",
    "Ser" = "S",
    "Tyr" = "Y",
    "Gln" = "Q",
    "Thr" = "T",
    "Glu" = "E",
    "Asp" = "D",
    "gap" = "-",
    "unknown" = "X",
    "Ile_or_Leu" = "J",
    "Asp_or_Asn" = "B",
    "Glu_or_Gln" = "Z"
)

AA_UNAMBIGUOUS <- c(
    "H",
    "R",
    "K",
    "I",
    "F",
    "L",
    "W",
    "A",
    "M",
    "P",
    "V",
    "N",
    "C",
    "G",
    "S",
    "Y",
    "Q",
    "T",
    "E",
    "D"
)

NT_COLORS <- c(
    "A" = "#5050ff",
    "C" = "#e00000",
    "G" = "#00c000",
    "T" = "#e6e600",
    "U" = "#cc9900",
    "-" = "#d3d3d3",
    "N" = "#000000",
    "X" = "#000000",
    "Y" = "#000000",
    "R" = "#000000",
    "W" = "#000000",
    "S" = "#000000",
    "K" = "#000000",
    "M" = "#000000",
    "D" = "#000000",
    "V" = "#000000",
    "H" = "#000000",
    "B" = "#000000"
)

NT_UNAMBIGUOUS <- c("A",
                    "C",
                    "G",
                    "T",
                    "U")
