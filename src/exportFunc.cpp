#include "treemer.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix getSimilarityMatrix(
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
) {
    int dim = alignedSeqs.size();
    Rcpp::NumericMatrix simMatrix(dim, dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = i; j < dim; ++j) {
            if (i == j) { simMatrix(i, j) = 1; } else {
                simMatrix(j, i) = simMatrix(i, j) = compare(
                    Rcpp::as<std::string>(alignedSeqs[i]),
                    Rcpp::as<std::string>(alignedSeqs[j])
                );
            }
        }
    }
    return simMatrix;
}

// [[Rcpp::export]]
SEXP runTreemer(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        Rcpp::NumericMatrix &simMatrixInput,
        const float similarity, const bool getTips
) {
    std::map<std::pair<int, int>, float> simMatrix;
    int nrow = simMatrixInput.nrow(), ncol = simMatrixInput.ncol();
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            if (!R_IsNA(simMatrixInput(i, j))) {
                simMatrix[std::make_pair(i + 1, j + 1)] = simMatrixInput(i, j);
            }
        }
    }
    Treemer::BySimilarity match(tipPaths, alignedSeqs, similarity, simMatrix);
    if (getTips) {
        return Rcpp::wrap(match.getTips());
    } else {
        return Rcpp::wrap(match.getPaths());
    }
}

// // [[Rcpp::export]]
// SEXP customTrimTree(
//         const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
//         const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
//         Rcpp::NumericMatrix &simMatrixInput,
//         const Rcpp::NumericMatrix &treeEdge,
//         const Rcpp::Function &customQualifyFunc,
//         const bool getTips
// ) {
//     std::map<std::pair<int, int>, float> simMatrix;
//     int nrow = simMatrixInput.nrow(), ncol = simMatrixInput.ncol();
//     for (int i = 0; i < nrow; ++i) {
//         for (int j = 0; j < ncol; ++j) {
//             if (!R_IsNA(simMatrixInput(i, j))) {
//                 simMatrix[std::make_pair(i + 1, j + 1)] = simMatrixInput(i, j);
//             }
//         }
//     }
//     std::map< int, std::vector<int> > nodeLink;
//     for (int i = 0; i < treeEdge.nrow(); ++i) {
//         nodeLink[treeEdge(i, 0)].push_back(treeEdge(i, 1));
//     }
//     CustomizablePruner match(
//             tipPaths, alignedSeqs, nodeLink,
//             simMatrix, customQualifyFunc
//     );
//     if (getTips) {
//         return Rcpp::wrap(match.getTips());
//     } else {
//         return Rcpp::wrap(match.getPaths());
//     }
// }

// [[Rcpp::export]]
Rcpp::IntegerVector divergentNode(
        const Rcpp::ListOf<Rcpp::IntegerVector> &paths
) {
    std::vector<int> res;
    for (int i = 0; i < paths.size() - 1; i++) {
        for (int j = i + 1; j < paths.size(); j++) {
            Rcpp::IntegerVector::const_iterator 
            q  = paths[i].begin(), s = paths[j].begin();
            do { q++, s++; } while (*q == *s);
            if (--q != paths[i].begin()) { res.push_back(*q); }
        }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::IntegerVector getReference(
        const std::string &refSeq, 
        const char gapChar
) {
    std::vector<int> res;
    for (unsigned int i = 0; i < refSeq.size(); i++) {
        if (refSeq[i] != gapChar) { res.push_back(i + 1); }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> ancestralPaths(
        const Rcpp::ListOf<Rcpp::IntegerVector> &paths,
        const int minLen
) {
    std::vector<Rcpp::IntegerVector> res;
    for (int i = 0; i < paths.size(); ++i) {
        if (paths[i].size() >= minLen) { 
            res.push_back(paths[i][Rcpp::Range(0, minLen - 1)]);
        }
    }
    return wrap(res);
}

// [[Rcpp::export]]
Rcpp::CharacterVector summarizeAA(
        const Rcpp::CharacterVector &seqs, 
        const int siteIndex, 
        const float tolerance
) {
    int nseq = seqs.size();
    std::map<char, int> aaSummary;
    for (int i = 0; i < nseq; ++i) { aaSummary[seqs[i][siteIndex]]++; }
    auto it = aaSummary.begin();
    int currentMax = it->second;
    char maxArg = it->first;
    for (++it; it != aaSummary.end(); ++it) {
        if (it->second > currentMax) {
            currentMax = it->second;
            maxArg = it->first;
        }
    }
    // the default tolerance would be a number
    int tt = static_cast<int>(tolerance);
    // get the tolerance number if receiving float
    if (tolerance < 0.5) {
        tt = tolerance * nseq;
    } else if (tolerance < 1) {
        tt = nseq - tolerance * nseq;
    }
    // test if the number of non-dominant AAs within tolerance
    nseq -= currentMax;
    if (nseq > tt) {
        return NA_STRING;
    } else {
        Rcpp::CharacterVector res = Rcpp::wrap(maxArg);
        res.attr("n") = Rcpp::IntegerVector::create(nseq);
        return res;
    }
}

// [[Rcpp::export]]
Rcpp::CharacterVector tip2colorEdge(
        Rcpp::CharacterVector &colorEdge,
        const std::string &color,
        const Rcpp::IntegerMatrix &treeEdge,
        const Rcpp::IntegerVector &tips,
        const int rootNode
) {
    std::map< int, std::pair<int, int> > nodeLink;
    for (int i = 0; i < treeEdge.nrow(); ++i) {
        nodeLink[treeEdge(i, 1)].first = treeEdge(i, 0);
        nodeLink[treeEdge(i, 1)].second = i;
    }
    for (int i = 0; i < tips.size(); ++i) {
        int cn = tips[i], an;
        do {
            an = nodeLink[cn].first;
            colorEdge[nodeLink[cn].second] = color;
            cn = an;
        } while (an != rootNode);
    }
    return colorEdge;
}