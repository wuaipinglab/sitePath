#include "pruner.h"

// [[Rcpp::export]]
NumericMatrix getSimilarityMatrix(const ListOf<CharacterVector> &alignedSeqs) {
  int dim = alignedSeqs.size();
  NumericMatrix simMatrix(dim, dim);
  for (int i = 0; i < dim; ++i) {
    for (int j = i; j < dim; ++j) {
      if (i == j) {
        simMatrix(i, j) = 1;
      } else {
        simMatrix(j, i) = simMatrix(i, j) = compare(
          as<std::string>(alignedSeqs[i]),
          as<std::string>(alignedSeqs[j])
        );
      }
    }
  }
  return simMatrix;
}

// [[Rcpp::export]]
SEXP trimTree(
    const ListOf<IntegerVector> &tipPaths, 
    const ListOf<CharacterVector> &alignedSeqs,
    NumericMatrix &simMatrixInput,
    const float similarity, const bool getTips
) {
  std::map<std::pair<int, int>, float> simMatrix;
  int nrow = simMatrixInput.nrow();
  int ncol = simMatrixInput.ncol();
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (!R_IsNA(simMatrixInput(i, j))) {
        simMatrix[std::make_pair(i + 1, j + 1)] = simMatrixInput(i, j);
      }
    }
  }
  Pruner match(tipPaths, alignedSeqs, similarity, simMatrix);
  if (getTips) {
    return wrap(match.getTips());
  } else {
    return wrap(match.getPaths());
  }
}

// [[Rcpp::export]]
SEXP customTrimTree(
    const ListOf<IntegerVector> &tipPaths, 
    const ListOf<CharacterVector> &alignedSeqs,
    NumericMatrix &simMatrixInput,
    const NumericMatrix &treeEdge,
    const Function &customQualifyFunc,
    const bool getTips
) {
  std::map<std::pair<int, int>, float> simMatrix;
  int nrow = simMatrixInput.nrow();
  int ncol = simMatrixInput.ncol();
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (!R_IsNA(simMatrixInput(i, j))) {
        simMatrix[std::make_pair(i + 1, j + 1)] = simMatrixInput(i, j);
      }
    }
  }
  std::map< int, std::vector<int> > nodeLink;
  for (int i = 0; i < treeEdge.nrow(); ++i) {
    nodeLink[treeEdge(i, 0)].push_back(treeEdge(i, 1));
  }
  CustomizablePruner match(tipPaths, alignedSeqs, nodeLink, simMatrix, customQualifyFunc);
  if (getTips) {
    return wrap(match.getTips());
  } else {
    return wrap(match.getPaths());
  }
}

// [[Rcpp::export]]
IntegerVector divergentNode(const ListOf<IntegerVector> &paths) {
  std::vector<int> res;
  for (int i = 0; i < paths.size() - 1; i++) {
    for (int j = i + 1; j < paths.size(); j++) {
      IntegerVector::const_iterator q  = paths[i].begin(), s = paths[j].begin();
      do {q++, s++;} while (*q == *s);
      if (--q != paths[i].begin()) res.push_back(*q);
    }
  }
  return wrap(res);
}

// [[Rcpp::export]]
IntegerVector getReference(const std::string &refSeq, const char gapChar) {
  std::vector<int> res;
  for (unsigned int i = 0; i < refSeq.size(); i++) {
    if (refSeq[i] != gapChar) {
      res.push_back(i + 1);
    }
  }
  return wrap(res);
}

// [[Rcpp::export]]
ListOf<IntegerVector> ancestralPaths(const ListOf<IntegerVector> &paths, const int minLen) {
  std::vector<IntegerVector> res;
  for (int i = 0; i < paths.size(); ++i) {
    if (paths[i].size() >= minLen) {
      res.push_back(paths[i][Range(0, minLen - 1)]);
    }
  }
  return wrap(res);
}

// [[Rcpp::export]]
CharacterVector summarizeAA(
    const CharacterVector &seqs, 
    const int siteIndex, 
    const float tolerance
) {
  int nseq = seqs.size();
  std::map<char, int> aaSummary;
  for (int i = 0; i < nseq; ++i) {
    aaSummary[seqs[i][siteIndex]]++;
  }
  int currentMax = 0;
  char maxArg;
  for (std::map<char, int>::iterator i = aaSummary.begin(); i != aaSummary.end(); ++i) {
    if (i->second > currentMax) {
      currentMax = i->second;
      maxArg = i->first;
    }
  }
  int tt;
  if (tolerance < 0) {
    std::invalid_argument("tolerance can only be positive numeric");
  } else if (tolerance < 0.5) {
    tt = tolerance * nseq;
  } else if (tolerance < 1) {
    tt = nseq - tolerance * nseq;
  } else {
    tt = tolerance;
  }
  nseq -= currentMax;
  if (nseq > tt) {
    return NA_STRING;
  } else {
    return wrap(maxArg);
  }
}

// [[Rcpp::export]]
CharacterVector tip2colorEdge(
    CharacterVector &colorEdge,
    const std::string &color,
    const IntegerMatrix &treeEdge,
    const IntegerVector &tips,
    const int rootNode
) {
  std::map< int, std::pair<int, int> > nodeLink;
  for (int i = 0; i < treeEdge.nrow(); ++i) {
    nodeLink[treeEdge(i, 1)].first = treeEdge(i, 0);
    nodeLink[treeEdge(i, 1)].second = i;
  }
  for (int i = 0; i < tips.size(); ++i) {
    int cn, an;
    cn = tips[i];
    do {
      an = nodeLink[cn].first;
      colorEdge[nodeLink[cn].second] = color;
      cn = an;
    } while (an != rootNode);
  }
  return colorEdge;
}
