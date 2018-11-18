#include "pruner.h"

// [[Rcpp::export]]
NumericMatrix similarityMatrix(const ListOf<CharacterVector> &alignedSeqs) {
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
  std::map<std::pair<int, int>, float> simMatrix = std::map<std::pair<int, int>, float>();
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

/*
 * The following functions might be useful
 * regarding the match.getPath() output from trimTree() function
 * but are not relevant for now so commented out
 */

// // [[Rcpp::export]]
// ListOf<IntegerVector> terminalNode(ListOf<IntegerVector> paths) {
//   std::map< int, std::vector<int> > res;
//   for (int i = 0; i < paths.size(); i++) {
//     res[paths[i].end()[-2]].push_back(paths[i].end()[-1]);
//   }
//   return wrap(res);
// }
// 
// // [[Rcpp::export]]
// ListOf<IntegerVector> pathBeforeDivergence(ListOf<IntegerVector> paths) {
//   std::map< std::vector<int>, std::set<int> > fusedPaths;
//   for (int i = 0; i < paths.size() - 1; i++) {
//     for (int j = i + 1; j < paths.size(); j++) {
//       std::vector<int> path;
//       IntegerVector::iterator q  = paths[i].begin(), s = paths[j].begin();
//       do {path.push_back(*q); q++, s++;} while (*q == *s);
//       if (path.size() != 1) {
//         while (q != paths[i].end()) {fusedPaths[path].insert(*q); q++;}
//         while (s != paths[j].end()) {fusedPaths[path].insert(*s); s++;}
//       }
//     }
//   }
//   std::vector<IntegerVector> res;
//   for (std::map< std::vector<int>, std::set<int> >::iterator i = fusedPaths.begin(); i != fusedPaths.end(); i++) {
//     IntegerVector path = wrap(i->first);
//     path.attr("endNodes") = wrap(i->second);
//     path.attr("class") = "fused";
//     res.push_back(path);
//   }
//   return wrap(res);
// }
