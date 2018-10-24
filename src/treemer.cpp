#include "pruner.h"

//[[Rcpp::export]]
const float compare(const std::string &query, const std::string &subject) {
  float match = 0.0, length = 0.0;
  for (
      std::string::const_iterator q = query.begin(), s = subject.begin();
      q != query.end(); ++q, ++s
  ) {
    if (*q != '-' && *s != '-') {
      length++;
      if (*q == *s) match++;
    }
  }
  return match / length;
}

// [[Rcpp::export]]
SEXP trimTree(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    const float &similarity, const bool &getTips
) {
  Pruner match(tipPaths, alignedSeqs, similarity);
  if (getTips) {
    return wrap(match.getTips());
  } else {
    return wrap(match.getPaths());
  }
}

// [[Rcpp::export]]
IntegerVector divergentNode(ListOf<IntegerVector> paths) {
  std::vector<int> res;
  for (int i = 0; i < paths.size() - 1; i++) {
    for (int j = i + 1; j < paths.size(); j++) {
      IntegerVector::iterator q  = paths[i].begin(), s = paths[j].begin();
      do {q++, s++;} while (*q == *s);
      if (--q != paths[i].begin()) res.push_back(*q);
    }
  }
  return wrap(res);
}

// [[Rcpp::export]]
IntegerVector getReference(std::string refSeq, const char &gapChar) {
  std::vector<int> res;
  for (int i = 0; i < refSeq.size(); i++) {
    if (refSeq[i] != gapChar) {
      res.push_back(i + 1);
    }
  }
  return wrap(res);
}

// [[Rcpp::export]]
ListOf<IntegerVector> ancestralPaths(ListOf<IntegerVector> paths, const int &n) {
  std::vector<IntegerVector> res;
  for (int i = 0; i < paths.size(); ++i) {
    if (paths[i].size() >= n) {
      res.push_back(paths[i][Range(0, n - 1)]);
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
