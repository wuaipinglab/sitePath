#include "pruner.h"

// [[Rcpp::export]]
SEXP trimTree(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    float similarity, bool getTips
) {
  Pruner match(tipPaths, alignedSeqs, similarity);
  if (getTips) {
    return wrap(match.getTips());
  } else {
    return wrap(match.getPaths());
  }
}

// [[Rcpp::export]]
IntegerVector getReference(std::string refSeq, char gapChar) {
  std::vector<int> res;
  for (int i = 0; i < refSeq.size(); i++) {
    if (refSeq[i] != gapChar) {
      res.push_back(i + 1);
    }
  }
  return wrap(res);
}

/*
 * The following functions might be useful in the future
 * regarding the match.getPath() output from trimTree function
 * but are not relevant for now so commented out
 */

// // [[Rcpp::export]]
// IntegerVector divergentNode(ListOf<IntegerVector> paths) {
//   std::vector<int> res;
//   for (int i = 0; i < paths.size() - 1; i++) {
//     for (int j = i + 1; j < paths.size(); j++) {
//       IntegerVector::iterator q = paths[i].begin(), s = paths[j].begin();
//       while (*q == *s) {q++, s++;}
//       if (--q != paths[i].begin()) res.push_back(*q);
//     }
//   }
//   return wrap(res);
// }
// 
// // [[Rcpp::export]]
// ListOf<IntegerVector> terminalNode(ListOf<IntegerVector> paths) {
//   std::map< int, std::vector<int> > res;
//   for (int i = 0; i < paths.size(); i++) {
//     if (paths[i].size() > 1) {
//       res[paths[i].end()[-2]].push_back(paths[i].end()[-1]);
//     }
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
//       IntegerVector::iterator q = paths[i].begin(), s = paths[j].begin();
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
