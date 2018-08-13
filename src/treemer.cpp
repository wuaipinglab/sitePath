#include "pruner.h"

// [[Rcpp::export]]
SEXP trimTree(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    NumericVector similarity,
    LogicalVector getTips
) {
  Pruner match(tipPaths, alignedSeqs, as<float>(similarity));
  if (as<bool>(getTips)) {
    return wrap(match.getTips());
  } else {
    return wrap(match.getPaths());
  }
}

// [[Rcpp::export]]
IntegerVector terminalNode(ListOf<IntegerVector> paths) {
  std::map<int, int> res;
  for (int i = 0; i < paths.size(); i++) {
    if (res.find(paths[i].end()[-2]) == res.end()) {
      res[paths[i].end()[-2]] = paths[i].end()[-1];
    } else {
      res[paths[i].end()[-2]] = paths[i].end()[-2];
    }
  }
  return wrap(res);
}

// [[Rcpp::export]]
ListOf<IntegerVector> pathBeforeDivergence(ListOf<IntegerVector> paths) {
  std::vector<IntegerVector> res;
  for (int i = 0; i < paths.size() - 1; i++) {
    for (int j = i + 1; j < paths.size(); j++) {
      IntegerVector::iterator q  = paths[i].begin(), s = paths[j].begin();
      while (*q == *s) {q++, s++;}
      if (q - 1 != paths[i].begin()) {
        res.push_back(IntegerVector(paths[i].begin(), q));
      }
    }
  }
  return wrap(res);
}
