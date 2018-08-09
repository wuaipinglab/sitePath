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
ListOf<IntegerVector> divergentNode(ListOf<IntegerVector> paths) {
  std::vector<IntegerVector> divPoints;
  for (int i = 0; i < paths.size() - 1; i++) {
    IntegerVector query = paths[i];
    for (int j = i + 1; j < paths.size(); j++) {
      IntegerVector subject = paths[j];
      IntegerVector::iterator q, s;
      for (q = query.begin(), s = subject.begin(); *q == *s; q++, s++) {
        continue;
      }
      if (q - 1 != query.begin()) {
        divPoints.push_back(IntegerVector(query.begin(), q));
      }
    }
  }
  return wrap(divPoints);
}
