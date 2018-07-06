#include "pruner.h"
#include "siteExplorer.h"

// [[Rcpp::export]]
ListOf<IntegerVector> trimTree(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    NumericVector similarity
) {
  Pruner match(tipPaths, alignedSeqs);
  match.setThreshold(as<float>(similarity));
  return wrap(match.groupTips());
}

// [[Rcpp::export]]
ListOf< ListOf<CharacterVector> > mutationPath(
    ListOf<IntegerVector> tipPaths,
    ListOf<CharacterVector> alignedSeqsAR,
    NumericVector similarity,
    IntegerVector siteMode
) {
  SiteExplorer match(tipPaths, alignedSeqsAR);
  match.setThreshold(as<float>(similarity));
  std::map< std::string, std::map< std::string, std::set<std::string> > > res;
  for (IntegerVector::iterator m = siteMode.begin(); m != siteMode.end(); ++m) {
    switch (*m) {
    case 1: res["fixed"] = match.getSitePath(0);
      break;
    case 2: res["alternate"] = match.getSitePath(1);
      break;
    case 3: res["coevolve"] = match.getSitePath(2);
      break;
    default: throw std::invalid_argument("Invalid siteMode argument");
    }
  }
  return wrap(res);
}
