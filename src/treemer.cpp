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
  std::map< std::string, std::map< std::string, std::set<std::string> > > sitePath;
  for (IntegerVector::iterator m = siteMode.begin(); m != siteMode.end(); ++m) {
    switch (*m) {
    case 1: sitePath["fixed"] = match.getSitePath(0);
      break;
    case 2: sitePath["alternate"] = match.getSitePath(1);
      break;
    case 3: sitePath["coevolve"] = match.getSitePath(2);
      break;
    default: throw std::invalid_argument("Invalid siteMode argument");
    }
  }
  List res = wrap(sitePath);
  res.attr("evolPath") = wrap(match.getPath());
  res.attr("divPoints") = wrap(match.getDivPoints());
  return res;
}

// [[Rcpp::export]]
ListOf<CharacterVector> mutationList(
    ListOf<IntegerVector> tipPaths,
    ListOf<CharacterVector> alignedSeqsAR,
    NumericVector similarity
) {
  MutationExplorer match(tipPaths, alignedSeqsAR);
  match.setThreshold(as<float>(similarity));
  List res = wrap(match.getMutationList());
  res.attr("divPoints") = wrap(match.getDivPoints());
  return res;
}
