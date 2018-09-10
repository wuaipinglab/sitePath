#ifndef UTIL_H
#define UTIL_H

#include <map>
#include <string>
#include <utility>
#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;

class TipSeqLinker {
public:
  TipSeqLinker(CharacterVector sequence, IntegerVector tipPath);
  void proceed();
  const float compare(TipSeqLinker *linker);
  const int nextClade();
  const int currentClade();
  const int getTip();
  const int getRoot();
  const int getSeqLen();
  IntegerVector getPath();
private:
  std::string seq;
  IntegerVector path;
  const int tipIndex;
  int cIndex;
};

class TreeAlignmentMatch {
public:
  TreeAlignmentMatch(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    const float simThreshold
  );
protected:
  const float simCut;
  const int root, seqLen;
  std::vector<TipSeqLinker*> linkers;
  std::map< int, std::vector<TipSeqLinker*> > clusters;
private:
  std::map<std::pair<int, int>, float> compared;
  void pruneTree();
  const bool qualified(const std::vector<TipSeqLinker*> &clstr);
};

#endif