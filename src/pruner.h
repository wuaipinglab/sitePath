#ifndef SITEPATH_PRUNER_H
#define SITEPATH_PRUNER_H

#include <map>
#include <utility>
#include "util.h"

class TreeAlignmentMatch {
public:
  TreeAlignmentMatch(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs
  );
  virtual ~TreeAlignmentMatch() {}
protected:
  const int root, seqLen;
  std::vector<TipSeqLinker*> linkers;
  std::map< int, std::vector<TipSeqLinker*> > clusters;
  void pruneTree();
  virtual const bool qualified(const std::vector<TipSeqLinker*> &clstr) = 0;
};

class Pruner: public TreeAlignmentMatch {
public:
  Pruner(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    const float simThreshold,
    std::map<std::pair<int, int>, float> &simMatrix
  );
  std::map< int, std::vector<int> > getTips();
  std::vector<IntegerVector> getPaths();
private:
  const float simCut;
  std::map<std::pair<int, int>, float> compared;
  const bool qualified(const std::vector<TipSeqLinker*> &clstr);
};

#endif // SITEPATH_PRUNER_H