#ifndef SITEPATH_PRUNER_H
#define SITEPATH_PRUNER_H

#include <map>
#include <utility>
#include "util.h"

class TreeAlignmentMatch {
public:
  TreeAlignmentMatch(
    const ListOf<IntegerVector> &tipPaths, 
    const ListOf<CharacterVector> &alignedSeqs
  );
  virtual ~TreeAlignmentMatch() {}
  std::map< int, std::vector<int> > getTips() const;
  std::vector<IntegerVector> getPaths() const;
protected:
  std::vector<TipSeqLinker*> linkers;
  std::map< int, std::vector<TipSeqLinker*> > clusters;
  std::map< int, std::vector<TipSeqLinker*> > trueCluster;
  void pruneTree();
  virtual const bool qualified(
      const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate
  ) = 0;
private:
  const int root, seqLen;
};

class Pruner: public TreeAlignmentMatch {
public:
  Pruner(
    const ListOf<IntegerVector> &tipPaths, 
    const ListOf<CharacterVector> &alignedSeqs,
    const float simThreshold,
    std::map<std::pair<int, int>, float> &simMatrix
  );
private:
  const float simCut;
  std::map<std::pair<int, int>, float> compared;
  const bool qualified(
      const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate
  );
};

class CustomizablePruner: public TreeAlignmentMatch {
public:
  CustomizablePruner(
    const ListOf<IntegerVector> &tipPaths, 
    const ListOf<CharacterVector> &alignedSeqs,
    const std::map< int, std::vector<int> > &treeEdge,
    std::map<std::pair<int, int>, float> &simMatrix,
    const Function &customQualifyFunc
  );
private:
  std::map< int, std::vector<int> > nodeLink;
  std::map<std::pair<int, int>, float> compared;
  Function qualifyFunc;
  const bool qualified(
      const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate
  );
};

#endif // SITEPATH_PRUNER_H