#include "pruner.h"

Pruner::Pruner(
  ListOf<IntegerVector> tipPaths, 
  ListOf<CharacterVector> alignedSeqs,
  const float simThreshold
): TreeAlignmentMatch(tipPaths, alignedSeqs, simThreshold) {}

std::map< int, std::vector<int> > Pruner::getTips() {
  std::map< int, std::vector<int> > tipCluster;
  for (std::vector<TipSeqLinker*>::iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
    tipCluster[(*tsLinker)->currentClade()].push_back((*tsLinker)->getTip());
  }
  return tipCluster;
}

std::vector<IntegerVector> Pruner::getPaths() {
  std::vector<IntegerVector> paths;
  for (std::vector<TipSeqLinker*>::iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
    paths.push_back((*tsLinker)->getPath());
  }
  return paths;
}