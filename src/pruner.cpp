#include "pruner.h"

Pruner::Pruner(
  ListOf<IntegerVector> tipPaths, 
  ListOf<CharacterVector> alignedSeqs
): TreeAlignmentMatch(tipPaths, alignedSeqs) {}

std::map< int, std::vector<int> > Pruner::groupTips () {
  pruneTree();
  std::map< int, std::vector<int> > tipCluster;
  for (auto tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
    clusters[(*tsLinker)->currentClade()].push_back(*tsLinker);
  }
  for (auto it = clusters.begin(); it != clusters.end(); it++) {
    for (auto tsLinker = it->second.begin(); tsLinker != it->second.end(); tsLinker++) {
      tipCluster[it->first].push_back((*tsLinker)->getTip());
    }
  }
  return tipCluster;
}
