#include "pruner.h"

Pruner::Pruner(
  ListOf<IntegerVector> tipPaths, 
  ListOf<CharacterVector> alignedSeqs
): TreeAlignmentMatch(tipPaths, alignedSeqs) {}

map< string, vector<int> > Pruner::groupTips () {
  pruneTree();
  map< string, vector<int> > tipCluster;
  for (vector<TipSeqLinker*>::iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
    clusters[(*tsLinker)->currentClade()].push_back(*tsLinker);
  }
  for (map< int, vector<TipSeqLinker*> >::iterator it = clusters.begin(); it != clusters.end(); it++) {
    for (vector<TipSeqLinker*>::iterator tsLinker = it->second.begin(); tsLinker != it->second.end(); tsLinker++) {
      tipCluster[to_string(it->first)].push_back((*tsLinker)->getTip());
    }
  }
  return tipCluster;
}
