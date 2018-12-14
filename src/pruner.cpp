#include "pruner.h"

TreeAlignmentMatch::TreeAlignmentMatch(
  const ListOf<IntegerVector> &tipPaths,
  const ListOf<CharacterVector> &alignedSeqs
):
  root(*(tipPaths[0].begin())),
  seqLen((as<std::string>(alignedSeqs[0])).size()) {
  TipSeqLinker *linker;
  for (int i = 0; i < tipPaths.size(); i++) {
    linker = new TipSeqLinker(alignedSeqs[i], tipPaths[i]);
    linkers.push_back(linker);
    clusters[linker->getTip()].push_back(linker);
    if (linker->getRoot() != root) {
      throw std::invalid_argument("Root in tree path not equal");
    } else if (linker->getSeqLen() != seqLen) {
      throw std::invalid_argument("Sequence length not equal");
    }
  }
}

std::map< int, std::vector<int> > TreeAlignmentMatch::getTips() const {
  std::map< int, std::vector<int> > tipCluster;
  for (std::vector<TipSeqLinker*>::const_iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
    tipCluster[(*tsLinker)->currentClade()].push_back((*tsLinker)->getTip());
  }
  return tipCluster;
}

std::vector<IntegerVector> TreeAlignmentMatch::getPaths() const {
  std::vector<IntegerVector> paths;
  for (std::vector<TipSeqLinker*>::const_iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
    paths.push_back((*tsLinker)->getPath());
  }
  return paths;
}

void TreeAlignmentMatch::pruneTree() {
  std::map< int, std::vector<TipSeqLinker*> > oldCluster;
  while (true) {
    for (std::vector<TipSeqLinker*>::const_iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
      trueCluster[(*tsLinker)->currentClade()].push_back(*tsLinker);
    }
    oldCluster = clusters;
    clusters.clear();
    // look down one more node (fake 'proceed') for each tip after new positioning
    for (std::vector<TipSeqLinker*>::iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
      clusters[(*tsLinker)->nextClade()].push_back(*tsLinker); //  group tips by fake 'proceed' node
    }
    // if no more group 'kissed' each other by a common ancestral node after fake 'proceed', then pruning is done
    if (clusters.size() == oldCluster.size()) {
      clusters.clear();
      break;
    }
    // only 'kissed' group can do real 'proceed'
    // if a grouping doesn't exist in 'oldCluster' then all tips in that group can 'proceed'
    for (std::map< int, std::vector<TipSeqLinker*> >::iterator it = clusters.begin(); it != clusters.end(); it++) {
      bool kissed = true; // assume a group is kissed with another (give it benefit of the doubt)
      for (std::map< int, std::vector<TipSeqLinker*> >::iterator it2 = oldCluster.begin(); it2 != oldCluster.end(); it2++) {
        // a group is 'non-kissed' after fake 'proceed' if it can be found in 'oldCluster'
        if (it->second == it2->second) {
          kissed = false;
          oldCluster.erase(it2); // a 'non-kissed' group won't appear twice in 'clusters' so deleted
          break;
        }
      }
      // candidate group needs to pass some requirement to be qualified 'kissed'
      if (kissed && qualified(it)) {
        for (std::vector<TipSeqLinker*>::iterator tsLinker = it->second.begin(); tsLinker != it->second.end(); tsLinker++) {
          (*tsLinker)->proceed();
        }
      }
    }
    trueCluster.clear();
    oldCluster.clear();
  }
}

Pruner::Pruner(
  const ListOf<IntegerVector> &tipPaths,
  const ListOf<CharacterVector> &alignedSeqs,
  const float simThreshold,
  std::map<std::pair<int, int>, float> &simMatrix
):
  TreeAlignmentMatch(tipPaths, alignedSeqs),
  simCut(simThreshold),
  compared(simMatrix) {
  if (simCut <= 0) {
    throw std::invalid_argument("Similarity cannot be lower or equal to 0");
  } else if (simCut > 1) {
    throw std::invalid_argument("Similarity cannot be greater than 1");
  }
  if (simCut != 1) {pruneTree();}
}

const bool Pruner::qualified(const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate) {
  for (std::vector<TipSeqLinker*>::const_iterator tsLinker = candidate->second.begin(); tsLinker != candidate->second.end() - 1; tsLinker++) {
    for (std::vector<TipSeqLinker*>::const_iterator tsLinker2 = tsLinker + 1; tsLinker2 != candidate->second.end(); tsLinker2++) {
      std::pair<int, int> pairing = std::make_pair((*tsLinker)->getTip(), (*tsLinker2)->getTip());
      float sim;
      if (compared.find(pairing) != compared.end()) {
        sim = compared[pairing];
      } else {
        sim = compare((*tsLinker)->getSeq(), (*tsLinker2)->getSeq());
        compared[pairing] = sim;
      }
      if (sim < simCut) {
        return false;
      }
    }
  }
  return true;
}

CustomizablePruner::CustomizablePruner(
  const ListOf<IntegerVector> &tipPaths,
  const ListOf<CharacterVector> &alignedSeqs,
  const std::map< int, std::vector<int> > &treeEdge,
  std::map<std::pair<int, int>, float> &simMatrix,
  const Function &customQualifyFunc
):
  TreeAlignmentMatch(tipPaths, alignedSeqs),
  nodeLink(treeEdge),
  compared(simMatrix),
  qualifyFunc(customQualifyFunc) {pruneTree();}

const bool CustomizablePruner::qualified(const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate) {
  std::vector<float> crossSim, combinedSim;
  int nChildren = nodeLink[candidate->first].size();
  for (int i = 0; i < nChildren; ++i) {
    std::vector<TipSeqLinker*> *oldClstr = &trueCluster[nodeLink[candidate->first][i]];
    if (oldClstr->size() == 1) {goto CROSSCOMPARE;}
    for (std::vector<TipSeqLinker*>::iterator it1 = oldClstr->begin(); it1 != oldClstr->end() - 1; ++it1) {
      for (std::vector<TipSeqLinker*>::iterator it2 = it1 + 1; it2 != oldClstr->end(); ++it2) {
        combinedSim.push_back(compared[std::make_pair((*it1)->getTip(), (*it2)->getTip())]);
      }
    }
    CROSSCOMPARE: for (int j = i + 1; j < nChildren; ++j) {
      std::vector<TipSeqLinker*> *oldClstr2 = &trueCluster[nodeLink[candidate->first][j]];
      for (std::vector<TipSeqLinker*>::iterator it1 = oldClstr->begin(); it1 != oldClstr->end(); ++it1) {
        for (std::vector<TipSeqLinker*>::iterator it2 = oldClstr2->begin(); it2 != oldClstr2->end(); ++it2) {
          crossSim.push_back(compared[std::make_pair((*it1)->getTip(), (*it2)->getTip())]);
        }
      }
    }
  }
  if (crossSim.empty() || combinedSim.empty()) {
    return true;
  }
  return as<bool>(qualifyFunc(wrap(crossSim), wrap(combinedSim)));
}
