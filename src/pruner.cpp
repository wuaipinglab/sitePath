#include "pruner.h"

TreeAlignmentMatch::TreeAlignmentMatch(
  ListOf<IntegerVector> tipPaths, 
  ListOf<CharacterVector> alignedSeqs
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
      throw std::invalid_argument("Sequene length not equal");
    }
  }
}

void TreeAlignmentMatch::pruneTree() {
  std::map< int, std::vector<TipSeqLinker*> > oldCluster;
  while (true) {
    oldCluster = clusters;
    clusters.clear();
    for (std::vector<TipSeqLinker*>::iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
      clusters[(*tsLinker)->nextClade()].push_back(*tsLinker);
    }
    if (clusters.size() == oldCluster.size()) {
      clusters.clear();
      break;
    }
    for (std::map< int, std::vector<TipSeqLinker*> >::iterator it = clusters.begin(); it != clusters.end(); it++) {
      bool exist = false;
      for (std::map< int, std::vector<TipSeqLinker*> >::iterator it2 = oldCluster.begin(); it2 != oldCluster.end(); it2++) {
        if (it->second == it2->second) {
          exist = true;
          oldCluster.erase(it2);
          break;
        }
      }
      if (!exist && qualified(it->second)) {
        for (std::vector<TipSeqLinker*>::iterator tsLinker = it->second.begin(); tsLinker != it->second.end(); tsLinker++) {
          (*tsLinker)->proceed();
        }
      }
    }
  }
}

Pruner::Pruner(
  ListOf<IntegerVector> tipPaths, 
  ListOf<CharacterVector> alignedSeqs,
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

const bool Pruner::qualified(const std::vector<TipSeqLinker*> &clstr) {
  for (std::vector<TipSeqLinker*>::const_iterator tsLinker = clstr.begin(); tsLinker != clstr.end() - 1; tsLinker++) {
    for (std::vector<TipSeqLinker*>::const_iterator tsLinker2 = tsLinker + 1; tsLinker2 != clstr.end(); tsLinker2++) {
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
