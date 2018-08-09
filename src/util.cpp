#include "util.h"

TipSeqLinker::TipSeqLinker(
  CharacterVector sequence,
  IntegerVector tipPath
):
  seq(as<std::string>(sequence)),
  path(tipPath),
  tipIndex(tipPath.size() - 1),
  cIndex(tipIndex) {}

void TipSeqLinker::proceed() {
  if (cIndex > 0) {cIndex--;}
};

const float TipSeqLinker::compare(TipSeqLinker *linker) {
  float match = 0, length = 0;
  for (
      std::string::iterator q = seq.begin(), s = linker->seq.begin();
      q != seq.end(); q++, s++
  ) {
    if (*q != '-' && *s != '-') {
      length++;
      if (*q == *s) {match++;}
    }
  }
  return match/length;
}

const int TipSeqLinker::currentClade() {
  return path[cIndex];
}

const int TipSeqLinker::nextClade() {
  if (cIndex > 0) {
    return path[cIndex - 1];
  } else {
    return currentClade();
  }
}

const int TipSeqLinker::getTip() {
  return path[tipIndex];
}

const int TipSeqLinker::getRoot() {
  return path[0];
}

const int TipSeqLinker::getSeqLen() {
  return seq.size();
}

IntegerVector TipSeqLinker::getPath() {
  return path[Range(0, cIndex)];
}

TreeAlignmentMatch::TreeAlignmentMatch(
  ListOf<IntegerVector> tipPaths, 
  ListOf<CharacterVector> alignedSeqs,
  const float simThreshold
):
  root(*(tipPaths[0].begin())),
  seqLen((as<std::string>(alignedSeqs[0])).size()),
  simCut(simThreshold) {
  if (simCut <= 0) {
    throw std::invalid_argument("Similarity cannot be lower or equal to 0");
  } else if (simCut > 1) {
    throw std::invalid_argument("Similarity cannot be greater than 1");
  }
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
  if (simCut != 1) {pruneTree();}
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

const bool TreeAlignmentMatch::qualified(const std::vector<TipSeqLinker*> &clstr) {
  float sim;
  std::pair<int, int> pairing;
  for (std::vector<TipSeqLinker*>::const_iterator tsLinker = clstr.begin(); tsLinker != clstr.end() - 1; tsLinker++) {
    for (std::vector<TipSeqLinker*>::const_iterator tsLinker2 = tsLinker + 1; tsLinker2 != clstr.end(); tsLinker2++) {
      pairing = std::make_pair((*tsLinker)->getTip(), (*tsLinker2)->getTip());
      if (compared.find(pairing) != compared.end()) {
        sim = compared[pairing];
      } else {
        sim = (*tsLinker)->compare(*tsLinker2);
        compared[pairing] = sim;
      }
      if (sim < simCut) {
        return false;
      }
    }
  }
  return true;
}
