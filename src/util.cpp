#include "util.h"

TipSeqLinker::TipSeqLinker(
  CharacterVector sequence,
  IntegerVector tipPath
): seq(sequence), path(tipPath), maxIndex(tipPath.size() - 1) {
  pIndex = 0;
}

void TipSeqLinker::proceed() {
  if (maxIndex > pIndex) {
    pIndex++;
  }
};

const float TipSeqLinker::compare(TipSeqLinker *linker) {
  float match = 0, length = 0;
  for (
      CharacterVector::iterator q = seq.begin(), s = linker->seq.begin();
      q != seq.end(); q++, s++
  ) {
    if (*q == *s && *q != '-') {match++;}
    length++;
  }
  return match/length;
};

const std::vector<std::string> TipSeqLinker::siteComp(const std::vector<int> &sites) {
  std::vector<std::string> comp;
  for (
      std::vector<int>::const_iterator pos = sites.begin(); 
      pos != sites.end(); pos++
  ) {
    comp.push_back(as<std::string>(seq[*pos]));
  }
  return comp;
}

const int TipSeqLinker::nextClade() {
  if (maxIndex > pIndex) {
    return path[pIndex + 1];
  } else {
    return path[pIndex];
  }
}

const int TipSeqLinker::currentClade() {
  return path[pIndex];
}

const int TipSeqLinker::getTip() {
  return path[0];
}
const int TipSeqLinker::getRoot() {
  return path[maxIndex];
}
const int TipSeqLinker::getSeqLen() {
  return seq.size();
}

std::deque<int> TipSeqLinker::getPath() {
  std::deque<int> tPath;
  for (IntegerVector::iterator node = path.begin() + pIndex; node != path.end(); node++) {
    tPath.push_front(*node);
  }
  return tPath;
}

TreeAlignmentMatch::TreeAlignmentMatch(
  ListOf<IntegerVector> tipPaths, 
  ListOf<CharacterVector> alignedSeqs
): root(*(tipPaths[0].end() - 1)), seqLen(alignedSeqs[0].size()) {
  TipSeqLinker *linker;
  simCut = 0.9;
  pruned = false;
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

void TreeAlignmentMatch::setThreshold(const float sim = 0.9) {
  if (sim <= 0) {
    throw std::invalid_argument("Similarity cannot be lower or equal to 0");
  } else if (sim > 1) {
    throw std::invalid_argument("Similarity cannot be greater than 1");
  } else {
    simCut = sim;
  }
}

void TreeAlignmentMatch::setSites(IntegerVector rSites = IntegerVector::create()) {
  sites.clear();
  for (IntegerVector::iterator site = rSites.begin(); site != rSites.end(); site++) {
    if ((*site) <= 0) {
      throw std::invalid_argument("Site can't be negative");
    } else if ((*site) > seqLen) {
      throw std::invalid_argument("Site out of index");
    } else {
      sites.push_back((*site) - 1);
    }
  }
}

void TreeAlignmentMatch::pruneTree() {
  if (pruned) {
    throw std::runtime_error("A tree can't be pruned twice");
  }
  bool exist;
  std::map< int, std::vector<TipSeqLinker*> > oldCluster;
  while (true) {
    oldCluster = clusters;
    clusters.clear();
    for (std::vector<TipSeqLinker*>::iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
      clusters[(*tsLinker)->nextClade()].push_back(*tsLinker);
    }
    if (clusters.size() == oldCluster.size()) {
      clusters.clear();
      pruned = true;
      break;
    }
    for (std::map< int, std::vector<TipSeqLinker*> >::iterator it = clusters.begin(); it != clusters.end(); it++) {
      exist = false;
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
  std::vector<std::string> qComp, sComp;
  std::pair<int, int> pairing;
  for (std::vector<TipSeqLinker*>::const_iterator tsLinker = clstr.begin(); tsLinker != clstr.end(); tsLinker++) {
    qComp = (*tsLinker)->siteComp(sites);
    for (std::vector<TipSeqLinker*>::const_iterator tsLinker2 = tsLinker + 1; tsLinker2 != clstr.end(); tsLinker2++) {
      sComp = (*tsLinker2)->siteComp(sites);
      pairing = std::make_pair((*tsLinker)->getTip(), (*tsLinker2)->getTip());
      if (compared.find(pairing) != compared.end()) {
        sim = compared[pairing];
      } else {
        sim = (*tsLinker)->compare(*tsLinker2);
        compared[pairing] = sim;
      }
      if (sim < simCut or qComp != sComp) {
        return false;
      }
    }
  }
  return true;
}
