#include <Rcpp.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <deque>
#include <set>
#include <algorithm>
using namespace Rcpp;
using namespace std;

class TipSeqLinker {
public:
  TipSeqLinker(
    CharacterVector sequence, 
    IntegerVector tipPath
  ) {
    seq = sequence;
    path = tipPath;
    maxIndex = tipPath.size() - 1;
    pIndex = 0;
  }
  void proceed() {
    if (maxIndex > pIndex) {
      pIndex++;
    }
  };
  float compare(TipSeqLinker *linker) {
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
  vector<string> siteComp(vector<int> sites) {
    vector<string> comp;
    for (
        vector<int>::iterator pos = sites.begin(); 
        pos != sites.end(); pos++
    ) {
      comp.push_back(as<string>(seq[*pos]));
    }
    return comp;
  }
  int nextClade() {
    if (maxIndex > pIndex) {
      return path[pIndex + 1];
    } else {
      return path[pIndex];
    }
  }
  int currentClade() {
    return path[pIndex];
  }
  const int getTip() {
    return path[0];
  }
  const int getRoot() {
    return path[maxIndex];
  }
  const int getSeqLen() {
    return seq.size();
  }
  deque<int> getPath() {
    deque<int> tPath;
    for (
        IntegerVector::iterator node = path.begin() + pIndex; 
        node != path.end(); node++
    ) {
      tPath.push_front(*node);
    }
    return tPath;
  }
private:
  CharacterVector seq;
  IntegerVector path;
  int pIndex;
  int maxIndex;
};

class TreeAlignmentMatch {
public:
  TreeAlignmentMatch(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs
  ) {
    TipSeqLinker *linker;
    simCut = 0.9;
    pruned = false;
    for (int i = 0; i < tipPaths.size(); i++) {
      linker = new TipSeqLinker(
        alignedSeqs[i],
        tipPaths[i]
      );
      linkers.push_back(linker);
      clusters[linker->getTip()].push_back(linker);
      if (i == 0) {
        root = linker->getRoot();
        seqLen = linker->getSeqLen();
      } else if (linker->getRoot() != root) {
        throw invalid_argument("Root in tree path not equal");
      } else if (linker->getSeqLen() != seqLen) {
        throw invalid_argument("Sequene length not equal");
      }
    }
  }
  virtual ~TreeAlignmentMatch() {}
  void setThreshold(float sim = 0.9) {
    if (sim <= 0) {
      throw invalid_argument("Similarity cannot be lower or equal to 0");
    } else if (sim > 1) {
      throw invalid_argument("Similarity cannot be greater than 1");
    } else {
      simCut = sim;
    }
  }
  void setSites(IntegerVector rSites = IntegerVector::create()) {
    sites.clear();
    for (
        IntegerVector::iterator site = rSites.begin(); 
        site != rSites.end(); site++
    ) {
      if ((*site) <= 0) {
        throw invalid_argument("Site can't be negative");
      } else if ((*site) > seqLen) {
        throw invalid_argument("Site out of index");
      } else {
        sites.push_back((*site) - 1);
      }
    }
  }
protected:
  bool pruned;
  float simCut;
  int root, seqLen;
  vector<int> sites;
  vector<TipSeqLinker*> linkers;
  vector<TipSeqLinker*>::iterator tsLinker;
  map< int, vector<TipSeqLinker*> > clusters;
  map< int, vector<TipSeqLinker*> >::iterator it;
  void pruneTree() {
    if (pruned) {
      throw runtime_error("A tree can't be pruned twice");
    }
    bool exist;
    map< int, vector<TipSeqLinker*> > oldCluster;
    map< int, vector<TipSeqLinker*> >::iterator it2;
    while (true) {
      oldCluster = clusters;
      clusters.clear();
      for (tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
        clusters[(*tsLinker)->nextClade()].push_back(*tsLinker);
      }
      if (clusters.size() == oldCluster.size()) {
        clusters.clear();
        pruned = true;
        break;
      }
      for (it = clusters.begin(); it != clusters.end(); it++) {
        exist = false;
        for (it2 = oldCluster.begin(); it2 != oldCluster.end(); it2++) {
          if (it->second == it2->second) {
            exist = true;
            oldCluster.erase(it2);
            break;
          }
        }
        if (!exist && qualified(it->second)) {
          for (tsLinker = it->second.begin(); tsLinker != it->second.end(); tsLinker++) {
            (*tsLinker)->proceed();
          }
        }
      }
    }
  }
private:
  float sim;
  map<pair<int, int>, float> compared;
  vector<TipSeqLinker*>::iterator tsLinker2;
  bool qualified(vector<TipSeqLinker*> clstr) {
    vector<string> qComp, sComp;
    pair<int, int> pairing;
    for (tsLinker = clstr.begin(); tsLinker != clstr.end(); tsLinker++) {
      qComp = (*tsLinker)->siteComp(sites);
      for (tsLinker2 = tsLinker + 1; tsLinker2 != clstr.end(); tsLinker2++) {
        sComp = (*tsLinker2)->siteComp(sites);
        pairing = make_pair((*tsLinker)->getTip(), (*tsLinker2)->getTip());
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
};

class Pruner: public TreeAlignmentMatch {
public:
  Pruner(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs
  ): TreeAlignmentMatch(tipPaths, alignedSeqs) {}
  map< string, vector<int> > groupTips () {
    pruneTree();
    map< string, vector<int> > tipCluster;
    for (tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
      clusters[(*tsLinker)->currentClade()].push_back(*tsLinker);
    }
    for (it = clusters.begin(); it != clusters.end(); it++) {
      for (tsLinker = it->second.begin(); tsLinker != it->second.end(); tsLinker++) {
        tipCluster[to_string(it->first)].push_back((*tsLinker)->getTip());
      }
    }
    return tipCluster;
  }
};

class SiteExplorer: public TreeAlignmentMatch {
public:
  SiteExplorer(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs
  ): TreeAlignmentMatch(tipPaths, alignedSeqs) {
    ancestralSeqs = alignedSeqs;
    divPoints.insert(root);
  }
  void getEvolPath() {
    for (tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
      evolPath.push_back((*tsLinker)->getPath());
    }
    if (simCut != 1) {
      pruneTree();
      map<deque<int>, int> pn;
      for (tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
        pn[(*tsLinker)->getPath()]++;
      }
      if (pn.size() != linkers.size()) {
        evolPath.clear();
        for (map<deque<int>, int>::iterator vp = pn.begin(); vp != pn.end(); vp++) {
          if (vp->first.size() != 1 && vp->second > 1) {
            evolPath.push_back(vp->first);
          }
        }
      }
    }
    vector< deque<int> >::iterator ePath2;
    deque<int>::iterator cNode2;
    for (ePath = evolPath.begin(); ePath != evolPath.end(); ePath++) {
      for (ePath2 = ePath + 1; ePath2 != evolPath.end(); ePath2++) {
        for (
            cNode = (*ePath).begin(), cNode2 = (*ePath2).begin();
            (*cNode) == (*cNode2); cNode++, cNode2++
        ) continue;
        if (cNode != (*ePath).begin()) divPoints.insert(*(cNode - 1));
      }
      for (pos = 0; pos < seqLen; pos++) {
        for (cNode = (*ePath).begin(); cNode != (*ePath).end(); cNode++) {
          cSite = ancestralSeqs[(*cNode) - 1][pos];
          allele[pos].insert(cSite);
          if (!pSite.empty() && cSite != pSite) {
            mutNodes[pos].insert(make_pair(pNode, *cNode));
          }
          pSite = cSite;
          pNode = *cNode;
        }
        pSite.clear();
      }
    }
    for (pos = 0; pos < seqLen; pos++) {
      if (allele[pos].size() == 1) {
        allele.erase(pos);
        mutNodes.erase(pos);
      } else {
        coEvol[mutNodes[pos]].push_back(pos);
      }
    }
  }
  map< string, set<string> > fixedMutation() {
    int node;
    vector<int> linked;
    map< string, set<string> > linkages;
    map< int, set< pair<int, int> > >::iterator p;
    for (p = mutNodes.begin(); p != mutNodes.end(); p++) {
      if ((*p).second.size() == 1) {
        pos = (*p).first;
        pNode = (*(*p).second.begin()).first;
        node = (*(*p).second.begin()).second;
        linked.push_back(pNode);
        linked.push_back(node);
        linkages[to_string(pNode) + "~" + to_string(node)].insert(
          as<string>(ancestralSeqs[pNode - 1][pos]) + 
            to_string(pos) +
            as<string>(ancestralSeqs[node - 1][pos])
        );
      }
    }
    int pLinked;
    string assumedLink;
    for (ePath = evolPath.begin(); ePath != evolPath.end(); ePath++) {
      pLinked = root;
      for (cNode = (*ePath).begin() + 1; cNode != (*ePath).end(); cNode++) {
        if (
            find(linked.begin(), linked.end(), *cNode) != linked.end() ||
              find(divPoints.begin(), divPoints.end(), *cNode) != divPoints.end()
        ) {
          assumedLink = to_string(pLinked) + "~" + to_string(*cNode);
          if (linkages.find(assumedLink) == linkages.end()) {
            linkages[assumedLink].clear();
          }
          pLinked = *cNode;
        }
      }
    }
    return linkages;
  }
private:
  int pos, pNode;
  string sitePos, pSite, cSite;
  ListOf<CharacterVector> ancestralSeqs;
  vector< deque<int> > evolPath;
  vector< deque<int> >::iterator ePath;
  deque<int>::iterator cNode;
  set<int> divPoints;
  map< int, set<string> > allele;
  map< int, set< pair<int, int> > > mutNodes;
  map< set< pair<int, int> >, vector<int> > coEvol;
};

// [[Rcpp::export]]
ListOf<IntegerVector> trimTree(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    NumericVector similarity
) {
  Pruner match(tipPaths, alignedSeqs);
  match.setThreshold(as<float>(similarity));
  return wrap(match.groupTips());
}

// [[Rcpp::export]]
ListOf<CharacterVector> mutationPath(
    ListOf<IntegerVector> tipPaths,
    ListOf<CharacterVector> alignedSeqsAR,
    NumericVector similarity
) {
  SiteExplorer match(tipPaths, alignedSeqsAR);
  match.setThreshold(as<float>(similarity));
  match.getEvolPath();
  return wrap(match.fixedMutation());
}
