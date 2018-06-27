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
    maxPathIndex = tipPath.size() - 1;
    pathIndex = 0;
    tip = tipPath[0];
  }
  void proceed() {
    if (maxPathIndex > pathIndex) {
      pathIndex++;
    }
  };
  float compare(TipSeqLinker *linker) {
    float match, length;
    match = 0, length = 0;
    CharacterVector::iterator query, subject;
    for (query = seq.begin(), subject = linker->seq.begin();
         query != seq.end(); query++, subject++) {
      if (*query == *subject && *query != '-') {match++;}
      length++;
    }
    return match/length;
  };
  vector<string> siteComp(vector<int> sites) {
    vector<string> comp;
    vector<int>::iterator pos;
    for (pos = sites.begin(); pos != sites.end(); pos++) {
      comp.push_back(as<string>(seq[*pos]));
    }
    return comp;
  }
  int nextClade() {
    if (maxPathIndex > pathIndex) {
      return path[pathIndex + 1];
    } else {
      return path[pathIndex];
    }
  }
  int currentClade() {
    return path[pathIndex];
  }
  int getTip() {
    return tip;
  }
  int getRoot() {
    return path[maxPathIndex];
  }
  int getSeqLen() {
    return seq.size();
  }
  deque<int> getPath() {
    deque<int> trimmedPath;
    IntegerVector::iterator it;
    for (it = path.begin() + pathIndex; it != path.end(); it++) {
      trimmedPath.push_front(*it);
    }
    return trimmedPath;
  }
private:
  CharacterVector seq;
  IntegerVector path;
  int pathIndex;
  int maxPathIndex;
  int tip;
};

class TreeAlignmentMatch {
public:
  TreeAlignmentMatch(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs
  ) {
    tipNum = tipPaths.size();
    simCut = 0.9;
    TipSeqLinker *linker;
    int refRoot, refSeqLen;
    for (int i = 0; i < tipNum; i++) {
      linker = new TipSeqLinker(
        alignedSeqs[i],
        tipPaths[i]
      );
      linkerList.push_back(linker);
      clusters[linker->getTip()].push_back(linker);
      root = linker->getRoot();
      seqLen = linker->getSeqLen();
      if (i == 0) {
        refRoot = root;
        refSeqLen = seqLen;
      } else if (refRoot != root) {
        throw invalid_argument("Root in tree path not equal");
      } else if (refSeqLen != seqLen) {
        throw invalid_argument("Sequene length not equal");
      }
    }
  }
  void setThreshold(float similarity = 0.9) {
    if (similarity <= 0) {
      throw invalid_argument("Similarity cannot be lower or equal to 0");
    } else if (similarity > 1) {
      throw invalid_argument("Similarity cannot be greater than 1");
    } else {
      simCut = similarity;
    }
  }
  void setSites(IntegerVector restrictSites = IntegerVector::create()) {
    sites.clear();
    IntegerVector::iterator it;
    for (it = restrictSites.begin(); it != restrictSites.end(); it++) {
      if ((*it) <= 0) {
        throw invalid_argument("Site can't be negative");
      } else if ((*it) > seqLen) {
        throw invalid_argument("Site out of index");
      } else {
        sites.push_back(*it);
      }
    }
  }
  map<string, vector<int> > groupTips () {
    pruneTree();
    map<string, vector<int> > tipCluster;
    map<int, deque<TipSeqLinker*> >::iterator it;
    deque<TipSeqLinker*>::iterator linker;
    for (it = clusters.begin(); it != clusters.end(); it++) {
      for (linker = it->second.begin(); linker != it->second.end(); linker++) {
        tipCluster[to_string(it->first)].push_back((*linker)->getTip());
      }
    }
    return tipCluster;
  }
  map< string, set<string> > getSitePath(
      ListOf<CharacterVector> ancestralSeqs
  ) {
    pruneTree();
    map< string, set<string> > linkages;
    string sitePos, previousSite, currentSite, assumedLink;
    int n, m, siteIndex, pathSize, 
    previousNode, currentNode, previousLinked, currentLinked;
    vector<int> linkedNodes, conserved;
    vector< vector<int> > mutMode;
    map< vector< vector<int> >, vector<string> > mutModeSum;
    map< vector< vector<int> >, vector<string> >::iterator mutModeIt;
    map< string, set<string> > siteMutSum;
    set<int> divPoints;
    deque<int> path;
    CharacterVector nodeSeq;
    evolPath = getEvolPath();
    divPoints = getDivPoints(evolPath);
    for (n = 0; n < evolPath.size(); n++) {
      path = evolPath[n];
      linkedNodes.clear();
      for (siteIndex = 0; siteIndex < seqLen; siteIndex++) {
        pathSize = path.size();
        sitePos = to_string(siteIndex + 1);
        for (m = 0; m < pathSize; m++) {
          currentNode = path[m];
          currentSite = ancestralSeqs[currentNode][siteIndex];
          siteMutSum[sitePos].insert(currentSite);
          if (!previousSite.empty() && currentSite != previousSite) {
            mutMode.push_back(conserved);
            conserved.clear();
            linkedNodes.push_back(previousNode);
            linkedNodes.push_back(currentNode);
            linkages[to_string(previousNode) + "~" + to_string(currentNode)].
            insert(previousSite + sitePos + currentSite);
          }
          conserved.push_back(currentNode);
          previousSite = currentSite;
          previousNode = currentNode;
        }
        mutMode.push_back(conserved);
        mutModeSum[mutMode].push_back(sitePos);
        previousSite.clear();
        conserved.clear();
        mutMode.clear();
      }
      // for (mutModeIt = mutModeSum.begin(); mutModeIt != mutModeSum.end(); mutModeIt++) {
      //   if (mutModeIt->first.size() == 2) {
      //     
      //   }
      // }
      previousLinked = root;
      for (m = 1; m < pathSize; m++) {
        currentLinked = path[m];
        if (
            find(linkedNodes.begin(), linkedNodes.end(), currentNode) != linkedNodes.end() ||
              find(divPoints.begin(), divPoints.end(), currentNode) != divPoints.end()
        ) {
          assumedLink = to_string(previousLinked) + "~" + to_string(currentLinked);
          if (linkages.find(assumedLink) == linkages.end()) {
            linkages[assumedLink].clear();
          }
        }
        previousLinked = currentLinked;
      }
    }
    return linkages;
  }
private:
  int tipNum, root, seqLen;
  float simCut;
  vector<int> sites;
  map<pair<int, int>, float> compared;
  vector<TipSeqLinker*> linkerList;
  map< int, deque<TipSeqLinker*> > clusters;
  vector< deque<int> > evolPath;
  bool qualified(deque<TipSeqLinker*> clstr) {
    float similarity;
    int i, j, length;
    TipSeqLinker *query, *subject;
    vector<string> queryComp, subjectComp;
    pair<int, int> pairing;
    length = clstr.size();
    for (i = 0; i < length; i++) {
      query = clstr[i];
      queryComp = (*query).siteComp(sites);
      for (j = i + 1; j < length; j++) {
        subject = clstr[j];
        subjectComp = (*subject).siteComp(sites);
        pairing = make_pair((*query).getTip(), (*subject).getTip());
        if (compared.find(pairing) != compared.end()) {
          similarity = compared[pairing];
        } else {
          similarity = (*query).compare(subject);
          compared[pairing] = similarity;
        }
        if (similarity < simCut or queryComp != subjectComp) {
          return false;
        }
      }
    }
    return true;
  }
  void pruneTree() {
    bool exist;
    map< int, deque<TipSeqLinker*> > oldCluster;
    vector<TipSeqLinker*>::iterator it;
    map< int, deque<TipSeqLinker*> >::iterator it1, it2;
    deque<TipSeqLinker*>::iterator linker;
    while (true) {
      oldCluster = clusters;
      clusters.clear();
      for (it = linkerList.begin(); it != linkerList.end(); it++) {
        clusters[(*it)->nextClade()].push_back(*it);
      }
      if (clusters.size() == oldCluster.size()) {
        clusters.clear();
        for (it = linkerList.begin(); it != linkerList.end(); it++) {
          clusters[(*it)->currentClade()].push_back(*it);
        }
        break;
      }
      for (it1 = clusters.begin(); it1 != clusters.end(); it1++) {
        exist = false;
        for (it2 = oldCluster.begin(); it2 != oldCluster.end(); it2++) {
          if (it1->second == it2->second) {
            exist = true;
            oldCluster.erase(it2);
            break;
          }
        }
        if (!exist && qualified(it1->second)) {
          for (linker = it1->second.begin(); linker != it1->second.end(); linker++) {
            (*linker)->proceed();
          }
        }
      }
    }
  }
  vector<deque<int> > getEvolPath() {
    vector<deque<int> > tipPath;
    vector<TipSeqLinker*>::iterator it;
    map<deque<int>, int> pathCount;
    for (it = linkerList.begin(); it != linkerList.end(); it++) {
      pathCount[(*it)->getPath()]++;
    }
    map<deque<int>, int>::iterator path;
    for (path = pathCount.begin(); path != pathCount.end(); path++) {
      if (path->first.size() != 1 && path->second > 1) {
        tipPath.push_back(path->first);
      }
    }
    return tipPath;
  }
  set<int> getDivPoints(vector<deque<int> > evolPath) {
    set<int> divPoints;
    divPoints.insert(root);
    deque<int> path1, path2;
    deque<int>::iterator it1, it2;
    int pathNum, divPoint;
    pathNum = evolPath.size();
    for (int i = 0; i < pathNum; i++) {
      path1 = evolPath[i];
      for (int j = i + 1; j < pathNum; j++) {
        path2 = evolPath[j];
        for (it1 = path1.begin(), it2 = path2.begin(); (*it1) == (*it2); it1++, it2++) {
          divPoint = *it1;
        }
        divPoints.insert(divPoint);
      }
    }
    return divPoints;
  }
};

// class MutationFilter {
// public:
//   virtual ~MutationFilter() {};
//   virtual map<string, set<string> > getFiltered() = 0;
// protected:
//   int site;
//   map<string, set<string> > linkages;
//   void cleanUpLinkages() {
//     for (it = linkages.begin(); it != linkages.end(); it++) {
//       if (it->second.empty()) {
//         unlinked.push_back(it->first);
//       }
//     }
//     unlinkedNum = unlinked.size();
//     for (i = 0; i < unlinkedNum; i++) {
//       for (j = i + 1; j < unlinkedNum; j++) {
//         continue;
//       }
//     }
//   }
// private:
//   int i, j, unlinkedNum;
//   string previousUnlinked, currentUnlinked;
//   vector<string> unlinked;
//   map<string, set<string> >::iterator it;
// };
// 
// class Fixation: MutationFilter {
// public:
//   Fixation(map<string, set<string> > mutations) {
//     linkages = mutations;
//   }
//   map<string, set<string> > getFiltered() {
//     for (int i = 0; i < mutPath.size(); i++) {
//       muts = mutPath[i];
//       for (s = muts.begin(); s != muts.end(); s++) {
//         site = stoi(as<string>(*s).substr(1, (*s).size()));
//         mutSummary[site]++;
//       }
//     }
//     for (l = links.begin(); l != links.end(); l++) {
//       muts = mutPath[as<string>(*l)];
//       set<string> qualified;
//       for (s = muts.begin(); s != muts.end(); s++) {
//         site = stoi(as<string>(*s).substr(1, (*s).size()));
//         if (mutSummary[site] == 1) {
//           qualified.insert(as<string>(*s));
//         }
//       }
//       linkages[as<string>(*l)] = qualified;
//     }
//     return linkages;
//   };
// private:
//   map<int, int> mutSummary;
// };

// [[Rcpp::export]]
ListOf<IntegerVector> trimTree(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    NumericVector similarity
) {
  TreeAlignmentMatch match(tipPaths, alignedSeqs);
  match.setThreshold(as<float>(similarity));
  return wrap(match.groupTips());
}

// [[Rcpp::export]]
ListOf<CharacterVector> mutationPath(
    ListOf<IntegerVector> tipPaths,
    ListOf<CharacterVector> alignedSeqsAR,
    NumericVector similarity
) {
  TreeAlignmentMatch match(tipPaths, alignedSeqsAR);
  match.setThreshold(as<float>(similarity));
  return wrap(match.getSitePath(alignedSeqsAR));
}
