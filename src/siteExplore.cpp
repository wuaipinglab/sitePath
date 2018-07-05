#include "siteExplorer.h"
#include <algorithm>

SiteExplorer::SiteExplorer(
  ListOf<IntegerVector> tipPaths,
  ListOf<CharacterVector> alignedSeqs
): TreeAlignmentMatch(tipPaths, alignedSeqs) {
  ancestralSeqs = alignedSeqs;
  divPoints.insert(root);
}

map< string, set<string> > SiteExplorer::getSitePath(int mode) {
  if (!pruned) {getEvolPath();}
  vector<int> linked;
  map< string, set<string> > linkages;
  switch (mode) {
  case 0: fixedMutationFilter(linked, linkages);
    break;
  case 1: nonFixedFilter(linked, linkages);
    break;
  case 2: coSiteFilter(linked, linkages);
    break;
  }
  int node;
  string assumedLink;
  for (pathIter ePath = evolPath.begin(); ePath != evolPath.end(); ePath++) {
    node = root;
    for (deque<int>::iterator cNode = (*ePath).begin() + 1; cNode != (*ePath).end(); cNode++) {
      if (
          find(linked.begin(), linked.end(), *cNode) != linked.end() ||
            find(divPoints.begin(), divPoints.end(), *cNode) != divPoints.end()
      ) {
        assumedLink = to_string(node) + "~" + to_string(*cNode);
        if (linkages.find(assumedLink) == linkages.end()) {
          linkages[assumedLink].clear();
        }
        node = *cNode;
      }
    }
  }
  return linkages;
}

void SiteExplorer::fixedMutationFilter(
    vector<int> &linked,
    map< string, set<string> > &linkages
) {
  for (modeIter p = mutNodes.begin(); p != mutNodes.end(); p++) {
    if ((*p).second.size() == 1) {
      addToLinkage(
        linked, linkages,
        (*(*p).second.begin()).first,
        (*(*p).second.begin()).second,
        (*p).first
      );
    }
  }
}

void SiteExplorer::coSiteFilter(
    vector<int> &linked,
    map< string, set<string> > &linkages
) {
  for (map< set< pair<int, int> >, vector<int> >::iterator m = coEvol.begin(); m != coEvol.end(); m++) {
    if ((*m).second.size() > 1) {
      for (vector<int>::iterator p = (*m).second.begin(); p != (*m).second.end(); p++) {
        for (set< pair<int, int> >::iterator l = (*m).first.begin(); l != (*m).first.end(); l++) {
          addToLinkage(linked, linkages, (*l).first, (*l).second, *p);
        }
      }
    }
  }
}

void SiteExplorer::nonFixedFilter(
    vector<int> &linked,
    map< string, set<string> > &linkages
) {
  for (modeIter p = mutNodes.begin(); p != mutNodes.end(); p++) {
    if ((*p).second.size() > 1 && allele[(*p).first].size() == 2) {
      for (set< pair<int, int> >::iterator l = (*p).second.begin(); l != (*p).second.end(); l++) {
        addToLinkage(linked, linkages, (*l).first, (*l).second, (*p).first);
      }
    }
  }
}

void SiteExplorer::addToLinkage(
    vector<int> &linked,
    map< string, set<string> > &linkages,
    const int &node1,
    const int &node2,
    const int &siteIndex
) {
  linked.push_back(node1);
  linked.push_back(node2);
  linkages[to_string(node1) + "~" + to_string(node2)].insert(
      as<string>(ancestralSeqs[node1 - 1][siteIndex]) +
        to_string(siteIndex) +
        as<string>(ancestralSeqs[node2 - 1][siteIndex])
  );
}

void SiteExplorer::getEvolPath() {
  for (vector<TipSeqLinker*>::iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
    evolPath.push_back((*tsLinker)->getPath());
  }
  if (simCut != 1) {
    pruneTree();
    map<deque<int>, int> pn;
    for (vector<TipSeqLinker*>::iterator tsLinker = linkers.begin(); tsLinker != linkers.end(); tsLinker++) {
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
  for (pathIter ePath = evolPath.begin(); ePath != evolPath.end(); ePath++) {
    for (pathIter ePath2 = ePath + 1; ePath2 != evolPath.end(); ePath2++) {
      deque<int>::iterator cNode, cNode2;
      for (
          cNode = (*ePath).begin(), cNode2 = (*ePath2).begin();
          (*cNode) == (*cNode2); cNode++, cNode2++
      ) continue;
      if (cNode != (*ePath).begin()) divPoints.insert(*cNode);
    }
    for (int pos = 0; pos < seqLen; pos++) {
      int pNode;
      string pSite, cSite;
      for (deque<int>::iterator cNode = (*ePath).begin(); cNode != (*ePath).end(); cNode++) {
        cSite = ancestralSeqs[(*cNode) - 1][pos];
        allele[pos].insert(cSite);
        if (!pSite.empty() && cSite != pSite) {
          mutNodes[pos].insert(make_pair(pNode, *cNode));
        }
        pSite = cSite;
        pNode = *cNode;
      }
    }
  }
  for (int pos = 0; pos < seqLen; pos++) {
    if (allele[pos].size() == 1) {
      allele.erase(pos);
      mutNodes.erase(pos);
    } else {
      coEvol[mutNodes[pos]].push_back(pos);
    }
  }
}
