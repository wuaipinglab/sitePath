/*
 * Here are all the functions exported to the R interface. This is all made
 * possible by the excellent Rcpp package. The cpp11 package has come out and as
 * one of the depending package ggplot2 uses cpp11, migrating from Rcpp to cpp11
 * might be considered in the future.
 */

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <Rcpp.h>

#include "lumpyCluster.h"
#include "minEntropy.h"
#include "treemer.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix getSimilarityMatrix(
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
) {
    int dim = alignedSeqs.size();
    Rcpp::NumericMatrix simMatrix(dim, dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = i; j < dim; ++j) {
            if (i == j) { simMatrix(i, j) = 1; } else {
                simMatrix(j, i) = simMatrix(i, j) = Treemer::compare(
                    Rcpp::as<std::string>(alignedSeqs[i]),
                    Rcpp::as<std::string>(alignedSeqs[j])
                );
            }
        }
    }
    return simMatrix;
}

// [[Rcpp::export]]
Rcpp::ListOf< Rcpp::ListOf<Rcpp::IntegerVector> > majorSNPtips(
        const Rcpp::CharacterVector &alignedSeqs,
        const Rcpp::IntegerVector &siteIndices,
        const int minSNPnum
) {
    const int totalTipNum = alignedSeqs.size();
    std::map<int, std::map<std::string, std::vector<int> > > res;
    // Select the major SNPs for each site (aa/nt)
    for (
            Rcpp::IntegerVector::const_iterator it = siteIndices.begin();
            it != siteIndices.end(); it++
    ) {
        int siteIndex = *it;
        std::map<char, int> snpSummary;
        for (int i = 0; i < totalTipNum; ++i) {
            snpSummary[alignedSeqs[i][siteIndex]]++;
        }
        // Select major SNPs
        for (
                std::map<char, int>::iterator it = snpSummary.begin();
                it != snpSummary.end(); it++
        ) {
            if (it->second > minSNPnum && it->second < totalTipNum/2) {
                // Find the tips with the SNP and record the site
                for (int i = 0; i < totalTipNum; ++i) {
                    if (alignedSeqs[i][siteIndex] == it->first) {
                        res[siteIndex+1][std::string(1, it->first)].push_back(i+1);
                    }
                }
            }
        }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::ListOf< Rcpp::ListOf<Rcpp::IntegerVector> > terminalTipsBySim(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const Rcpp::NumericMatrix &metricMatrix,
        const Rcpp::IntegerVector &siteIndices,
        const int zValue
) {
    using namespace LumpyCluster;
    std::map<int, tipNodes> res = terminalTips<BySimMatrix>(
        tipPaths,
        alignedSeqs,
        metricMatrix,
        siteIndices,
        zValue
    );
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector> > terminalTipsByDist(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const Rcpp::NumericMatrix &metricMatrix,
        const Rcpp::IntegerVector &siteIndices,
        const int zValue
) {
    using namespace LumpyCluster;
    std::map<int, tipNodes> res = terminalTips<ByDistMatrix>(
        tipPaths,
        alignedSeqs,
        metricMatrix,
        siteIndices,
        zValue
    );
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> mergePaths(
        const Rcpp::ListOf<Rcpp::IntegerVector> &paths
) {
    // The list of path node to be output
    std::vector<Rcpp::IntegerVector> res;
    // The first path node is assumed to be added
    res.push_back(paths[0]);
    // Iterate the list to remove path node that is subset to other path node
    for (int i = 1; i < paths.size(); ++i) {
        // Assume the coming path node is not a subset to any
        bool toAddNew = true;
        // Iterate each existing path node and compare with the coming path node
        Rcpp::IntegerVector::const_iterator q, s;
        for (
                std::vector<Rcpp::IntegerVector>::iterator it = res.begin();
                it != res.end(); ++it
        ) {
            // Assume the current existing path node is subset to the coming
            // path node
            bool toRemoveOld = false;
            q = paths[i].begin(), s = it->begin();
            while (*q == *s) {
                // The path node is considered subset when the part beside the
                // terminal node has overlap with other
                ++q, ++s;
                if (s == it->end()) {
                    toRemoveOld = true;
                    break;
                }
                if (q == paths[i].end()) {
                    toAddNew = false;
                    break;
                }
            }
            if (toRemoveOld) {
                res.erase(it);
                break;
            }
            if (!toAddNew) {
                break;
            }
        }
        if (toAddNew) {
            res.push_back(paths[i]);
        }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::IntegerVector divergentNode(
        const Rcpp::ListOf<Rcpp::IntegerVector> &paths
) {
    std::set<int> res;
    for (int i = 0; i < paths.size() - 1; i++) {
        for (int j = i + 1; j < paths.size(); j++) {
            Rcpp::IntegerVector::const_iterator
            q  = paths[i].begin(), s = paths[j].begin();
            do { q++, s++; } while (*q == *s);
            if (--q != paths[i].begin()) { res.insert(*q); }
        }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::IntegerVector getReference(
        const std::string &refSeq,
        const char gapChar
) {
    std::vector<int> res;
    for (unsigned int i = 0; i < refSeq.size(); i++) {
        if (refSeq[i] != gapChar) { res.push_back(i + 1); }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::IntegerVector tableAA(
        const Rcpp::CharacterVector &seqs,
        const int siteIndex
) {
    // Summarize the AAs for tree tips of a node
    std::map<std::string, int> res;
    for (int i = 0; i < seqs.size(); ++i) {
        res[std::string(1, seqs[i][siteIndex])]++;
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> minEntropyByInserting(
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries,
        const unsigned int minEffectiveSize,
        const unsigned int searchDepth
) {
    using namespace MinEntropy;
    SearchTree<Segmentor> iSearch(
            minEffectiveSize,
            searchDepth,
            nodeSummaries
    );
    iSearch.search();
    return updatedSegmentation(nodeSummaries, iSearch.getFinal());
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> minEntropyByDeleting(
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries,
        const unsigned int minEffectiveSize,
        const unsigned int searchDepth
) {
    using namespace MinEntropy;
    SearchTree<Amalgamator> dSearch(
            minEffectiveSize,
            searchDepth,
            nodeSummaries
    );
    dSearch.search();
    return updatedSegmentation(nodeSummaries, dSearch.getFinal());
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> minEntropyByComparing(
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries,
        const unsigned int minEffectiveSize,
        const unsigned int searchDepth
) {
    using namespace MinEntropy;
    SearchTree<Segmentor> iSearch(
            minEffectiveSize,
            searchDepth,
            nodeSummaries
    );
    iSearch.search();
    SearchTree<Amalgamator> dSearch(
            minEffectiveSize,
            searchDepth,
            nodeSummaries
    );
    dSearch.search();
    segment iFinal = iSearch.getFinal(), dFinal = dSearch.getFinal();
    segment final;
    if (iFinal.size() > dFinal.size()) {
        final = iFinal;
    } else if (iFinal.size() == dFinal.size()) {
        final = iSearch.getFinal();
    } else {
        final = dFinal;
    }
    return updatedSegmentation(nodeSummaries, final);
}
