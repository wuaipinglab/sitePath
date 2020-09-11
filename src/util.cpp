/*
 * Here are all the functions exported to the R interface. This is all made
 * possible by the excellent Rcpp package. The cpp11 package has come out and as
 * one of the depending package ggplot2 uses cpp11, migrating from Rcpp to cpp11
 * might be considered in the future.
 */

#include <algorithm>
#include <map>
#include <set>
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
Rcpp::ListOf< Rcpp::ListOf<Rcpp::IntegerVector> > runTreemerBySite(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const Rcpp::IntegerVector &loci
) {
    std::map< int, std::map< int, std::vector<int> > > res;
    for (
            Rcpp::IntegerVector::const_iterator it = loci.begin();
            it != loci.end(); it++
    ) {
        Treemer::BySite match(tipPaths, alignedSeqs, *it);
        res[*it] = match.getTips();
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> majorSNPtips(
        const Rcpp::CharacterVector &alignedSeqs,
        const int minSNPnum
) {
    const int totalTipNum = alignedSeqs.size();

    std::vector< std::vector<int> > res;
    // Select the major SNPs for each site (aa/nt)
    for (int siteIndex = 0; siteIndex < alignedSeqs[0].size(); ++siteIndex) {
        std::map<char, int> snpSummary;
        for (int i = 0; i < totalTipNum; ++i) {
            snpSummary[alignedSeqs[i][siteIndex]]++;
        }
        // Select major SNPs
        for (
                std::map<char, int>::iterator it = snpSummary.begin();
                it != snpSummary.end(); it++
        ) {
            if (it->second > minSNPnum && it->second != totalTipNum) {
                std::vector<int> tips;
                // Find the tips with the SNP and record the site
                for (int i = 0; i < totalTipNum; ++i) {
                    if (alignedSeqs[i][siteIndex] == it->first) {
                        tips.push_back(i+1);
                    }
                }
                res.push_back(tips);
            }
        }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> mergePaths(
        const Rcpp::ListOf<Rcpp::IntegerVector> &paths
) {
    std::vector<Rcpp::IntegerVector> res;
    res.push_back(paths[0]);
    for (int i = 1; i < paths.size(); ++i) {
        bool toAddNew = true;
        Rcpp::IntegerVector::const_iterator q, s;
        for (
                std::vector<Rcpp::IntegerVector>::iterator it = res.begin();
                it != res.end(); ++it
        ) {
            bool toRemoveOld = false;
            q = paths[i].begin(), s = it->begin();
            while (*q == *s) {
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
