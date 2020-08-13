#include <set>
#include <algorithm>
#include "treemer.h"
#include "minEntropy.h"

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

// [[Rcpp::export]]
Rcpp::ListOf< Rcpp::ListOf<Rcpp::IntegerVector> > summarizeAA(
    const Rcpp::List &allMutations,
    const Rcpp::ListOf<Rcpp::CharacterVector> &allSampledTips,
    const Rcpp::ListOf<Rcpp::CharacterVector> &originalNodeTips,
    const Rcpp::Function &setTxtProgressBar,
    const Rcpp::List &pb
) {
    unsigned int nSampling = allMutations.size();
    const Rcpp::CharacterVector ancestralNodes = originalNodeTips.names();
    std::map< std::string, std::map< std::string, std::map<std::string, int> > > res;
    // Assign amino acid for each ancestral nodes from the sampling result
    for (unsigned int i = 0; i < nSampling;) {
        const Rcpp::CharacterVector sampledTips = allSampledTips[i];
        const Rcpp::List mutations = allMutations[i];
        const Rcpp::CharacterVector sites = mutations.names();
        for (int j = 0; j < mutations.size(); ++j) {
            const Rcpp::List sitePath = mutations[j];
            const std::string site = Rcpp::as<std::string>(sites[j]);
            for (int n = 0; n < sitePath.size(); ++n) {
                const Rcpp::List mutPath = sitePath[n];
                for (int m = 0; m < mutPath.size(); ++m) {
                    // Grab the tips during a fixation and relevant info
                    const Rcpp::IntegerVector tips = mutPath[m];
                    const std::string fixedAA = Rcpp::as<std::string>(tips.attr("AA"));
                    Rcpp::CharacterVector stn = sampledTips[tips - 1];
                    bool qualified = true;
                    // The tips during a fixation should not disappear due to the
                    // sampling in case
                    for (unsigned int it = 0; it < nSampling; ++it) {
                        const Rcpp::CharacterVector overlap = Rcpp::intersect(
                            stn,
                            allSampledTips[it]
                        );
                        if (overlap.size() == 0) {
                            qualified = false;
                            break;
                        }
                    }
                    if (qualified) {
                        // Get the ancestral node on the original tree for each tip
                        // and assign the amino acid.
                        // Duplicated nodes only count once.
                        for (int it = 0; it < ancestralNodes.size(); ++it) {
                            const std::string node = Rcpp::as<std::string>(ancestralNodes[it]);
                            bool nodeIncluded = false;
                            Rcpp::CharacterVector otn = originalNodeTips[node];
                            Rcpp::CharacterVector::iterator stn_itr, otn_itr;
                            // There will be shared tips if node is to be included
                            for (stn_itr = stn.begin(); stn_itr != stn.end();) {
                                otn_itr = std::find(otn.begin(), otn.end(), *stn_itr);
                                if (otn_itr == otn.end()) {
                                    ++stn_itr;
                                } else {
                                    nodeIncluded = true;
                                    // Remove the shared tip for both
                                    otn.erase(otn_itr);
                                    stn_itr = stn.erase(stn_itr);
                                }
                            }
                            // Increase the count for the site-node pair
                            if (nodeIncluded) {
                                ++res[site][node][fixedAA];
                            }
                        }
                    }
                }
            }
        }
        setTxtProgressBar(pb, ++i);
    }
    return Rcpp::wrap(res);
}
