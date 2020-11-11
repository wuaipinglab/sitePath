#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "lumpyCluster.h"
#include "treemer.h"

LumpyCluster::Base::Base(
    const Rcpp::NumericMatrix &metricMatrix,
    const int maxSNPnum
):
    m_metricMatrix(metricMatrix),
    m_maxSNPnum(maxSNPnum) {}

LumpyCluster::tipGrouping LumpyCluster::Base::finalClusters() const {
    typedef std::vector<Rcpp::IntegerVector> pathList;
    // The clustering is already known. Need to extract the tip node so can be
    // used in R
    tipGrouping res;
    res.reserve(m_merged.size());
    // Iterate each cluster
    for (
            lumpingTips::const_iterator it = m_merged.begin();
            it != m_merged.end(); ++it
    ) {
        std::vector<int> tipNodes;
        tipNodes.reserve(it->size());
        // The node path for each tip needs to be collected to infer the lineage
        // path of the tip cluster
        pathList nodepaths;
        // Iterate each tip to get its node
        for (
                Treemer::tips::const_iterator tips_itr = it->begin();
                tips_itr != it->end(); ++tips_itr
        ) {
            tipNodes.push_back((**tips_itr).getTip());
            nodepaths.push_back((**tips_itr).getPath());
        }
        // Turn into R object so node path can be set as attribute
        Rcpp::IntegerVector toAdd = Rcpp::wrap(tipNodes);
        if (tipNodes.size() == 1) {
            toAdd.attr("nodepath") = nodepaths[0];
        } else {
            int an = 0;
            while (true) {
                // The node path of the first tip as reference
                pathList::iterator np_itr = nodepaths.begin();
                int refNode = np_itr->at(an);
                // Iterate node path of every tip to find the divergent point
                for (++np_itr; np_itr != nodepaths.end(); ++np_itr) {
                    if (np_itr->at(an) != refNode) {
                        toAdd.attr("nodepath") = (*np_itr)[Rcpp::Range(0, an-1)];
                        goto NODEPATH_FOUND;
                    }
                }
                an++;
            }
        }
        NODEPATH_FOUND: res.push_back(toAdd);
    }
    return res;
}

void LumpyCluster::Base::mergeClusters(
        const Treemer::clusters &clusters,
        const int zValue
) {
    // The fist raw cluster
    Treemer::clusters::const_iterator it = clusters.begin();
    m_merged.push_back(it->second);
    // There is no need to merge when there is only one cluster
    if (clusters.size() == 1) {
        return;
    }
    // The reformatted raw clusters input
    lumpingTips rawClusters;
    // Initiate the merged clusters with the first raw clusters so the size of
    // the reformatted raw cluster will be one size less
    rawClusters.reserve(clusters.size()-1);
    // Collect all tips for later calculating the average paired metric
    Treemer::tips allTips = it->second;
    // Calculate the threshold value and reformat the data structure from map to
    // vector
    for (++it; it != clusters.end(); ++it) {
        // Reformat
        rawClusters.push_back(it->second);
        // Collect tips
        allTips.insert(
            allTips.begin(),
            it->second.begin(),
            it->second.end()
        );
    }
    if (allTips.size() >= m_maxSNPnum) {
        return;
    }
    // Collect sequence similarity between all pairs of tips
    std::vector<float> allMetric;
    for (
            Treemer::tips::iterator it1 = allTips.begin();
            it1 != allTips.end()-1; ++it1
    ) {
        for (
                Treemer::tips::iterator it2 = it1 + 1;
                it2 != allTips.end(); ++it2
        ) {
            float pairMetric = m_metricMatrix(
                (**it1).getTip()-1,
                (**it2).getTip()-1
            );
            allMetric.push_back(pairMetric);
        }
    }
    // Find the median as the threshold
    std::nth_element(
        allMetric.begin(),
        allMetric.begin()+allMetric.size()/2,
        allMetric.end()
    );
    m_metricThreshold = allMetric[allMetric.size()/2];
    // Iterate the raw clusters from treemerBySite as the candidate for merging
    // or adding
    for (
            lumpingTips::const_iterator candidate = rawClusters.begin();
            candidate != rawClusters.end(); ++candidate
    ) {
        // Iterate the merged clusters to find the best one for the current raw
        // candidate cluster. But the merging metric has to pass the threshold.
        // Assume the first one is the best to be merged with the candidate.
        lumpingTips::iterator toMerge = m_merged.begin();
        float bestMetric = clusterCompare(*candidate, *toMerge);
        for (
                lumpingTips::iterator it = m_merged.begin()+1;
                it != m_merged.end(); ++it
        ) {
            float metric = clusterCompare(*candidate, *it);
            // Replace the best if the current is better to be merged with the
            // candidate
            if (betterMetric(metric, bestMetric)) {
                toMerge = it;
                bestMetric = metric;
            }
        }
        // The best metric has to pass the threshold or a new cluster is formed
        if (qualifiedMetric(bestMetric)) {
            // The raw candidate cluster is merged
            toMerge->insert(
                    toMerge->end(),
                    candidate->begin(),
                    candidate->end()
            );
        } else {
            // A cluster is formed
            m_merged.push_back(*candidate);
        }
    }
}

float LumpyCluster::Base::clusterCompare(
        const Treemer::tips &query,
        const Treemer::tips &subject
) {
    int count = 0;
    float metricSum = 0.0;
    // The average metric between tip pairs from each two clusters to
    // approximate the metric between the two clusters
    for (
            Treemer::tips::const_iterator it1 = query.begin();
            it1 != query.end(); ++it1
    ) {
        for (
                Treemer::tips::const_iterator it2 = subject.begin();
                it2 != subject.end(); ++it2
        ) {
            metricSum += m_metricMatrix(
                (**it1).getTip()-1,
                (**it2).getTip()-1
            );
            count++;
        }
    }
    return metricSum/count;
}

LumpyCluster::BySimMatrix::BySimMatrix(
    const Rcpp::NumericMatrix &simMatrix,
    const Treemer::clusters &clusters,
    const int maxSNPnum,
    const int zValue
):
    Base(simMatrix, maxSNPnum) {
    mergeClusters(clusters, zValue);
}

void LumpyCluster::BySimMatrix::thresholdOffset(
        const float stdev,
        const int zValue
) {
    // When using similarity, the number goes up to indicate closer
    // relationship.
    m_metricThreshold += stdev*zValue;
}

bool LumpyCluster::BySimMatrix::betterMetric(
        const float query,
        const float subject
) const {
    return query > subject;
}

bool LumpyCluster::BySimMatrix::qualifiedMetric(const float metric) const {
    return metric > m_metricThreshold;
}

LumpyCluster::ByDistMatrix::ByDistMatrix(
    const Rcpp::NumericMatrix &distMatrix,
    const Treemer::clusters &clusters,
    const int maxSNPnum,
    const int zValue
):
    Base(distMatrix, maxSNPnum) {
    mergeClusters(clusters, zValue);
}

void LumpyCluster::ByDistMatrix::thresholdOffset(
        const float stdev,
        const int zValue
) {
    // When using distance, the number goes down to indicate closer
    // relationship.
    m_metricThreshold -= stdev*zValue;
}

bool LumpyCluster::ByDistMatrix::betterMetric(
        const float query,
        const float subject
) const {
    return query < subject;
}

bool LumpyCluster::ByDistMatrix::qualifiedMetric(const float metric) const {
    return metric < m_metricThreshold;
}

template<class T>
std::map<int, LumpyCluster::tipGrouping> LumpyCluster::terminalTips(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const Rcpp::NumericMatrix &metricMatrix,
        const Rcpp::IntegerVector &siteIndices,
        const int zValue
) {
    const int totalNum = tipPaths.size();
    Treemer::tips tips;
    Treemer::clusters initClusters;
    // Iterate tipPaths and alignedSeqs to construct a list of TipSeqLinkers
    for (int i = 0; i < totalNum; i++) {
        Treemer::TipSeqLinker *tip = new Treemer::TipSeqLinker(
            alignedSeqs[i], tipPaths[i]
        );
        tips.push_back(tip);
        // The initial clustering is each tip as a cluster
        initClusters[tip->getTip()].push_back(tip);
        // The root of each tipPath should be the same. The sequences should be
        // of the same length
        if (tips[i]->getRoot() != tips[0]->getRoot()) {
            throw std::invalid_argument("Root in tree paths not equal");
        } else if (tips[i]->getSeqLen() != tips[0]->getSeqLen()) {
            throw std::invalid_argument("Sequence length not equal");
        }
    }
    std::map<int, tipGrouping> res;
    const int maxSNPnum = totalNum/2;
    for (
            Rcpp::IntegerVector::const_iterator it = siteIndices.begin();
            it != siteIndices.end(); it++
    ) {
        // Use treemer to group tips having monophyletic relationship
        Treemer::BySite matched(tips, initClusters, *it);
        // Group the treemer result by aa/nt of the site
        Treemer::siteClusters clsByAA = matched.getSiteClusters();
        for (
                Treemer::siteClusters::iterator clusters_itr = clsByAA.begin();
                clusters_itr != clsByAA.end(); ++clusters_itr
        ) {
            T merger(
                    metricMatrix,
                    clusters_itr->second,
                    maxSNPnum,
                    0
            );
            tipGrouping mergedCls = merger.finalClusters();
            for (
                    tipGrouping::iterator cls_itr = mergedCls.begin();
                    cls_itr != mergedCls.end(); ++cls_itr
            ) {
                if (cls_itr->size() > 1) {
                    res[*it].push_back(*cls_itr);
                }
            }
        }
    }
    // Manually free the memory allocated by new
    for (Treemer::tips::iterator it = tips.begin(); it != tips.end(); ++it) {
        delete *it;
    }
    return res;
}

template std::map<int, LumpyCluster::tipGrouping>
    LumpyCluster::terminalTips<LumpyCluster::BySimMatrix>(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const Rcpp::NumericMatrix &metricMatrix,
        const Rcpp::IntegerVector &siteIndices,
        const int zValue
    );

template std::map<int, LumpyCluster::tipGrouping>
    LumpyCluster::terminalTips<LumpyCluster::ByDistMatrix>(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const Rcpp::NumericMatrix &metricMatrix,
        const Rcpp::IntegerVector &siteIndices,
        const int zValue
    );
