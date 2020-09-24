/*
 * When finding the terminal of each phylogenetic lineages, the major SNP is
 * thought to be unique to the lineage but sometimes the same SNP could appear
 * in other lineages. In this case, the ancestral node of the tips with the SNP
 * cannot locate to the desired terminal of the lineages. A clustering method
 * similar to CD-HIT http://weizhongli-lab.org/cd-hit/ is used to detect the
 * outlier tip with the SNP or identify the multiple groups having the same SNP.
 */

#ifndef SITEPATH_LUMPCLUSTER_H
#define SITEPATH_LUMPCLUSTER_H

#include <map>
#include <vector>
#include <Rcpp.h>

#include "treemer.h"

namespace LumpyCluster {

typedef std::vector<Treemer::tips> lumpingTips;

class Base {
public:
    Base(
        const Rcpp::NumericMatrix &metricMatrix,
        const int maxSNPnum
    );
    std::vector< std::vector<int> > finalClusters() const;
protected:
    void mergeClusters(
            const Treemer::clusters &clusters,
            const int zValue
    );
    virtual void thresholdOffset(
            const float stdev,
            const int zValue
    ) = 0;
    float clusterCompare(
            const Treemer::tips &query,
            const Treemer::tips &subject
    );
    virtual bool betterMetric(
            const float query,
            const float subject
    ) const = 0;
    virtual bool qualifiedMetric(const float metric) const = 0;
protected:
    const Rcpp::NumericMatrix m_metricMatrix;
    // The merged clusters output
    const int m_maxSNPnum;
    // The threshold for clustering and the metric standard deviation of all
    // tips pairs
    lumpingTips m_merged;
    // The maximum number of SNP
    float m_metricThreshold;
};

/*
 * The clusters tend to merge when their metric is larger.
 */

class BySimMatrix: public Base {
public:
    BySimMatrix(
        const Rcpp::NumericMatrix &simMatrix,
        const Treemer::clusters &clusters,
        const int maxSNPnum,
        const int zValue
    );
protected:
    void thresholdOffset(
            const float stdev,
            const int zValue
    );
    bool betterMetric(
            const float query,
            const float subject
    ) const;
    bool qualifiedMetric(const float metric) const;
};

/*
 * The initiation of clustering by distance is the same as clustering by
 * similarity. Only the test of better metric and qualification is reversed. The
 * clusters tend to merge when their metric is smaller.
 */
class ByDistMatrix: public Base {
public:
    ByDistMatrix(
        const Rcpp::NumericMatrix &distMatrix,
        const Treemer::clusters &clusters,
        const int maxSNPnum,
        const int zValue
    );
protected:
    void thresholdOffset(
            const float stdev,
            const int zValue
    );
    bool betterMetric(
            const float query,
            const float subject
    ) const;
    bool qualifiedMetric(const float metric) const;
};

typedef std::vector< std::vector<int> > tipNodes;

template<class T>
std::map<int, tipNodes> terminalTips(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const Rcpp::NumericMatrix &metricMatrix,
        const Rcpp::IntegerVector &siteIndices,
        const int zValue
);

}

#endif /* SITEPATH_LUMPCLUSTER_H */
