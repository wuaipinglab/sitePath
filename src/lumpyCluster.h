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

namespace LumpyCluster {

class Tips {
public:
    Tips(const int tipIndex);
    void fakeAddTip(const int tipIndex);
    float averageSim(const Rcpp::NumericMatrix &simMatrix) const;
    int popOutlier();
private:
    std::vector<int> m_tips;
};

class Clusters {
public:
    Clusters(
        const Rcpp::CharacterVector &alignedSeqs,
        const Rcpp::NumericMatrix &simMatrix,
        const float similarity
    );
    std::vector<int> getClusters();
    void addTip(const int tipIndex);
private:
    const Rcpp::CharacterVector m_alignedSeqs;
    const Rcpp::NumericMatrix m_simMatrix;
    const float m_similarity;
    std::vector<Tips *> m_clusters;
};

}

#endif /* SITEPATH_LUMPCLUSTER_H */
