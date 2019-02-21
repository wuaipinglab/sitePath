#ifndef SITEPATH_PRUNER_H
#define SITEPATH_PRUNER_H

#include <map>
#include <utility>
#include "util.h"

class TreeAlignmentMatch {
public:
    TreeAlignmentMatch(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
    );
    virtual ~TreeAlignmentMatch();
    std::map< int, std::vector<int> > getTips() const;
    std::vector<Rcpp::IntegerVector> getPaths() const;
protected:
    std::vector<TipSeqLinker*> m_linkers;
    std::map< int, std::vector<TipSeqLinker*> > m_clusters;
    std::map< int, std::vector<TipSeqLinker*> > m_trueCluster;
    void pruneTree();
    virtual const bool qualified(
            const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate
    ) = 0;
private:
    const int m_root, m_seqLen;
};

class Pruner: public TreeAlignmentMatch {
public:
    Pruner(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const float simThreshold,
        std::map<std::pair<int, int>, float> &simMatrix
    );
private:
    const float m_simCut;
    std::map<std::pair<int, int>, float> m_compared;
    const bool qualified(
            const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate
    );
};

class CustomizablePruner: public TreeAlignmentMatch {
public:
    CustomizablePruner(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const std::map< int, std::vector<int> > &treeEdge,
        std::map<std::pair<int, int>, float> &simMatrix,
        const Rcpp::Function &customQualifyFunc
    );
private:
    std::map< int, std::vector<int> > m_nodeLink;
    std::map<std::pair<int, int>, float> m_compared;
    Rcpp::Function m_qualifyFunc;
    const bool qualified(
            const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate
    );
};

#endif // SITEPATH_PRUNER_H