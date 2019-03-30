#ifndef SITEPATH_TREEMER_H
#define SITEPATH_TREEMER_H

#include <utility>
#include "util.h"

namespace Treemer {

typedef std::vector<TipSeqLinker *> tips;
typedef std::map< int, std::vector<TipSeqLinker *> > clusters;

class Base {
public:
    Base(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
    );
    virtual ~Base();
    std::map< int, std::vector<int> > getTips() const;
    std::vector<Rcpp::IntegerVector> getPaths() const;
protected:
    void pruneTree();
    virtual bool qualified(const clusters::iterator &clusters_it) const = 0;
protected:
    tips m_tips;
    clusters m_clusters;
private:
    const int m_root, m_seqLen;
};

class BySite: public Base {
public:
    BySite(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const int site
    );
private:
    bool qualified(const clusters::iterator &clusters_it) const;
private:
    const int m_site;
};

class BySimilarity: public Base {
public:
    BySimilarity(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const float simThreshold,
        std::map<std::pair<int, int>, float> &simMatrix
    );
protected:
    const float m_simCut;
    std::map<std::pair<int, int>, float> *m_compared;
private:
    bool qualified(const clusters::iterator &clusters_it) const;
};

class ByAverageSimilarity: public BySimilarity {
public:
    ByAverageSimilarity(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        const float simThreshold,
        std::map<std::pair<int, int>, float> &simMatrix
    );
private:
    bool qualified(const clusters::iterator &clusters_it) const;
};

}

#endif /* SITEPATH_TREEMER_H */