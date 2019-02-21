#ifndef SITEPATH_TREEMER_H
#define SITEPATH_TREEMER_H

#include <map>
#include <utility>
#include "util.h"

namespace Treemer {

typedef std::vector<TipSeqLinker> tips;
typedef std::map<int, std::vector<TipSeqLinker *>> clusters;

class Base {
public:
    Base(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
    );
    std::map<int, std::vector<int>> getTips() const;
    std::vector<Rcpp::IntegerVector> getPaths() const;
protected:
    void pruneTree();
    virtual bool qualified(const clusters::iterator &candidate) const = 0;
protected:
    tips m_tips;
    clusters m_clusters;
    clusters m_trueClusters;
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
    bool qualified(const clusters::iterator &candidate) const;
private:
    const int m_site;
};

}

#endif /* SITEPATH_TREEMER_H */