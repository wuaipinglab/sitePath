#ifndef SITEPATH_UTIL_H
#define SITEPATH_UTIL_H

#include <string>
#include <iostream>
#include <Rcpp.h>

const float compare(const std::string &query, const std::string &subject);

class TipSeqLinker {
public:
    TipSeqLinker(
        const Rcpp::CharacterVector &sequence,
        const Rcpp::IntegerVector &tipPath
    );
    void proceed();
    const int nextClade() const;
    const int currentClade() const;
    const int getTip() const;
    const int getRoot() const;
    const int getSeqLen() const;
    Rcpp::IntegerVector getPath() const;
    std::string getSeq() const;
private:
    const std::string m_seq;
    const Rcpp::IntegerVector m_path;
    const int m_tipIndex;
    int m_cIndex;
};

#endif // SITEPATH_UTIL_H