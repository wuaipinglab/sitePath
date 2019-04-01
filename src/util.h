#ifndef SITEPATH_UTIL_H
#define SITEPATH_UTIL_H

#include <map>
#include <string>
#include <vector>
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

typedef unsigned int segIndex;
typedef std::vector<segIndex> segment;
typedef std::map<std::string, int> aaSummary;

float shannonEntropy(const aaSummary &values);

// TODO: Need a way to actually get the segmented list
// so "minEffectiveSize" can come in to play the role.
class Segmentor {
public:
    static std::vector<aaSummary> aaSummaries;
public:
    segment m_used;
    segment m_open;
    float m_entropy;
public:
    Segmentor(
        const segment all,
        const segIndex terminal
    );
    Segmentor(
        const Segmentor *parent,
        const unsigned int i
    );
private:
    const segment getUsed(
            const Segmentor *parent,
            const unsigned int i
    ) const;
    const segment getOpen(
            const Segmentor *parent,
            const unsigned int i
    ) const;
    const float totalEntropy() const;
};

#endif // SITEPATH_UTIL_H
