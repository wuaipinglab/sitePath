#include <cmath>
#include <algorithm>
#include "util.h"

const float compare(const std::string &query, const std::string &subject) {
    // Get the similarity between two aligned sequences.
    // Site is skipped for total length if both are gaps in alignment
    float match = 0.0, length = 0.0;
    for (
            std::string::const_iterator q = query.begin(), s = subject.begin();
            q != query.end(); ++q, ++s
    ) {
        if (*q != '-' && *s != '-') {
            length++;
            if (*q == *s) { match++; }
        }
    }
    return match / length;
}

TipSeqLinker::TipSeqLinker(
    const Rcpp::CharacterVector &sequence,
    const Rcpp::IntegerVector &tipPath
):
    m_seq(Rcpp::as<std::string>(sequence)),
    m_path(tipPath),
    m_tipIndex(tipPath.size() - 1), // last node of a path is the tip node
    m_cIndex(m_tipIndex) {}

void TipSeqLinker::proceed() {
    // Proceed towards the root node along the path
    // as the current index should decrement.
    // An index of 0 means reaching the root node.
    // Setting index greater than 1 is to prevent trivial trimming
    if (m_cIndex > 1) { m_cIndex--; }
}

const int TipSeqLinker::currentClade() const {
    // Current clade is the node at the current index of path
    return m_path[m_cIndex];
}

const int TipSeqLinker::nextClade() const {
    // Look up the immediate ancestral node
    // aka fake 'proceed'
    if (m_cIndex > 1) {
        return m_path[m_cIndex - 1];
    } else {
        return currentClade(); // The current clade would be root
    }
}

const int TipSeqLinker::getTip() const {
    // Tip node is the last node in the path
    return m_path[m_tipIndex];
}

const int TipSeqLinker::getRoot() const {
    // Root node is the first node in the path
    return m_path[0];
}

const int TipSeqLinker::getSeqLen() const {
    // The length of the aligned sequence
    // This is going to be done for every sequence to make sure they're aligned
    return m_seq.size();
}

Rcpp::IntegerVector TipSeqLinker::getPath() const {
    // The path from current clade towards root node
    return m_path[Rcpp::Range(0, m_cIndex)];
}

std::string TipSeqLinker::getSeq() const {
    // The aligned sequence
    return m_seq;
}

float shannonEntropy(const aaSummary &values) {
    int total = 0;
    for (
            aaSummary::const_iterator it = values.begin();
            it != values.end(); ++it
    ) {
        total += it->second;
    }
    float res = 0;
    for (
            aaSummary::const_iterator it = values.begin();
            it != values.end(); ++it
    ) {
        float p = it->second / static_cast<float>(total);
        res -= p * std::log(p);
    }
    return res;
}

Segmentor::Segmentor(
    const segment all,
    const segment terminal,
    const std::vector<aaSummary> &aaSummaries
):
    m_used(terminal),
    m_open(all)
{
    m_entropy = this->totalEntropy(aaSummaries);
}

Segmentor::Segmentor(
    const Segmentor *parent,
    const unsigned int i,
    const std::vector<aaSummary> &aaSummaries
) {
    m_used = this->getUsed(parent, i);
    m_open = this->getOpen(parent, i);
    m_entropy = this->totalEntropy(aaSummaries);
}

const segment Segmentor::getUsed(
        const Segmentor *parent,
        const unsigned int i
) const {
    segment res = parent->m_used;
    res.push_back(parent->m_open.at(i));
    std::sort(res.begin(), res.end());
    return res;
}

const segment Segmentor::getOpen(
        const Segmentor *parent,
        const unsigned int i
) const {
    segment res = parent->m_open;
    res.erase(res.begin() + i);
    return res;
}

const float Segmentor::totalEntropy(
        const std::vector<aaSummary> &aaSummaries
) const {
    float res = 0;
    segIndex start = 0;
    for (
            segment::const_iterator m_used_itr = m_used.begin();
            m_used_itr != m_used.end(); ++m_used_itr
    ) {
        aaSummary values;
        for (unsigned int i = start; i < *m_used_itr; ++i) {
            const aaSummary toBeCombined = aaSummaries.at(i);
            for (
                    aaSummary::const_iterator it = toBeCombined.begin();
                    it != toBeCombined.end(); ++it
            ) {
                values[it->first] += it->second;
            }
        }
        res += shannonEntropy(values);
        start = *m_used_itr;
    }
    return res;
}
