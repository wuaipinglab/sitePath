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
