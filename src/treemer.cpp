#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <Rcpp.h>

#include "treemer.h"

float Treemer::compare(const std::string &query, const std::string &subject) {
    // Get the similarity between two aligned sequences. Site is skipped for
    // total length if both are gaps in alignment
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
    return match/length;
}

Treemer::TipSeqLinker::TipSeqLinker(
    const Rcpp::CharacterVector &sequence,
    const Rcpp::IntegerVector &tipPath
):
    m_seq(Rcpp::as<std::string>(sequence)),
    m_path(tipPath),
    m_tipIndex(tipPath.size() - 1), // last node of a path is the tip node
    m_cIndex(m_tipIndex) {}

void Treemer::TipSeqLinker::reset() {
    m_cIndex = m_tipIndex;
}

void Treemer::TipSeqLinker::proceed() {
    // Proceed towards the root node along the path as the current index should
    // decrement. An index of 0 means reaching the root node. Setting index
    // greater than 1 is to prevent trivial trimming
    if (m_cIndex > 1) { m_cIndex--; }
}

int Treemer::TipSeqLinker::currentClade() const {
    // Current clade is the node at the current index of path
    return m_path[m_cIndex];
}

int Treemer::TipSeqLinker::nextClade() const {
    // Look up the immediate ancestral node aka fake 'proceed'
    if (m_cIndex > 1) {
        return m_path[m_cIndex - 1];
    } else {
        return currentClade(); // The current clade would be root
    }
}

int Treemer::TipSeqLinker::getTip() const {
    // Tip node is the last node in the path
    return m_path[m_tipIndex];
}

int Treemer::TipSeqLinker::getRoot() const {
    // Root node is the first node in the path
    return m_path[0];
}

int Treemer::TipSeqLinker::getSeqLen() const {
    // The length of the aligned sequence. This is going to be done for every
    // sequence to make sure they're aligned
    return m_seq.size();
}

Rcpp::IntegerVector Treemer::TipSeqLinker::getPath() const {
    // The path from current clade towards root node
    return m_path[Rcpp::Range(0, m_cIndex)];
}

std::string Treemer::TipSeqLinker::getSeq() const {
    // The aligned sequence
    return m_seq;
}

char Treemer::TipSeqLinker::siteChar(const int siteIndex) const {
    return m_seq[siteIndex];
}

Treemer::Base::Base(
    const tips &tips,
    const clusters initClusters
):
    m_tips(tips),
    m_clusters(initClusters) {}

Treemer::Base::~Base() {
    // Reset all TipSeqLinkers
    for (tips::const_iterator it = m_tips.begin(); it != m_tips.end(); ++it) {
        (**it).reset();
    }
}

std::map< int, std::vector<int> > Treemer::Base::getTips() const {
    std::map< int, std::vector<int> > res;
    // Get tip from each TipSeqLinker in the list
    for (tips::const_iterator it = m_tips.begin(); it != m_tips.end(); ++it) {
        res[(*it)->currentClade()].push_back((*it)->getTip());
    }
    return res;
}

std::vector<Rcpp::IntegerVector> Treemer::Base::getPaths() const {
    std::vector<Rcpp::IntegerVector> res;
    // Get path from root to trimming-node from each TipSeqLinker in the list
    for (tips::const_iterator it = m_tips.begin(); it != m_tips.end(); ++it) {
        res.push_back((*it)->getPath());
    }
    return res;
}


void Treemer::Base::pruneTree() {
    while (true) {
        clusters oldClusters = m_clusters;
        m_clusters.clear();
        // look down one more node (fake 'proceed')
        // group each tip after new positioning
        for (tips::const_iterator it = m_tips.begin(); it != m_tips.end(); ++it) {
            m_clusters[(*it)->nextClade()].push_back(*it);
        }
        // if no more group 'kissed' each other by a common ancestral node after
        // fake 'proceed', then pruning is done
        if (m_clusters.size() == oldClusters.size()) {
            m_clusters.clear();
            break;
        }
        // only 'kissed' group can do real 'proceed' if a grouping doesn't exist
        // in 'oldClusters' then all tips in that group can 'proceed'
        for (
                clusters::iterator it = m_clusters.begin();
                it != m_clusters.end(); ++it
        ) {
            // assume a group is kissed with another (give it benefit of the
            // doubt)
            bool kissed = true;
            for (
                    clusters::iterator it2 = oldClusters.begin();
                    it2 != oldClusters.end(); ++it2
            ) {
                // a group is 'non-kissed' after fake 'proceed' if it can be
                // found in 'oldClusters'
                if (it->second == it2->second) {
                    kissed = false;
                    // a 'non-kissed' group won't appear twice in 'clusters' so
                    // deleted
                    oldClusters.erase(it2);
                    break;
                }
            }
            // clusters_it group needs to pass some requirement to be qualified
            // 'kissed'
            if (kissed && qualified(it)) {
                for (
                        tips::iterator tips_itr = it->second.begin();
                        tips_itr != it->second.end(); ++tips_itr
                ) { (*tips_itr)->proceed(); }
            }
        }
    }
}

Treemer::BySite::BySite(
    const tips &tips,
    const clusters initClusters,
    const int siteIndex
):
    Base(tips, initClusters),
    m_siteIndex(siteIndex) { pruneTree(); }

std::map<char, Treemer::clusters> Treemer::BySite::siteClusters() const {
    std::map<char, clusters> res;
    for (tips::const_iterator it = m_tips.begin(); it != m_tips.end(); ++it) {
        res[(**it).siteChar(m_siteIndex)][(**it).currentClade()].push_back(*it);
    }
    return res;
}

bool Treemer::BySite::qualified(const clusters::iterator &clusters_it) const {
    tips::const_iterator it  = clusters_it->second.begin();
    char ref_value = (**it).siteChar(m_siteIndex);
    for (++it; it != clusters_it->second.end(); ++it) {
        if ((**it).siteChar(m_siteIndex) != ref_value) {
            return false;
        }
    }
    return true;
}

Treemer::BySimilarity::BySimilarity(
    const tips &tips,
    const clusters initClusters,
    const float simThreshold,
    std::map<std::pair<int, int>, float> &simMatrix
):
    Base(tips, initClusters),
    m_simCut(simThreshold),
    m_compared(&simMatrix) {
    // The similarity threshold should be between 0 and 1
    if (m_simCut <= 0) {
        throw std::invalid_argument("Similarity cannot be lower or equal to 0");
    } else if (m_simCut > 1) {
        throw std::invalid_argument("Similarity cannot be greater than 1");
    }
    // Trim the tree when the object is instantiated
    if (m_simCut != 1) { pruneTree(); }
}

bool Treemer::BySimilarity::qualified(const clusters::iterator &clusters_it) const {
    // Iterate all pairs of tips in the cluster
    for (
            tips::const_iterator it = clusters_it->second.begin();
            it != clusters_it->second.end() - 1; ++it
    ) {
        for (
                tips::const_iterator it2 = it + 1;
                it2 != clusters_it->second.end(); ++it2
        ) {
            std::pair<int, int> pairing = std::make_pair(
                (*it)->getTip(),
                (*it2)->getTip()
            );
            // Retrieve similarity from existing values. Calculate the value
            // otherwise
            float sim = 0;
            std::map<std::pair<int, int>, float>::iterator
                pos = m_compared->find(pairing);
            if (pos != m_compared->end()) {
                sim = pos->second;
            } else {
                sim = compare((*it)->getSeq(), (*it2)->getSeq());
                (*m_compared)[pairing] = sim;
            }
            // Break the loop and disqualified the cluster if applied
            if (sim < m_simCut) {
                return false;
            }
        }
    }
    return true;
}
