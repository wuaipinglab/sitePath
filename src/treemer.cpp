#include "treemer.h"

Treemer::Base::Base(
    const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
    const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
):
    m_root(*(tipPaths[0].begin())),
    m_seqLen((Rcpp::as<std::string>(alignedSeqs[0])).size()) {
    // Iterate tipPaths and alignedSeqs to construct a list of TipSeqLinkers
    for (int i = 0; i < tipPaths.size(); i++) {
        TipSeqLinker *tip = new TipSeqLinker(alignedSeqs[i], tipPaths[i]);
        m_tips.push_back(tip);
        // The initial clustering is each tip as a cluster
        m_clusters[tip->getTip()].push_back(tip);

        // The root of each tipPath should be the same
        // The sequences should be of the same length
        if (m_tips[i]->getRoot() != m_root) {
            throw std::invalid_argument("Root in tree paths not equal");
        } else if (m_tips[i]->getSeqLen() != m_seqLen) {
            throw std::invalid_argument("Sequence length not equal");
        }
    }
}

Treemer::Base::~Base() {
    // Release the memory used by TipSeqLinkers
    for (tips::iterator it = m_tips.begin(); it != m_tips.end(); ++it) {
        delete *it;
    }
    m_tips.clear();
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
        for (tips::iterator it = m_tips.begin(); it != m_tips.end(); ++it) {
            m_clusters[(*it)->nextClade()].push_back(*it);
        }
        // if no more group 'kissed' each other by a common ancestral node
        // after fake 'proceed', then pruning is done
        if (m_clusters.size() == oldClusters.size()) {
            m_clusters.clear();
            break;
        }
        // only 'kissed' group can do real 'proceed'
        // if a grouping doesn't exist in 'oldClusters'
        // then all tips in that group can 'proceed'
        for (
                clusters::iterator it = m_clusters.begin();
                it != m_clusters.end(); ++it
        ) {
            // assume a group is kissed with another
            // (give it benefit of the doubt)
            bool kissed = true;
            for (
                    clusters::iterator it2 = oldClusters.begin();
                    it2 != oldClusters.end(); ++it2
            ) {
                // a group is 'non-kissed' after fake 'proceed'
                // if it can be found in 'oldClusters'
                if (it->second == it2->second) {
                    kissed = false;
                    // a 'non-kissed' group won't appear twice in
                    // 'clusters' so deleted
                    oldClusters.erase(it2);
                    break;
                }
            }
            // clusters_it group needs to pass some requirement to
            // be qualified 'kissed'
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
    const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
    const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
    const int site
):
    Base(tipPaths, alignedSeqs),
    m_site(site - 1) { pruneTree(); }

bool Treemer::BySite::qualified(const clusters::iterator &clusters_it) const {
    tips::const_iterator it  = clusters_it->second.begin();
    char ref_value = (*it)->getSeq()[m_site];
    for (++it; it != clusters_it->second.end(); ++it) {
        if ((*it)->getSeq()[m_site] != ref_value) { return false; }
    }
    return true;
}

Treemer::BySimilarity::BySimilarity(
    const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
    const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
    const float simThreshold,
    std::map<std::pair<int, int>, float> &simMatrix
):
    Base(tipPaths, alignedSeqs),
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
            // Retrieve similarity from existing values
            // Calculate the value otherwise
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
