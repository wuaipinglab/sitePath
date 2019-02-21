#include "treemer.h"

Treemer::Base::Base(
    const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
    const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
): 
    m_root(*(tipPaths[0].begin())),
    m_seqLen((Rcpp::as<std::string>(alignedSeqs[0])).size()) 
{
    for (int i = 0; i < tipPaths.size(); i++) {
        m_tips.push_back(TipSeqLinker(alignedSeqs[i], tipPaths[i]));
        m_clusters[m_tips[i].getTip()].push_back(&m_tips[i]);
        
        if (m_tips[i].getRoot() != m_root) {
            throw std::invalid_argument("Root in trees path not equal");
        } else if (m_tips[i].getSeqLen() != m_seqLen) {
            throw std::invalid_argument("Sequence length not equal");
        }
    }
}

std::map<int, std::vector<int>> Treemer::Base::getTips() const {
    std::map<int, std::vector<int>> res;
    for (const auto &tip: m_tips) { 
        res[tip.currentClade()].push_back(tip.getTip());
    }
    return res;
}

std::vector<Rcpp::IntegerVector> Treemer::Base::getPaths() const {
    std::vector<Rcpp::IntegerVector> res;
    for (const auto &tip: m_tips) { res.push_back(tip.getPath()); }
    return res;
}

void Treemer::Base::pruneTree() {
    clusters oldClusters;
    while (true) {
        for (auto &tip: m_tips) {
            m_trueClusters[tip.currentClade()].push_back(&tip);
        }
        oldClusters = m_clusters;
        m_clusters.clear();
        // look down one more node (fake 'proceed') 
        // for each tip after new positioning
        for (auto tip: m_tips) {
            //  group tips by fake 'proceed' node
            m_clusters[tip.nextClade()].push_back(&tip);
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
        for (auto it = m_clusters.begin(); it != m_clusters.end(); ++it) {
            // assume a group is kissed with another (give it benefit of the doubt)
            bool kissed = true;
            for (auto it2 = oldClusters.begin(); it2 != oldClusters.end(); ++it2) {
                // a group is 'non-kissed' after fake 'proceed' 
                // if it can be found in 'oldClusters'
                if (it->second == it2->second) {
                    kissed = false;
                    // a 'non-kissed' group won't appear twice in 'clusters' so deleted
                    oldClusters.erase(it2);
                    break;
                }
            }
            // candidate group needs to pass some requirement to be qualified 'kissed'
            if (kissed and qualified(it)) {
                for (auto tip: it->second) { tip->proceed(); }
            }
        }
        m_trueClusters.clear();
        oldClusters.clear();
    }
}

Treemer::BySite::BySite(
    const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
    const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
    const int site
):
    Base(tipPaths, alignedSeqs), 
    m_site(site - 1) { pruneTree(); }

bool Treemer::BySite::qualified(const clusters::iterator &candidate) const {
    auto it  = candidate->second.begin();
    char ref_value = (*it)->getSeq()[m_site];
    for (++it; it != candidate->second.end(); ++it) {
        if ((*it)->getSeq()[m_site] != ref_value) { return false; }
    }
    return true;
}
