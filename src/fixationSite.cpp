#include <map>
#include <cmath>
#include "fixationSite.h"

FixationSite::NodePath::NodePath(
    const Rcpp::IntegerVector &path,
    const std::vector<int> &tips,
    const char siteChar
):
    m_path(path),
    m_tips(tips),
    m_siteChar(siteChar) {
    m_monophyleticNode = m_path[m_path.size() - 2];
    m_paraphyleticNode = m_path.size() >= 3? m_path[m_path.size() - 2]: m_monophyleticNode;
}

Rcpp::IntegerVector FixationSite::NodePath::getPaths() const {
    return m_path;
}

std::vector<int> FixationSite::NodePath::getTips() const {
    return m_tips;
}

char FixationSite::NodePath::getSiteChar() const {
    return m_siteChar;
}

int FixationSite::NodePath::tipNum() const {
    return m_tips.size();
}

int FixationSite::NodePath::getMonophyleticNode() const {
    return m_monophyleticNode;
}

int FixationSite::NodePath::getParaphyleticNode() const {
    return m_paraphyleticNode;
}

FixationSite::TreeSearchNode::TreeSearchNode(
    const std::vector<NodePath *> &allNodePaths
) {
    pathsCombine(allNodePaths);
    fixationScore();
}

FixationSite::TreeSearchNode::TreeSearchNode(
    const TreeSearchNode &parent,
    const FixationSite::duoIndices &toCombine
) {
    pathsCombine(parent, toCombine);
    fixationScore();
}

FixationSite::sitePath FixationSite::TreeSearchNode::getSitePath() const {
    return m_sitePath;
}

float FixationSite::TreeSearchNode::getScore() const {
    return m_fixationScore;
}

std::vector<FixationSite::duoIndices> FixationSite::TreeSearchNode::allToCombine() const {
    // All possible pairs of adjacent paraphyletic clades
    std::vector<duoIndices> res;
    // Iterate all grouped paths
    for (unsigned int first_i = 0; first_i < m_sitePath.size() - 1; first_i++) {
        // All mono- and para-phyletic nodes of the first grouped paths
        std::vector<int> monoNodes;
        std::vector<int> paraNodes;
        for (
                mutPath::const_iterator mp_itr = m_sitePath[first_i].begin();
                mp_itr != m_sitePath[first_i].end(); mp_itr++
        ) {
            monoNodes.push_back((**mp_itr).getMonophyleticNode());
            paraNodes.push_back((**mp_itr).getParaphyleticNode());
        }
        // Iterate to get all possible pairs of grouped paths
        for (unsigned int second_i = first_i+1; second_i < m_sitePath.size(); second_i++) {
            // If any monophyletic node of the second grouped paths can
            // be found in the mono- and para-phyletic of the first group
            for (
                    mutPath::const_iterator mp_itr = m_sitePath[second_i].begin();
                    mp_itr != m_sitePath[second_i].end(); mp_itr++
            ) {
                bool found = false;
                for (
                        std::vector<int>::iterator nodes_itr = monoNodes.begin();
                        nodes_itr != monoNodes.end(); nodes_itr++
                ) {
                    if (*nodes_itr == (**mp_itr).getParaphyleticNode()) {
                        found = true;
                        res.push_back(std::make_pair(first_i, second_i));
                        break;
                    }
                }
                for (
                        std::vector<int>::iterator nodes_itr = paraNodes.begin();
                        nodes_itr != paraNodes.end(); nodes_itr++
                ) {
                    if (*nodes_itr == (**mp_itr).getMonophyleticNode()) {
                        found = true;
                        res.push_back(std::make_pair(second_i, first_i));
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
        }
    }
    return res;
}

void FixationSite::TreeSearchNode::pathsCombine(
        const std::vector<NodePath *> &allNodePaths
) {
    for (
            std::vector<NodePath *>::const_iterator it = allNodePaths.begin();
            it != allNodePaths.end(); it++
    ) {
        std::vector<NodePath *> toAdd;
        toAdd.push_back(*it);
        m_sitePath.push_back(toAdd);
    }
}

void FixationSite::TreeSearchNode::pathsCombine(
        const TreeSearchNode &parent,
        const duoIndices &toCombine
) {
    m_sitePath = parent.getSitePath();
    m_sitePath[toCombine.first].insert(
            m_sitePath[toCombine.first].begin(),
            m_sitePath[toCombine.second].begin(),
            m_sitePath[toCombine.second].end()
    );
    m_sitePath[toCombine.second] = m_sitePath.back();
    m_sitePath.pop_back();
}

void FixationSite::TreeSearchNode::fixationScore() {
    // The score equals to the exponential of sum of each group's
    // Shannon Entropy times total number of groups
    float totalScore = static_cast<float>(m_sitePath.size());
    float entropyExpInv = 1.0;
    for (
            sitePath::const_iterator sp_itr = m_sitePath.begin();
            sp_itr != m_sitePath.end(); sp_itr++
    ) {
        std::map<char, int> siteSummary;
        int totalNum = 0;
        for (
            mutPath::const_iterator mp_itr = sp_itr->begin();
            mp_itr != sp_itr->end(); mp_itr++
        ) {
            siteSummary[(*mp_itr)->getSiteChar()] += (*mp_itr)->tipNum();
            totalNum += (*mp_itr)->tipNum();
        }
        for (
                std::map<char, int>::iterator it = siteSummary.begin();
                it != siteSummary.end(); it++
        ) {
            float p = it->second/static_cast<float>(totalNum);
            entropyExpInv *= std::pow(p, p);
        }
    }
    m_fixationScore = totalScore/entropyExpInv;
}

FixationSite::TreeSearch::TreeSearch(
    const Treemer::clusters &grouping,
    const int site
) {
    initSearch(grouping, site - 1);
    search();
}

FixationSite::TreeSearch::~TreeSearch() {
    for (
            std::vector<NodePath *>::const_iterator it = m_allNodePaths.begin();
            it != m_allNodePaths.end(); it++
    ) {
        delete *it;
    }
    for (
            std::vector<TreeSearchNode *>::const_iterator it = m_searchNodes.begin();
            it != m_searchNodes.end(); it++
    ) {
        delete *it;
    }
}

FixationSite::sitePath FixationSite::TreeSearch::getFinal() const {
    return m_bestSitePath;
}

void FixationSite::TreeSearch::initSearch(
        const Treemer::clusters &grouping,
        const int siteIndex
) {
    for (
            Treemer::clusters::const_iterator clusters_itr = grouping.begin();
            clusters_itr != grouping.end(); clusters_itr++
    ) {
        std::vector<int> tipNodes;
        for (
                Treemer::tips::const_iterator tips_itr = clusters_itr->second.begin();
                tips_itr != clusters_itr->second.end(); tips_itr++
        ) {
            tipNodes.push_back((*tips_itr)->getTip());
        }
        NodePath *p = new NodePath(
            clusters_itr->second[0]->getPath(),
            tipNodes,
            clusters_itr->second[0]->getSiteChar(siteIndex)
        );
        m_allNodePaths.push_back(p);
    }
}

void FixationSite::TreeSearch::search() {
    m_parentNode = new TreeSearchNode(m_allNodePaths);
    m_bestSitePath = m_parentNode->getSitePath();
    m_bestScore = m_parentNode->getScore();
    while (true) {
        const std::vector<duoIndices> allToCombine = m_parentNode->allToCombine();
        // Create children nodes from the parent node and add to the search list
        for (
                std::vector<duoIndices>::const_iterator it = allToCombine.begin();
                it != allToCombine.end(); it++
        ) {
            TreeSearchNode *child = new TreeSearchNode(*m_parentNode, *it);
            m_searchNodes.push_back(child);
        }
        // Free the root search node or if during the search, the old parent's
        // pointer already removed from search list and needs to be freed
        delete m_parentNode;
        std::vector<TreeSearchNode *>::iterator
            searchNodes_itr = m_searchNodes.begin(), toRemove = searchNodes_itr;
        // Assume the first search node in the search list is the new parent node.
        m_parentNode = *searchNodes_itr;
        // Update the search node with better (less) fixationScore as new parent
        // Remove the new parent from the search list
        for (++searchNodes_itr; searchNodes_itr != m_searchNodes.end(); searchNodes_itr++) {
            if ((*searchNodes_itr)->getScore() < m_parentNode->getScore()) {
                m_parentNode = *searchNodes_itr;
                toRemove = searchNodes_itr;
            }
        }
        m_searchNodes.erase(toRemove);
        // Terminate when all groups are combined or the new parent is worse than the old parent
        if (m_parentNode->getSitePath().size() == 1 || m_parentNode->getScore() > m_bestScore) {
            delete m_parentNode;
            break;
        } else {
            m_bestSitePath = m_parentNode->getSitePath();
            m_bestScore = m_parentNode->getScore();
        }
    }
}
