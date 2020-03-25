#include <map>
#include <cmath>
#include "fixationSite.h"

fixationSite::NodePath::NodePath(
    const Rcpp::IntegerVector &path,
    const int tipNum,
    const char siteChar
):
    m_path(path),
    m_tipNum(tipNum),
    m_siteChar(siteChar) {}

Rcpp::IntegerVector fixationSite::NodePath::getPaths() const {
    return m_path;
}

int fixationSite::NodePath::getTipNum() const {
    return m_tipNum;
}

char fixationSite::NodePath::getSiteChar() const {
    return m_siteChar;
}

fixationSite::TreeSearchNode::TreeSearchNode(
    const std::vector<NodePath *> &allNodePaths
) {
    pathsCombine(allNodePaths);
    fixationScore();
}

fixationSite::TreeSearchNode::TreeSearchNode(
    const TreeSearchNode &parent,
    const fixationSite::duoIndices &toCombine
) {
    pathsCombine(parent, toCombine);
    fixationScore();
}

fixationSite::sitePath fixationSite::TreeSearchNode::getSitePath() const {
    return m_sitePath;
}

float fixationSite::TreeSearchNode::getScore() const {
    return m_fixationScore;
}

std::vector<fixationSite::duoIndices> fixationSite::TreeSearchNode::allToCombine() const {
    std::vector<fixationSite::duoIndices> res;
    // TODO: make it work
    return res;
}

void fixationSite::TreeSearchNode::pathsCombine(
        const std::vector<fixationSite::NodePath *> &allNodePaths
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

void fixationSite::TreeSearchNode::pathsCombine(
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

void fixationSite::TreeSearchNode::fixationScore() const {
    float totalScore = 0.0;
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
            siteSummary[(*mp_itr)->getSiteChar()] += (*mp_itr)->getTipNum();
            totalNum += (*mp_itr)->getTipNum();
        }
        for (
                std::map<char, int>::iterator it = siteSummary.begin();
                it != siteSummary.end(); it++
        ) {}
    }
}

// TODO: SearchTree implementation
