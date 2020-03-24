#include "minEntropy.h"
#include "fixationSite.h"

fixationSite::MutPath::MutPath(
    const std::vector<Rcpp::IntegerVector> &path,
    const std::map<char, int> siteSummary
):
    m_path(path),
    m_siteSummary(siteSummary) {}

std::map<char, int> fixationSite::MutPath::getSiteSummary() const {
    return m_siteSummary;
}

std::vector<Rcpp::IntegerVector> fixationSite::MutPath::getPaths() const {
    return m_path;
}

void fixationSite::MutPath::addPaths(const MutPath &other) {
    // TODO: make it work
}

fixationSite::TreeSearchNode::TreeSearchNode(
    const sitePath &paths
):
    m_combinedPaths(paths) {
    fixationScore();
}

fixationSite::TreeSearchNode::TreeSearchNode(
    const TreeSearchNode &parent,
    const fixationSite::pairPaths &toCombine
) {
    pathsCombine(parent, toCombine);
    fixationScore();
}

fixationSite::TreeSearchNode::~TreeSearchNode() {
    for (
            sitePath::iterator it = m_combinedPaths.begin();
            it != m_combinedPaths.end(); it++
    ) {
        delete *it;
    }
}

float fixationSite::TreeSearchNode::getScore() const {
    return m_fixationScore;
}

fixationSite::sitePath fixationSite::TreeSearchNode::getCombinedPaths() const {
    return m_combinedPaths;
}

std::vector<fixationSite::pairPaths> fixationSite::TreeSearchNode::allPairPaths() const {
    std::vector<fixationSite::pairPaths> res;
    // TODO: make it work
    return res;
}

void fixationSite::TreeSearchNode::pathsCombine(
        const TreeSearchNode &parent,
        const pairPaths &toCombine
) {
    // m_combinedPaths.push_back();
    // const sitePath parentPaths = parent.getCombinedPaths();
}


// TODO: SearchTree implementation
