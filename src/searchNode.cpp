#include <cmath>
#include <algorithm>
#include "minEntropy.h"

// Empty constructor for derived class
MinEntropy::TreeSearchNode::TreeSearchNode(): m_qualified(true) {}

MinEntropy::TreeSearchNode::~TreeSearchNode() {}

MinEntropy::TreeSearchNode::TreeSearchNode(
    const segment &used,
    const std::vector<aaSummary> &aaSummaries,
    const unsigned int minEffectiveSize
):
    m_used(used),
    m_qualified(true) {
    // Calculate the total entropy upon instantiating
    m_entropy = this->totalEntropy(aaSummaries, minEffectiveSize);
}

MinEntropy::segment MinEntropy::TreeSearchNode::getUsed() const {
    return m_used;
}

float MinEntropy::TreeSearchNode::getEntropy() const {
    return m_entropy;
}

bool MinEntropy::TreeSearchNode::isQualified() const {
    return m_qualified;
}

float MinEntropy::TreeSearchNode::getScore() const {
    return m_score;
}

float MinEntropy::TreeSearchNode::totalEntropy(
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
) {
    float res = 0.0;
    segIndex start = 0;
    // Iterate through all the used segment points
    for (
            segment::const_iterator m_used_itr = m_used.begin();
            m_used_itr != m_used.end(); ++m_used_itr
    ) {
        unsigned int tipNum = 0;
        aaSummary values;
        // Iterate throught the tree nodes between the segment points
        for (unsigned int i = start; i < *m_used_itr; ++i) {
            // Get the tree node and its summary on amino acid of tips
            // And combine them into a segment
            const aaSummary toBeCombined = aaSummaries.at(i);
            for (
                    aaSummary::const_iterator it = toBeCombined.begin();
                    it != toBeCombined.end(); ++it
            ) {
                values[it->first] += it->second;
                tipNum += it->second;
            }
        }
        // The search node will be disqualified if the number of tips
        // in the segment is lower than the constrain
        if (tipNum < minEffectiveSize) { m_qualified = false; }
        // TODO: tipNum can be used as total in calculating entropy
        res += shannonEntropy(values, tipNum);
        // Update the starting segment point
        start = *m_used_itr;
    }
    return res;
}

float MinEntropy::TreeSearchNode::fixationScore(
        const std::vector<aaSummary> &aaSummaries
) {
    float res = 0.0;
    segIndex start = 0;
    // Iterate through all the used segment points
    for (
            segment::const_iterator m_used_itr = m_used.begin();
            m_used_itr != m_used.end(); ++m_used_itr
    ) {
        unsigned int tipNum = 0;
        aaSummary values;
        // Iterate throught the tree nodes between the segment points
        for (unsigned int i = start; i < *m_used_itr; ++i) {
            // Get the tree node and its summary on amino acid of tips
            // And combine them into a segment
            const aaSummary toBeCombined = aaSummaries.at(i);
            for (
                    aaSummary::const_iterator it = toBeCombined.begin();
                    it != toBeCombined.end(); ++it
            ) {
                values[it->first] += it->second;
                tipNum += it->second;
            }
        }
        res += std::exp(shannonEntropy(values, tipNum)) / tipNum;
        // Update the starting segment point
        start = *m_used_itr;
    }
    return res;
}

MinEntropy::Segmentor::Segmentor(
    const segment &all,
    const segment &terminal,
    const std::vector<aaSummary> &aaSummaries,
    const unsigned int minEffectiveSize
):
    MinEntropy::TreeSearchNode(terminal, aaSummaries, minEffectiveSize),
    m_open(all) {}

MinEntropy::Segmentor::Segmentor(
    const Segmentor *parent,
    const unsigned int i,
    const std::vector<aaSummary> &aaSummaries,
    const unsigned int minEffectiveSize
):
    MinEntropy::TreeSearchNode() {
        // List-initialization is only supported since c++11
        m_used = this->newUsed(parent, i);
        m_open = this->newOpen(parent, i);
        m_entropy = this->totalEntropy(aaSummaries, minEffectiveSize);
    }

unsigned int MinEntropy::Segmentor::getOpenSize() const {
    // This is just for iterating the open list
    return m_open.size();
}

bool MinEntropy::Segmentor::isEndNode() const {
    // The search should be forced to end when the open list is empty
    return (m_open.empty()) ? true : false;
}

MinEntropy::segment MinEntropy::Segmentor::newUsed(
        const Segmentor *parent,
        const unsigned int i
) const {
    segment res = parent->m_used;
    // Add the "i"th segment point in the parent's open list to the parent's
    // used list
    res.push_back(parent->m_open.at(i));
    // Sort to make sure the segment points are in order
    std::sort(res.begin(), res.end());
    return res;
}

MinEntropy::segment MinEntropy::Segmentor::newOpen(
        const Segmentor *parent,
        const unsigned int i
) const {
    segment res = parent->m_open;
    // Erase the "i"th segment point in the parent's open list
    res.erase(res.begin() + i);
    return res;
}

MinEntropy::Amalgamator::Amalgamator(
    const segment &withTerminal,
    const std::vector<aaSummary> &aaSummaries,
    const unsigned int minEffectiveSize
):
    MinEntropy::TreeSearchNode(
        withTerminal,
        aaSummaries,
        minEffectiveSize
    ) {}

MinEntropy::Amalgamator::Amalgamator(
    const Amalgamator *parent,
    const unsigned int i,
    const std::vector<aaSummary> &aaSummaries,
    const unsigned int minEffectiveSize
):
    MinEntropy::TreeSearchNode() {
        // List-initialization is only supported since c++11
        m_used = this->newUsed(parent, i);
        m_entropy = this->totalEntropy(aaSummaries, minEffectiveSize);
    }

unsigned int MinEntropy::Amalgamator::getOpenSize() const {
    // This is just for iterating the open list
    return m_used.size() - 1;
}

bool MinEntropy::Amalgamator::isEndNode() const {
    // The search should be forced to end when the open list is empty
    return (m_used.size() == 1) ? true : false;
}

MinEntropy::segment MinEntropy::Amalgamator::newUsed(
        const Amalgamator *parent,
        const unsigned int i
) const {
    segment res = parent->getUsed();
    res.erase(res.begin() + i);
    // There is no need to sort as the order won't change
    return res;
}
