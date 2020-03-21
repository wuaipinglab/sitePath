#include <algorithm>
#include "minEntropy.h"

template <class T>
MinEntropy::SearchTree<T>::SearchTree(
    const unsigned int minEffectiveSize,
    const unsigned int searchDepth,
    const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries
):
    m_minTipNum(minEffectiveSize),
    m_searchDepth(searchDepth),
    m_enclosed(nodeSummaries.size()) {
    // Transform R list to a vector of AA and mapped frequency as
    // 'm_aaSummaries'. Get all the possible segment points
    for (segIndex i = 0; i < m_enclosed; ++i) {
        m_all.push_back(i);
        Rcpp::IntegerVector summary = nodeSummaries[i].attr("aaSummary");
        Rcpp::CharacterVector aa = summary.names();
        aaSummary node;
        for (unsigned int j = 0; j < aa.size(); ++j) {
            node[Rcpp::as<std::string>(aa[j])] = summary[j];
        }
        m_aaSummaries.push_back(node);
    }
    // The segment point of "0" is removed
    m_all.erase(m_all.begin());
    // The final used list contains the last segment point to enclose
    // the last segment
    m_final.push_back(m_enclosed);
    // Generate the first parent node to initialize the search
    initSearch();
}

template <class T>
MinEntropy::SearchTree<T>::~SearchTree() {
    // Release the memory used by search nodes
    typedef typename std::vector<T *>::iterator iter;
    for (iter it = m_segList.begin(); it != m_segList.end(); ++it) {
        delete *it;
    }
    m_segList.clear();
}

template <>
void MinEntropy::SearchTree<MinEntropy::Segmentor>::initSearch() {
    // The adding search starts with only the last segment point
    // and all the rest are open for children nodes
    m_parent = new MinEntropy::Segmentor(
        m_all,
        m_final,
        m_aaSummaries,
        m_minTipNum
    );
    // Initialize minimum entropy as the initial parent's entropy
    m_minEntropy = m_parent->getEntropy();
}

template <>
void MinEntropy::SearchTree<MinEntropy::Amalgamator>::initSearch() {
    // Use the initial entropy is the starting parent node of the adding search
    // because the starting parent node of the removing search can be invalid
    MinEntropy::Segmentor noSeg(
            m_all,
            m_final,
            m_aaSummaries,
            m_minTipNum
    );
    m_final = noSeg.getUsed();
    m_minEntropy = noSeg.getEntropy();
    // Add the last segment point to the used list for the staring parent node
    m_all.push_back(m_enclosed);
    m_parent = new MinEntropy::Amalgamator(
        m_all,
        m_aaSummaries,
        m_minTipNum
    );
}

template <>
void MinEntropy::SearchTree<MinEntropy::Segmentor>::growTree(
        MinEntropy::Segmentor *seg
) {
    // Only the qualified search node in the adding search is added to
    // the active list. Deleted otherwise
    if (seg->isQualified()) {
        m_segList.push_back(seg);
    } else {
        delete seg;
    }
}

template <>
void MinEntropy::SearchTree<MinEntropy::Amalgamator>::growTree(
        MinEntropy::Amalgamator *seg
) {
    // The node will not be included in the active list and deleted if
    // a node with the same used list has already been evaluated.
    //
    // The node doesn't have to be qualified because we don't want to rule
    // out the possibly qualifed children node
    segment x = seg->getUsed();
    if (std::find(
            m_segListHistory.begin(),
            m_segListHistory.end(),
            x
    ) == m_segListHistory.end()) {
        m_segListHistory.push_back(x);
        m_segList.push_back(seg);
    } else {
        delete seg;
    }
}

template <class T>
MinEntropy::segment MinEntropy::SearchTree<T>::getFinal() const {
    return m_final;
}

template <class T>
float MinEntropy::SearchTree<T>::getMinEntropy() const {
    return m_minEntropy;
}

template <class T>
T *MinEntropy::SearchTree<T>::updateParent() {
    typedef typename std::vector<T *>::iterator iter;
    // Assume the first search node in the active list is the new parent node.
    iter it = m_segList.begin(), rm = it;
    // The new parent node has the minimum entropy among the active list
    T *tempMin = *it;
    for (++it; it != m_segList.end(); ++it) {
        if ((*it)->getEntropy() < tempMin->getEntropy()) {
            tempMin = *it;
            rm = it;
        }
    }
    // Remove the pointer of the new parent from the search list.
    m_segList.erase(rm);
    // Return the pointer to the new parent node
    return tempMin;
}

template <class T>
void MinEntropy::SearchTree<T>::resumeSearch() {
    if (!m_segList.empty()) {
        m_parent = updateParent();
        search();
    }
}

template <class T>
void MinEntropy::SearchTree<T>::search() {
    unsigned int depth = 0;
    const unsigned int maxDepth = m_enclosed * m_searchDepth;
    // Find the search node with minimum entropy in the active list and
    // make it the growing parent node for the next round. The current minimum
    // entropy should be decreasing but increasing is allowed. The used list
    // of a node is returned when its entropy stays minimum for a long time.
    while (true) {
        // Stop when the search reaches the end and delete the new parent node
        if (m_parent->isEndNode()) {
            delete m_parent;
            break;
        }
        // Generate children node from the parent node
        for (unsigned int i = 0; i < m_parent->getOpenSize(); ++i) {
            T *seg = new T(m_parent, i, m_aaSummaries, m_minTipNum);
            // This is to decide whether the child is valid and can be included
            // in the search list
            growTree(seg);
        }
        // Delete the parent node as its pointer already removed
        delete m_parent;
        // Stop when there is no search node in the list
        if (m_segList.empty()) { break; }
        // Get a new parent node from the search list.
        T *tempMin = updateParent();
        // Increment the number (depth) of consecutive times the candidate
        // node being minimum entropy or update to a new candidate node and
        // re-count the depth
        if (tempMin->getEntropy() > m_minEntropy) {
            // Increment the depth when the candidate is not beaten
            ++depth;
            // The search stops when the depth reaches the threshold
            if (depth >= maxDepth) { break; }
            // The candidate node stays unchanged if the new parent node
            // cannot beat the it.
        } else {
            // The new parent node will be the new candidate node
            // if the previous candidate is beaten by it.
            if (tempMin->isQualified()) {
                m_final = tempMin->getUsed();
                m_minEntropy = tempMin->getEntropy();
            }
            // Stop when the entropy of the new parent node is 0
            if (m_minEntropy == 0) { break; }
            // Re-count the depth for the new candidate
            depth = 0;
        }
        // Update the new parent node no matter whether it becomes the new
        // candidate or not.
        m_parent = tempMin;
    }
}

template class MinEntropy::SearchTree<MinEntropy::Segmentor>;
template class MinEntropy::SearchTree<MinEntropy::Amalgamator>;
