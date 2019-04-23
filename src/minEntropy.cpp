#include <cmath>
#include "minEntropy.h"

MinEntropy::TreeSearchNode::TreeSearchNode() {}

MinEntropy::TreeSearchNode::TreeSearchNode(
    const segment &used,
    const std::vector<aaSummary> &aaSummaries,
    const unsigned int minEffectiveSize
):
    m_used(used) {
    m_entropy = this->totalEntropy(aaSummaries, minEffectiveSize);
}

segment MinEntropy::TreeSearchNode::getUsed() const { return m_used; }

float MinEntropy::TreeSearchNode::getEntropy() const { return m_entropy; }

bool MinEntropy::TreeSearchNode::isQualified() const { return m_qualified; }

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
            const aaSummary toBeCombined = aaSummaries.at(i);
            for (
                    aaSummary::const_iterator it = toBeCombined.begin();
                    it != toBeCombined.end(); ++it
            ) {
                values[it->first] += it->second;
                tipNum += it->second;
            }
        }
        if (tipNum < minEffectiveSize) { m_qualified = false; }
        // TODO: tipNum can be used as total in calculating entropy
        res += shannonEntropy(values);
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
        m_used = this->newUsed(parent, i);
        m_open = this->newOpen(parent, i);
        m_entropy = this->totalEntropy(aaSummaries, minEffectiveSize);
    }

unsigned int MinEntropy::Segmentor::getOpenSize() const {
    return m_open.size();
}

bool MinEntropy::Segmentor::isEndNode() const {
    return (m_open.empty()) ? true : false;
}

segment MinEntropy::Segmentor::newUsed(
        const Segmentor *parent,
        const unsigned int i
) const {
    segment res = parent->m_used;
    res.push_back(parent->m_open.at(i));
    // Sort the "res" to make sure the segment points
    // are in order
    std::sort(res.begin(), res.end());
    return res;
}

segment MinEntropy::Segmentor::newOpen(
        const Segmentor *parent,
        const unsigned int i
) const {
    segment res = parent->m_open;
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
        m_used = this->newUsed(parent, i);
        m_entropy = this->totalEntropy(aaSummaries, minEffectiveSize);
    }

unsigned int MinEntropy::Amalgamator::getOpenSize() const {
    return m_used.size() - 1;
}

bool MinEntropy::Amalgamator::isEndNode() const {
    return (m_used.size() == 1) ? true : false;
}

segment MinEntropy::Amalgamator::newUsed(
        const Amalgamator *parent,
        const unsigned int i
) const {
    segment res = parent->getUsed();
    res.erase(res.begin() + i);
    // There is no need to sort as the order won't change
    // And no need to create the list for available segment
    // points as "m_used" here already contains them.
    return res;
}

template <class T>
MinEntropy::SearchTree<T>::SearchTree(
    const unsigned int minEffectiveSize,
    const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries
):
    m_minTipNum(minEffectiveSize),
    m_enclosed(nodeSummaries.size())
{
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
    // The segment point of "0" should be removed
    m_all.erase(m_all.begin());
    // The "m_final" must have "m_enclosed" as the last element
    m_final.push_back(m_enclosed);
    // Generate the first parent node for the search
    initSearch();
}

template <class T>
MinEntropy::SearchTree<T>::~SearchTree() {
    typedef typename std::vector<T *>::iterator iter;
    for (iter it = m_segList.begin(); it != m_segList.end(); ++it) {
        delete *it;
    }
    m_segList.clear();
}

template<>
void MinEntropy::SearchTree<MinEntropy::Segmentor>::initSearch() {
    // The search by adding segment points starts with only "m_enclosed"
    // as used and "m_all" open for children node
    m_parent = new MinEntropy::Segmentor(
        m_all,
        m_final,
        m_aaSummaries,
        m_minTipNum
    );
    // Initialize minimum entropy as initial parent's entropy
    m_minEntropy = m_parent->getEntropy();
}

template<>
void MinEntropy::SearchTree<MinEntropy::Amalgamator>::initSearch() {
    // Use the "m_enclosed" segment point and its entropy as the
    // initial values as the initial "Amalgamator" can be unqualified
    MinEntropy::Segmentor noSeg(
            m_all,
            m_final,
            m_aaSummaries,
            m_minTipNum
    );
    m_final = noSeg.getUsed();
    m_minEntropy = noSeg.getEntropy();
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
    segment x = seg->getUsed();
    if (std::find(
            m_segListHistory.begin(),
            m_segListHistory.end(),
            x
    ) == m_segListHistory.end()) {
        m_segListHistory.push_back(x);
        m_segList.push_back(seg);
    }
}

template <class T>
segment MinEntropy::SearchTree<T>::getFinal() const { return m_final; }

template <class T>
float MinEntropy::SearchTree<T>::getMinEntropy() const { return m_minEntropy; }

template <class T>
T *MinEntropy::SearchTree<T>::updateParent() {
    typedef typename std::vector<T *>::iterator iter;
    // Assume the first "SearchNode" is the best hit.
    iter it = m_segList.begin(), rm = it;
    // This "tempMin" is always one of the "SearchNode" in the search
    // list. And once it's found, it's going to be the new "parent"
    // and its pointer is removed from the search list.
    T *tempMin = *it;
    for (++it; it != m_segList.end(); ++it) {
        if ((*it)->getEntropy() < tempMin->getEntropy()) {
            tempMin = *it;
            rm = it;
        }
    }
    m_segList.erase(rm);
    // Return the pointer to the current best hit in the "SearchNode" list
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
    unsigned int maxDepth = m_enclosed;
    // Find the "minEntropy" among "SearchNode" in "segList" and make
    // the "SearchNode" the "parent" for next round. The new "minEntropy"
    // should be decreasing but increasing is allowed. So "final" is
    // returned when its "minEntropy" stays for "m_enclosed" consecutive
    // of times.
    while (true) {
        // for (int i = 0; i < m_parent->getUsed().size(); ++i) {
        //     std::cout << m_parent->getUsed()[i] << " ";
        // }
        // std::cout << std::endl;
        // Stop when the search reaches the end and delete "m_parent"
        if (m_parent->isEndNode()) {
            delete m_parent;
            break;
        }
        // Children "SearchNode" of the current "parent"
        for (unsigned int i = 0; i < m_parent->getOpenSize(); ++i) {
            T *seg = new T(m_parent, i, m_aaSummaries, m_minTipNum);
            // This is to decide whether the child is included
            // in the search list
            growTree(seg);
        }
        // delete "parent" as its pointer already removed from "m_segList"
        delete m_parent;
        // Stop when there is no search node in the list
        if (m_segList.empty()) { break; }
        // Re-calculate the best hit from the search list and assign
        // to "tempMin".
        T *tempMin = updateParent();
        // Increment the number ("depth") of consecutive times last
        // "minEntropy" being best hit or update "minEntropy" and
        // re-count "depth"
        if (tempMin->getEntropy() > m_minEntropy) {
            // Increment the 'consecutive times'
            ++depth;
            // The search stops when none of the children "SearchNode"
            // can beat the "minEntropy"
            if (depth >= maxDepth) { break; }
            // "final" and "minEntropy" stay unchanged if "tempMin"
            // cannot beat the previous "minEntropy".
        } else {
            // "tempMin" will be the new "final" and "minEntropy"
            // if the previous "minEntropy" is beaten by it.
            if (tempMin->isQualified()) {
                m_final = tempMin->getUsed();
                m_minEntropy = tempMin->getEntropy();
            }
            // Stop when entropy is 0, which couldn't be better
            if (m_minEntropy == 0) { break; }
            // "depth" is reset to 0
            depth = 0;
        }
        // "tempMin" is going to be the new "parent" for getting new
        // children "SearchNode" no matter what. So "tempMin" doesn't
        // necessarily have to give the new "final" and "minEntropy"
        m_parent = tempMin;
    }
}

Rcpp::ListOf<Rcpp::IntegerVector> MinEntropy::updatedSegmentation(
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries,
        const segment &final
) {
    // Group tips by the segment indices in "final"
    std::vector<Rcpp::IntegerVector> res;
    segIndex start = 0;
    std::string prevFixedAA = "";
    aaSummary combNode;
    std::vector<int> combTips;
    for (
            segment::const_iterator final_itr = final.begin();
            final_itr != final.end(); ++final_itr
    ) {
        aaSummary node;
        std::vector<int> tips;
        for (unsigned int i = start; i < *final_itr; ++i) {
            Rcpp::IntegerVector nodeTips = nodeSummaries[i];
            tips.insert(tips.end(), nodeTips.begin(), nodeTips.end());
            Rcpp::IntegerVector summary = nodeTips.attr("aaSummary");
            Rcpp::CharacterVector aa = summary.names();
            for (unsigned int j = 0; j < aa.size(); ++j) {
                node[Rcpp::as<std::string>(aa.at(j))] += summary.at(j);
            }
        }
        // The most dominant AA in the current segment
        aaSummary::iterator node_itr = node.begin();
        std::string fixedAA = node_itr->first;
        int maxFreq = node_itr->second;
        for (++node_itr; node_itr != node.end(); ++node_itr) {
            if (node_itr->second > maxFreq) {
                fixedAA = node_itr->first;
                maxFreq = node_itr->second;
            }
        }
        // Combine with the previous segment if the most dominant
        // AA is the same as the previous.
        if (fixedAA == prevFixedAA) {
            for (node_itr = node.begin(); node_itr != node.end(); ++node_itr) {
                combNode[node_itr->first] += node_itr->second;
            }
            combTips.insert(combTips.end(), tips.begin(), tips.end());
        } else {
            Rcpp::IntegerVector combined = Rcpp::wrap(combTips);
            combined.attr("aaSummary") = Rcpp::wrap(combNode);
            combined.attr("AA") = prevFixedAA;
            res.push_back(combined);
            combTips = tips;
            combNode = node;
        }
        start = *final_itr;
        prevFixedAA = fixedAA;
    }
    Rcpp::IntegerVector combined = Rcpp::wrap(combTips);
    combined.attr("aaSummary") = Rcpp::wrap(combNode);
    combined.attr("AA") = prevFixedAA;
    res.push_back(combined);
    res.erase(res.begin());
    return Rcpp::wrap(res);
}

template class MinEntropy::SearchTree<MinEntropy::Segmentor>;
template class MinEntropy::SearchTree<MinEntropy::Amalgamator>;
