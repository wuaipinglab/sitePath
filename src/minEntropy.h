/*
 * Try to use branch and bound algorithm for finding the best segmenting
 * against tree tips along a lineage path. The basic idea is that a path
 * contains several segmenting points alongside. And by adding or removing
 * those points, the tree tips will be segmented into groups.
 *
 * By calculating Shannon Entropy, the purity of amino acid in all groups
 * amount to the total entropy. If using brute force, there is a max
 * number of 2^n calculations to find the minimum entropy.
 *
 * In branch and bound algorithm, each node represent a combination of
 * segment points. We start from no points and all points and converge
 * the result of two search trees.
 */

#ifndef SITEPATH_MINENTROPY_H
#define SITEPATH_MINENTROPY_H

#include "util.h"

namespace MinEntropy {

// The base class for node in tree search for minimum entropy
class TreeSearchNode {
public:
    // Pure virtual methods
    virtual unsigned int getOpenSize() const = 0;
    virtual bool isEndNode() const = 0;
    // Getters for private/protected member variables
    segment getUsed() const;
    float getEntropy() const;
    bool isQualified() const;
protected:
    // Empty constructor
    TreeSearchNode();
    // Initiate used segment points and calculate the entropy
    TreeSearchNode(
        const segment &used,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Calculate the total entropy
    float totalEntropy(
            const std::vector<aaSummary> &aaSummaries,
            const unsigned int minEffectiveSize
    );
protected:
    // The segment points being used and this list must
    // include the enclosed point
    segment m_used;
    // The Shannon Entropy of the search node
    float m_entropy;
    // Assume the search node is qualified
    bool m_qualified = true;
};

class Segmentor: public TreeSearchNode {
public:
    // Initially all segment points are open and
    // only terminal point is used.
    Segmentor(
        const segment &all,
        const segment &terminal,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Genetrat a child node from the parent node
    // The new used segment point is based on
    // open points in parent node.
    Segmentor(
        const Segmentor *parent,
        const unsigned int i,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Provide size of open points for generating children node
    unsigned int getOpenSize() const;
    // If the node has reach the end of the search
    bool isEndNode() const;
private:
    // The available segment points for children node
    segment m_open;
private:
    // The new used segment points
    segment newUsed(
            const Segmentor *parent,
            const unsigned int i
    ) const;
    // The available segment points for children
    segment newOpen(
            const Segmentor *parent,
            const unsigned int i
    ) const;
};

class Amalgamator: public TreeSearchNode {
public:
    // All segment points are being used initially
    Amalgamator(
        const segment &withTerminal,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Genetrat a child node from the parent node
    // The new used segment point is based on
    // used points in parent node.
    Amalgamator(
        const Amalgamator *parent,
        const unsigned int i,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Provide size of open points for generating children node
    unsigned int getOpenSize() const;
    // If the node has reach the end of the search
    bool isEndNode() const;
private:
    // The new used segment points
    segment newUsed(
            const Amalgamator *parent,
            const unsigned int i
    ) const;
};

template <class T>
class SearchTree {
public:
    SearchTree(
        const unsigned int minEffectiveSize,
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries
    );
    virtual ~SearchTree();
    // Pure virtual method for minimum entropy search
    void search();
    // Getter for final list of segment points
    segment getFinal() const;
private:
    // The minimum number of tips within a segmented group
    const unsigned int m_minTipNum;
    // The terminal segment point necessary for enclosing the segmenting
    const segIndex m_enclosed;
    // Store all possible segment points (except the enclosed point).
    // Track final list of segment points which gives minimum entropy
    segment m_all, m_final;
    // The transformed AA summaries for each node
    std::vector<aaSummary> m_aaSummaries;
    // Track the parent node
    T *m_parent;
    // To keep track of current minimum entropy of the segmenting
    float m_minEntropy;
    // The search list for resuming the search
    std::vector<T *> m_segList;
private:
    void initSearch();
    void growTree();
};

template<> void SearchTree<Segmentor>::initSearch();

template<> void SearchTree<Amalgamator>::initSearch();

template<> void SearchTree<Segmentor>::growTree();

template<> void SearchTree<Amalgamator>::growTree();

Rcpp::ListOf<Rcpp::IntegerVector> updatedSegmentation(
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries,
        const segment &final
);

}

#endif /* SITEPATH_MINENTROPY_H */
