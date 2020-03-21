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

#include <string>
#include <vector>
#include <map>
#include <Rcpp.h>

namespace MinEntropy {

typedef std::map<std::string, int> aaSummary;
typedef unsigned int segIndex;
typedef std::vector<segIndex> segment;

float shannonEntropy(const aaSummary &values, const unsigned int tipNum);

Rcpp::ListOf<Rcpp::IntegerVector> updatedSegmentation(
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries,
        const segment &final
);

class TreeSearchNode {
    /*
     * The base class for implementing node in tree search for minimum entropy
     * The node stores the segment points and calculates the total entropy of
     * the current segmentation
     */
public:
    // Pure virtual methods
    virtual unsigned int getOpenSize() const = 0;
    virtual bool isEndNode() const = 0;
    virtual ~TreeSearchNode();
    // Getters for private/protected member variables
    segment getUsed() const;
    float getEntropy() const;
    bool isQualified() const;
    float getScore() const;
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
    // Fixation score integrating entropy and tip numbers
    float fixationScore(
            const std::vector<aaSummary> &aaSummaries
    );
protected:
    // A list of the segment points which must include the enclosed point
    segment m_used;
    // The Shannon Entropy of the search node
    float m_entropy;
    // Whether the search node is qualified (minEffectiveSize)
    bool m_qualified;
    // The fixation score integrating entropy and tip numbers
    float m_score;
};

class Segmentor: public TreeSearchNode {
    /*
     * The class to implement the node in the adding search
     *
     * Store segment indices of a nodePath and the remaining unused indices
     * The growing of the tree search is by adding one segment point from
     * the open list to the used list for each children node
     */
public:
    // Initially only the terminal enclosed point is used.
    Segmentor(
        const segment &all,
        const segment &terminal,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Genetrat a children node from the parent node
    // Simply pick the "i"th segment point from the parent's open list
    Segmentor(
        const Segmentor *parent,
        const unsigned int i,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Provide size of open list for generating children node
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
    /*
     * The class to implemenet the node in the removing search
     *
     * Only the used list is needed tracking as it's also the open list.
     * The growing of the search tree is by removing one segment point
     * from the used list
     */
public:
    // All segment points are being used initially
    Amalgamator(
        const segment &withTerminal,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Genetrat a child node from the parent node
    Amalgamator(
        const Amalgamator *parent,
        const unsigned int i,
        const std::vector<aaSummary> &aaSummaries,
        const unsigned int minEffectiveSize
    );
    // Provide size of used list for generating children node
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
    /*
     * The template class for implementing the tree search
     *
     * Store the original nodePath segmentation and search constrain.
     * It's gonna carry the heuristic search for minimum entropy
     */
public:
    SearchTree(
        const unsigned int minEffectiveSize,
        const unsigned int searchDepth,
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries
    );
    virtual ~SearchTree();
    // Getters for the private/protected member variables
    segment getFinal() const;
    float getMinEntropy() const;
    // Minimum entropy search
    void search();
    void resumeSearch();
private:
    // The minimum number of tips within a segmented group
    const unsigned int m_minTipNum, m_searchDepth;
    // The terminal segment point essential for enclosing the segmenting
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
    // The search list containing all the active search nodes
    std::vector<T *> m_segList;
    // The existing or dropped segment
    std::vector<segment> m_segListHistory;
private:
    // Initialize the tree search
    void initSearch();
    // Grow the tree from the current parent node
    void growTree(T *seg);
    void updateFinal(T *tempMin);
    // Decide the new parent node
    T *updateParent();
};

template <> void SearchTree<Segmentor>::initSearch();
template <> void SearchTree<Amalgamator>::initSearch();

template <> void SearchTree<Segmentor>::growTree(Segmentor *seg);
template <> void SearchTree<Amalgamator>::growTree(Amalgamator *seg);

// Not yet implemented
template <> void SearchTree<Segmentor>::updateFinal(Segmentor *tempMin);
template <> void SearchTree<Amalgamator>::updateFinal(Amalgamator *tempMin);

}

#endif /* SITEPATH_MINENTROPY_H */
