/*
 * Branch and bound algorithm for finding path combination with minimum
 * entropy regarding one site.
 */

#ifndef SITEPATH_FIXATIONSITE_H
#define SITEPATH_FIXATIONSITE_H

#include <vector>
#include <utility>
#include <Rcpp.h>
#include "treemer.h"

namespace FixationSite {

class NodePath {
    /*
     * For each site, tip clustered by treemerBySite each have a nodePath.
     * The adjacent clusters sometimes can be combined but usually not due to
     * tree topology. They will be further clustered as predicted fixation.
     */
public:
    NodePath(
        const Rcpp::IntegerVector &path,
        const std::vector<int> &tips,
        const char siteChar
    );
    // Getters
    Rcpp::IntegerVector getPaths() const;
    std::vector<int> getTips() const;
    char getSiteChar() const;
    int getMonophyleticNode() const;
    int getParaphyleticNode() const;
    int tipNum() const;
protected:
    const Rcpp::IntegerVector m_path;
    const std::vector<int> m_tips;
    int m_monophyleticNode;
    int m_paraphyleticNode;
    const char m_siteChar;
};

typedef std::pair<int, int> duoIndices;

class FixedTips {
public:
    FixedTips(const NodePath *nodePath);
    // Getters
    const std::vector<NodePath *> *getNodePaths() const;
    std::map<char, int> getSiteSummary() const;
    std::vector<int> getMonophyleticNodes() const;
    std::vector<int> getParaphyleticNodes() const;
    // Functions
    void merge(const FixedTips &other);
    std::vector<int> allTips() const;
private:
    std::vector<const NodePath *> m_nodePaths;
    std::map<char, int> m_siteSummary;
    std::vector<int> m_monophyleticNodes;
    std::vector<int> m_paraphyleticNodes;
};

typedef std::vector<NodePath *> mutPath;
typedef std::vector<mutPath> sitePath;

class TreeSearchNode {
    /*
     * The node in the search tree, storing the combination of paths.
     * Each node will bind two paths (or combined paths) from the parent.
     * But only the adjacent paths (paraphyletic) can be combined.
     */
public:
    // The starting node of the search tree
    TreeSearchNode(const std::vector<NodePath *> &allNodePaths);
    // The node in the progression of search
    TreeSearchNode(
        const TreeSearchNode &parent,
        const duoIndices &toCombine
    );
    // Getters
    sitePath getSitePath() const;
    float getScore() const;
    // For children node of the search
    std::vector<duoIndices> allToCombine() const;
private:
    // The starting node
    void pathsCombine(const std::vector<NodePath *> &allNodePaths);
    // The node in progression of search
    void pathsCombine(
            const TreeSearchNode &parent,
            const duoIndices &toCombine
    );
    // Calculate the score
    void fixationScore();
private:
    sitePath m_sitePath;
    float m_fixationScore;
};

class TreeSearch {
    /*
     * The branch and bound method for finding the best sitePath
     */
public:
    TreeSearch(
        const Treemer::clusters &grouping,
        const int site
    );
    virtual ~TreeSearch();
    sitePath getFinal() const;
    void initSearch(
            const Treemer::clusters &grouping,
            const int siteIndex
    );
    void search();
private:
    std::vector<NodePath *> m_allNodePaths;
    sitePath m_bestSitePath;
    std::vector<TreeSearchNode *> m_searchNodes;
    const TreeSearchNode *m_parentNode;
    float m_bestScore;
};

}

#endif /* SITEPATH_FIXATIONSITE_H */
