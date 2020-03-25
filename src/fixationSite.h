/*
 * Branch and bound algorithm for finding path combination with minimum
 * entropy regarding one site.
 */

#ifndef SITEPATH_FIXATIONSITE_H
#define SITEPATH_FIXATIONSITE_H

#include <vector>
#include <utility>
#include <Rcpp.h>

namespace fixationSite {

class NodePath {
    /*
     * For each site, tip clustered by treemerBySite each have a nodePath.
     * The adjacent clusters sometimes can be combined but usually not due to
     * tree topology. They will be further clustered as predicted fixation.
     */
public:
    NodePath(
        const Rcpp::IntegerVector &path,
        const int tipNum,
        const char siteChar
    );
    // Getters
    Rcpp::IntegerVector getPaths() const;
    int getTipNum() const;
    char getSiteChar() const;
protected:
    const Rcpp::IntegerVector m_path;
    const int m_tipNum;
    const char m_siteChar;
};

typedef std::vector<NodePath *> mutPath;
typedef std::vector<mutPath> sitePath;
typedef std::pair<int, int> duoIndices;

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
    virtual ~TreeSearchNode();
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
    void fixationScore() const;
private:
    sitePath m_sitePath;
    float m_fixationScore;
};

class TreeSearch {
public:
    // TODO: not yet decided
    TreeSearch(const std::vector<NodePath *> &allNodePaths);
    virtual ~TreeSearch();
    // sitePath getFinal() const;
    void search();
private:
    const std::vector<NodePath *> m_allNodePaths;
    std::vector<TreeSearchNode *> m_searchNodes;
    const TreeSearchNode *m_parentNode;
    float m_bestScore;
};

}

#endif /* SITEPATH_FIXATIONSITE_H */
