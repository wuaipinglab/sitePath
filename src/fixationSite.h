/*
 * Branch and bound algorithm for finding path combination with minimum
 * entropy regarding one site.
 */

#ifndef SITEPATH_FIXATIONSITE_H
#define SITEPATH_FIXATIONSITE_H

#include <vector>
#include <map>
#include <utility>
#include <Rcpp.h>

namespace fixationSite {

class MutPath {
    /*
     * For each site, tip clustered by treemerBySite each have a nodePath.
     * The adjacent clusters sometimes can be combined but usually not due to
     * tree topology. They will be further clustered as predicted fixation.
     */
public:
    MutPath(
        const std::vector<Rcpp::IntegerVector> &path,
        const std::map<char, int> siteSummary
    );
    // Getters
    std::map<char, int> getSiteSummary() const;
    std::vector<Rcpp::IntegerVector> getPaths() const;
    // Merge another NodePath
    void addPaths(const MutPath &other);
protected:
    std::vector<Rcpp::IntegerVector> m_path;
    const std::map<char, int> m_siteSummary;
};

typedef std::vector<MutPath *> sitePath;
typedef std::pair<MutPath *, MutPath *> pairPaths;

class TreeSearchNode {
    /*
     * The node in the search tree, storing the combination of paths.
     * Each node will bind two paths (or combined paths) from the parent.
     * But only the adjacent paths (monophyletic) can be combined.
     */
public:
    // The starting node of the search tree
    TreeSearchNode(const sitePath &allPaths);
    // The node in the progression of search
    TreeSearchNode(
        const TreeSearchNode &parent,
        const pairPaths &toCombine
    );
    virtual ~TreeSearchNode();
    // Getters
    float getScore() const;
    sitePath getCombinedPaths() const;
    // For children node of the search
    std::vector<pairPaths> allPairPaths() const;
private:
    // The starting node
    void pathsCombine(const sitePath &allPaths);
    // The node in progression of search
    void pathsCombine(
            const TreeSearchNode &parent,
            const pairPaths &toCombine
    );
    // Calculate the score
    void fixationScore() const;
private:
    sitePath m_combinedPaths;
    float m_fixationScore;
};

class TreeSearch {
public:
    // TODO: not yet decided
    TreeSearch(const sitePath &allPaths);
    virtual ~TreeSearch();
    // sitePath getFinal() const;
    void search();
private:
    std::vector<TreeSearchNode *> m_searchNodes;
    const TreeSearchNode *m_parentNode;
    float m_bestScore;
};

}

#endif /* SITEPATH_FIXATIONSITE_H */
