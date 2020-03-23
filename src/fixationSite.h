#include <vector>
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
        const char siteChar,
        const int tipNum
    );
    char getSiteChar() const;
    int getTipNum() const;
    Rcpp::IntegerVector getPath() const;
    void updateDivIndices(std::vector<int> divIndices);
private:
    const Rcpp::IntegerVector m_path;
    std::vector<int> m_divIndices;
    const char m_siteChar;
    const int m_tipNum;
};

typedef std::vector<NodePath *> paths;

class TreeSearchNode {
public:
    TreeSearchNode(const paths &allPaths);
    float getScore() const;
private:
    float fixationScore() const;
    const paths m_selectedPaths;
};

class TreeSearch {
public:
    TreeSearch(const paths &allPaths);
    virtual ~TreeSearch();
    // sitePath getFinal() const;
    void search();
private:
    float m_bestScore;
    TreeSearchNode *m_parentNode;
    std::vector<TreeSearchNode *> m_searchNodes;
    const paths m_allPaths;
};

}
