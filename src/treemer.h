/*
 * Achieve a tip-to-root trimming behavior for a phylogenetic tree. The basic
 * idea is to start at the terminal node for each tree tip and move the node
 * towards the root node. Only when more than two nodes meet each other along
 * the way do the nodes make a real move otherwise the node stays.
 *
 * The tree is represented by tipPaths (from the R function "nodepath" in "ape"
 * package). The algorithm is performed on the tipPaths along with the
 * corresponding aligned sequences for each tip.
 */

#ifndef SITEPATH_TREEMER_H
#define SITEPATH_TREEMER_H

#include <map>
#include <string>
#include <utility>
#include <vector>
#include <Rcpp.h>

namespace Treemer {

float compare(const std::string &query, const std::string &subject);

/*
 * To hold all necessary data for a tree tip, including sequence, nodePath to
 * the tree root and the current node during treemer. Provides way to store
 * index of the 'kissed' node on the nodePath.
 */
class TipSeqLinker {
public:
    TipSeqLinker(
        const Rcpp::CharacterVector &sequence,
        const Rcpp::IntegerVector &tipPath
    );
    void reset();
    void proceed();
    int nextClade() const;
    int currentClade() const;
    int getTip() const;
    int getRoot() const;
    int getSeqLen() const;
    Rcpp::IntegerVector getPath() const;
    std::string getSeq() const;
    char siteChar(const int siteIndex) const;
private:
    const std::string m_seq;
    const Rcpp::IntegerVector m_path;
    const int m_tipIndex;
    int m_cIndex;
};

typedef std::vector<TipSeqLinker *> tips;
typedef std::map<int, tips> clusters;

/*
 * The Base class provides the structure for performing the trimming algorithm.
 * Depending on the extra constrain for cluster of tips, the class can be
 * derived to achieve a variety of behaviors.
 *
 * The extra constrain could be the largest similarity difference within a
 * cluster, the amino acid of the site and the average similarity.
 */
class Base {
public:
    Base(const tips &tips, const clusters initClusters);
    virtual ~Base();
    // Cluster of tree tips after trimming
    std::map< int, std::vector<int> > getTips() const;
    // Paths from root to trimming-node for all tips. Tips in the same cluster
    // will have the same path
    std::vector<Rcpp::IntegerVector> getPaths() const;
protected:
    // The actual trimming process happens here
    void pruneTree();
    // Whether the cluster of tips satisfy the extra constrain to be valid
    virtual bool qualified(const clusters::iterator &clusters_it) const = 0;
protected:
    // The container to hold all the tips
    const tips m_tips;
    // The clustering of tips during trimming process
    clusters m_clusters;
};

/*
 * Trim the tree by aa/nt of a site. The trimming stops when the non-dominant
 * aa/nt in a group is greater than the SNP percentage threshold
 */
class BySite: public Base {
public:
    BySite(
        const tips &tips,
        const clusters initClusters,
        const int siteIndex
    );
    // Cluster of TipSeqLinker grouped by aa/nt after trimming
    std::map<char, clusters> siteClusters() const;
private:
    bool qualified(const clusters::iterator &clusters_it) const;
private:
    const int m_siteIndex;
};

/*
 * Trim the tree and stop the process when the largest similarity difference in
 * each tip cluster is about to exceed the constrain
 */
class BySimilarity: public Base {
public:
    BySimilarity(
        const tips &tips,
        const clusters initClusters,
        const float simThreshold,
        std::map<std::pair<int, int>, float> &simMatrix
    );
protected:
    // The constrain of the largest similarity difference
    const float m_simCut;
    // Keep a record of pair-wise similarity to avoid repeating computation
    std::map<std::pair<int, int>, float> *m_compared;
private:
    bool qualified(const clusters::iterator &clusters_it) const;
};

}

#endif /* SITEPATH_TREEMER_H */
