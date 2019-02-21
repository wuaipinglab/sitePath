#include "pruner.h"

TreeAlignmentMatch::TreeAlignmentMatch(
    const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
    const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
):
    m_root(*(tipPaths[0].begin())),
    m_seqLen((Rcpp::as<std::string>(alignedSeqs[0])).size()) 
{
    for (int i = 0; i < tipPaths.size(); i++) {
        TipSeqLinker *linker = new TipSeqLinker(alignedSeqs[i], tipPaths[i]);
        m_linkers.push_back(linker);
        m_clusters[linker->getTip()].push_back(linker);
        if (linker->getRoot() != m_root) {
            throw std::invalid_argument("Root in tree path not equal");
        } else if (linker->getSeqLen() != m_seqLen) {
            throw std::invalid_argument("Sequence length not equal");
        }
    }
}

TreeAlignmentMatch::~TreeAlignmentMatch()
{
    for(
        std::vector<TipSeqLinker*>::iterator it = m_linkers.begin();
        it != m_linkers.end(); ++it
    ) {
        delete (*it);
    }
}

std::map< int, std::vector<int> > TreeAlignmentMatch::getTips() const 
{
    std::map< int, std::vector<int> > tipCluster;
    for (
            std::vector<TipSeqLinker*>::const_iterator tsLinker = m_linkers.begin();
            tsLinker != m_linkers.end(); tsLinker++
    ) {
        tipCluster[(*tsLinker)->currentClade()].push_back((*tsLinker)->getTip());
    }
    return tipCluster;
}

std::vector<Rcpp::IntegerVector> TreeAlignmentMatch::getPaths() const 
{
    std::vector<Rcpp::IntegerVector> paths;
    for (
            std::vector<TipSeqLinker*>::const_iterator tsLinker = m_linkers.begin(); 
            tsLinker != m_linkers.end(); tsLinker++
    ) {
        paths.push_back((*tsLinker)->getPath());
    }
    return paths;
}

void TreeAlignmentMatch::pruneTree() 
{
    std::map< int, std::vector<TipSeqLinker*> > oldCluster;
    while (true) {
        for (
                std::vector<TipSeqLinker*>::const_iterator tsLinker = m_linkers.begin(); 
                tsLinker != m_linkers.end(); tsLinker++
        ) {
            m_trueCluster[(*tsLinker)->currentClade()].push_back(*tsLinker);
        }
        oldCluster = m_clusters;
        m_clusters.clear();
        // look down one more node (fake 'proceed') 
        // for each tip after new positioning
        for (
                std::vector<TipSeqLinker*>::iterator tsLinker = m_linkers.begin(); 
                tsLinker != m_linkers.end(); tsLinker++
        ) {
            //  group tips by fake 'proceed' node
            m_clusters[(*tsLinker)->nextClade()].push_back(*tsLinker);
        }
        // if no more group 'kissed' each other by a common ancestral node 
        // after fake 'proceed', then pruning is done
        if (m_clusters.size() == oldCluster.size()) {
            m_clusters.clear();
            break;
        }
        // only 'kissed' group can do real 'proceed'
        // if a grouping doesn't exist in 'oldCluster' 
        // then all tips in that group can 'proceed'
        for (
                std::map< int, std::vector<TipSeqLinker*> >::iterator it = m_clusters.begin();
                it != m_clusters.end(); it++
        ) {
            // assume a group is kissed with another (give it benefit of the doubt)
            bool kissed = true;
            for (
                    std::map< int, std::vector<TipSeqLinker*> >::iterator it2 = oldCluster.begin(); 
                    it2 != oldCluster.end(); it2++
            ) {
                // a group is 'non-kissed' after fake 'proceed' if it can be found in 'oldCluster'
                if (it->second == it2->second) {
                    kissed = false;
                    // a 'non-kissed' group won't appear twice in 'clusters' so deleted
                    oldCluster.erase(it2);
                    break;
                }
            }
            // candidate group needs to pass some requirement to be qualified 'kissed'
            if (kissed && qualified(it)) {
                for (
                        std::vector<TipSeqLinker*>::iterator tsLinker = it->second.begin(); 
                        tsLinker != it->second.end(); tsLinker++
                ) { (*tsLinker)->proceed(); }
            }
        }
        m_trueCluster.clear();
        oldCluster.clear();
    }
}

Pruner::Pruner(
    const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
    const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
    const float simThreshold,
    std::map<std::pair<int, int>, float> &simMatrix
):
    TreeAlignmentMatch(tipPaths, alignedSeqs),
    m_simCut(simThreshold),
    m_compared(simMatrix) 
{
    if (m_simCut <= 0) {
        throw std::invalid_argument("Similarity cannot be lower or equal to 0");
    } else if (m_simCut > 1) {
        throw std::invalid_argument("Similarity cannot be greater than 1");
    }
    if (m_simCut != 1) { pruneTree(); }
}

const bool Pruner::qualified(
        const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate
) {
    
    for (
            std::vector<TipSeqLinker*>::const_iterator tsLinker = candidate->second.begin(); 
            tsLinker != candidate->second.end() - 1; tsLinker++
    ) {
        for (
                std::vector<TipSeqLinker*>::const_iterator tsLinker2 = tsLinker + 1; 
                tsLinker2 != candidate->second.end(); tsLinker2++
        ) {
            std::pair<int, int> pairing = std::make_pair(
                (*tsLinker)->getTip(), 
                (*tsLinker2)->getTip()
            );
            float sim;
            if (m_compared.find(pairing) != m_compared.end()) {
                sim = m_compared[pairing];
            } else {
                sim = compare((*tsLinker)->getSeq(), (*tsLinker2)->getSeq());
                m_compared[pairing] = sim;
            }
            if (sim < m_simCut) {
                return false;
            }
        }
    }
    return true;
}

CustomizablePruner::CustomizablePruner(
    const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths,
    const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
    const std::map< int, std::vector<int> > &treeEdge,
    std::map<std::pair<int, int>, float> &simMatrix,
    const Rcpp::Function &customQualifyFunc
):
    TreeAlignmentMatch(tipPaths, alignedSeqs),
    m_nodeLink(treeEdge),
    m_compared(simMatrix),
    m_qualifyFunc(customQualifyFunc) { pruneTree(); }

const bool CustomizablePruner::qualified(
        const std::map< int, std::vector<TipSeqLinker*> >::iterator candidate
) {
    std::vector<float> crossSim, combinedSim;
    int nChildren = m_nodeLink[candidate->first].size();
    for (int i = 0; i < nChildren; ++i) {
        std::vector<TipSeqLinker*> *oldClstr = &m_trueCluster[m_nodeLink[candidate->first][i]];
        if (oldClstr->size() == 1) {goto CROSSCOMPARE;}
        for (
                std::vector<TipSeqLinker*>::iterator it1 = oldClstr->begin(); 
                it1 != oldClstr->end() - 1; ++it1
        ) {
            for (
                    std::vector<TipSeqLinker*>::iterator it2 = it1 + 1; 
                    it2 != oldClstr->end(); ++it2
            ) {
                combinedSim.push_back(m_compared[std::make_pair(
                        (*it1)->getTip(), (*it2)->getTip()
                )]);
            }
        }
        CROSSCOMPARE: for (int j = i + 1; j < nChildren; ++j) {
            std::vector<TipSeqLinker*> *oldClstr2 = &m_trueCluster[m_nodeLink[candidate->first][j]];
            for (
                    std::vector<TipSeqLinker*>::iterator it1 = oldClstr->begin(); 
                    it1 != oldClstr->end(); ++it1
            ) {
                for (
                        std::vector<TipSeqLinker*>::iterator it2 = oldClstr2->begin(); 
                        it2 != oldClstr2->end(); ++it2
                ) {
                    crossSim.push_back(m_compared[std::make_pair(
                            (*it1)->getTip(), (*it2)->getTip()
                    )]);
                }
            }
        }
    }
    if (crossSim.empty() || combinedSim.empty()) {
        return true;
    }
    return Rcpp::as<bool>(m_qualifyFunc(Rcpp::wrap(crossSim), Rcpp::wrap(combinedSim)));
}
