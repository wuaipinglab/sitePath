#include "treemer.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix getSimilarityMatrix(
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs
) {
    int dim = alignedSeqs.size();
    Rcpp::NumericMatrix simMatrix(dim, dim);
    for (int i = 0; i < dim; ++i) {
        for (int j = i; j < dim; ++j) {
            if (i == j) { simMatrix(i, j) = 1; } else {
                simMatrix(j, i) = simMatrix(i, j) = compare(
                    Rcpp::as<std::string>(alignedSeqs[i]),
                    Rcpp::as<std::string>(alignedSeqs[j])
                );
            }
        }
    }
    return simMatrix;
}

// [[Rcpp::export]]
SEXP runTreemer(
        const Rcpp::ListOf<Rcpp::IntegerVector> &tipPaths, 
        const Rcpp::ListOf<Rcpp::CharacterVector> &alignedSeqs,
        Rcpp::NumericMatrix &simMatrixInput,
        const float similarity, const bool getTips
) {
    std::map<std::pair<int, int>, float> simMatrix;
    int nrow = simMatrixInput.nrow(), ncol = simMatrixInput.ncol();
    for (int i = 0; i < nrow; ++i) {
        for (int j = 0; j < ncol; ++j) {
            if (!R_IsNA(simMatrixInput(i, j))) {
                simMatrix[std::make_pair(i + 1, j + 1)] = simMatrixInput(i, j);
            }
        }
    }
    Treemer::BySimilarity match(tipPaths, alignedSeqs, similarity, simMatrix);
    if (getTips) {
        return Rcpp::wrap(match.getTips());
    } else {
        return Rcpp::wrap(match.getPaths());
    }
}

// [[Rcpp::export]]
Rcpp::IntegerVector divergentNode(
        const Rcpp::ListOf<Rcpp::IntegerVector> &paths
) {
    std::vector<int> res;
    for (int i = 0; i < paths.size() - 1; i++) {
        for (int j = i + 1; j < paths.size(); j++) {
            Rcpp::IntegerVector::const_iterator 
            q  = paths[i].begin(), s = paths[j].begin();
            do { q++, s++; } while (*q == *s);
            if (--q != paths[i].begin()) { res.push_back(*q); }
        }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::IntegerVector getReference(
        const std::string &refSeq, 
        const char gapChar
) {
    std::vector<int> res;
    for (unsigned int i = 0; i < refSeq.size(); i++) {
        if (refSeq[i] != gapChar) { res.push_back(i + 1); }
    }
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> ancestralPaths(
        const Rcpp::ListOf<Rcpp::IntegerVector> &paths,
        const int minLen
) {
    std::vector<Rcpp::IntegerVector> res;
    for (int i = 0; i < paths.size(); ++i) {
        if (paths[i].size() >= minLen) { 
            res.push_back(paths[i][Rcpp::Range(0, minLen - 1)]);
        }
    }
    return wrap(res);
}

// [[Rcpp::export]]
Rcpp::CharacterVector summarizeAA(
        const Rcpp::CharacterVector &seqs, 
        const int siteIndex, 
        const float tolerance
) {
    int nseq = seqs.size();
    std::map<char, int> aaSummary;
    for (int i = 0; i < nseq; ++i) { aaSummary[seqs[i][siteIndex]]++; }
    auto it = aaSummary.begin();
    int currentMax = it->second;
    char maxArg = it->first;
    for (++it; it != aaSummary.end(); ++it) {
        if (it->second > currentMax) {
            currentMax = it->second;
            maxArg = it->first;
        }
    }
    // The default tolerance would be a number
    int tt = static_cast<int>(tolerance);
    // Get the tolerance number if receiving float
    if (tolerance < 0.5) {
        tt = tolerance * nseq;
    } else if (tolerance < 1) {
        tt = nseq - tolerance * nseq;
    }
    // Test if the number of non-dominant AAs within tolerance
    nseq -= currentMax;
    if (nseq > tt) {
        return NA_STRING;
    } else {
        // The functin will return the fixed AA if the number
        // of non-dominant AA doesn't exceed the tolerance "tt"
        Rcpp::CharacterVector res = Rcpp::wrap(maxArg);
        res.attr("n") = Rcpp::IntegerVector::create(nseq);
        return res;
    }
}

// [[Rcpp::export]]
Rcpp::IntegerVector tableAA(
        const Rcpp::CharacterVector &seqs,
        const int siteIndex
) {
    std::map<std::string, int> aaSummary;
    for (int i = 0; i < seqs.size(); ++i) {
        aaSummary[std::string(1, seqs[i][siteIndex])]++;
    }
    return Rcpp::wrap(aaSummary);
}

// [[Rcpp::export]]
Rcpp::ListOf<Rcpp::IntegerVector> minimizeEntropy(
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries
) {
    const segIndex terminal = nodeSummaries.size();
    segment all;
    // Transform R list to a vector of AA and mapped frequency.
    // Get all the possible segment points
    for (segIndex i = 0; i < terminal; ++i) {
        all.push_back(i);
        Rcpp::IntegerVector summary = nodeSummaries[i].attr("aaSummary");
        Rcpp::CharacterVector aa = summary.names();
        aaSummary node;
        for (unsigned int j = 0; j < aa.size(); ++j) {
            node.emplace(aa[j], summary[j]);
        }
        Segmentor::aaSummaries.push_back(node);
    }
    // The segment point of "0" should be removed
    all.erase(all.begin());
    // The ultimate "parent" of all "Segmentor". Assume the "entropy"
    // of "parent" is minimum
    Segmentor *parent = new Segmentor(all, terminal);
    segment final = parent->m_used;
    float minEntropy = parent->m_entropy;
    // Initialize the search list and search depth
    std::vector<Segmentor *> segList;
    const unsigned int maxDepth = terminal / 3;
    unsigned int depth = 0;
    // Find the "minEntropy" among "Segmentor" in "segList" and make
    // the "Segmentor" the "parent" for next round. The new "minEntropy"
    // should be decreasing but increasing is allowed. So "final" is
    // returned when its "minEntropy" is not beaten for "maxDepth" of
    // consecutive times.
    while (true) {
        // Children "Segmentor" of the current "parent"
        for (unsigned int i = 0; i < parent->m_open.size(); ++i) {
            Segmentor *seg = new Segmentor(parent, i);
            segList.push_back(seg);
        }
        delete parent;
        // Re-calculate the best result from the search list and assign
        // to "tempMin". Temporarily assume the first "Segmentor" in the
        // list is the best result.
        auto it = segList.begin();
        // This "tempMin" is always one of the "Segmentors" in the search
        // list. And once it's found, it's going to be the new "parent"
        // and removed from the search list.
        auto rm = it;
        Segmentor *tempMin = *it;
        for (++it; it != segList.end(); ++it) {
            if ((*it)->m_entropy < tempMin->m_entropy) {
                tempMin = *it;
                rm = it;
            }
        }
        segList.erase(rm);
        // The search stops when none of the children "Segmentor" can
        // beat the "minEntropy" for "maxDepth" of consecutive times
        if (tempMin->m_entropy > minEntropy) {
            ++depth;
            if (tempMin->m_used.size() == terminal || depth == maxDepth) {
                for (const auto i: segList) {
                    delete i;
                }
                break;
            }
            // "final" and "minEntropy" stay unchanged if "tempMin"
            // cannot beat the previous "minEntropy".
        } else {
            // "tempMin" will be the new "final" and "minEntropy"
            // if the previous "minEntropy" is beaten by it.
            final = tempMin->m_used;
            minEntropy = tempMin->m_entropy;
            // And "depth" is reset to 0
            depth = 0;
        }
        // "tempMin" is going to be the new "parent" for getting new
        // children "Segmentor" no matter what. So "tempMin" doesn't
        // necessarily have to give the new "final" and "minentropy"
        parent = tempMin;
    }
    // Group tips by the segment indices in "final"
    std::vector<Rcpp::IntegerVector> res;
    segIndex start = 0;
    std::string prevFixedAA = "";
    aaSummary combNode;
    std::vector<int> combTips;
    for (const segIndex &end: final) {
        aaSummary node;
        std::vector<int> tips;
        for (unsigned int i = start; i < end; ++i) {
            Rcpp::IntegerVector nodeTips = nodeSummaries[i];
            tips.insert(tips.end(), nodeTips.begin(), nodeTips.end());
            Rcpp::IntegerVector summary = nodeTips.attr("aaSummary");
            Rcpp::CharacterVector aa = summary.names();
            for (unsigned int j = 0; j < aa.size(); ++j) {
                node[Rcpp::as<std::string>(aa.at(j))] += summary.at(j);
            }
        }
        // The most dominant AA in the current segment
        auto it = node.begin();
        std::string fixedAA = it->first;
        int maxFreq = it->second;
        for (++it; it != node.end(); ++it) {
            if (it->second > maxFreq) {
                fixedAA = it->first;
                maxFreq = it->second;
            }
        }
        // Combine with the previous segment if the most dominant
        // AA is the same as the previous. 
        if (fixedAA == prevFixedAA) {
            for (const auto &i: node) {
                combNode[i.first] += i.second;
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
        start = end;
        prevFixedAA = fixedAA;
    }
    Rcpp::IntegerVector combined = Rcpp::wrap(combTips);
    combined.attr("aaSummary") = Rcpp::wrap(combNode);
    combined.attr("AA") = prevFixedAA;
    res.push_back(combined);
    res.erase(res.begin());
    return Rcpp::wrap(res);
}

// [[Rcpp::export]]
Rcpp::CharacterVector tip2colorEdge(
        Rcpp::CharacterVector &colorEdge,
        const std::string &color,
        const Rcpp::IntegerMatrix &treeEdge,
        const Rcpp::IntegerVector &tips,
        const int rootNode
) {
    std::map< int, std::pair<int, int> > nodeLink;
    for (int i = 0; i < treeEdge.nrow(); ++i) {
        nodeLink[treeEdge(i, 1)].first = treeEdge(i, 0);
        nodeLink[treeEdge(i, 1)].second = i;
    }
    for (int i = 0; i < tips.size(); ++i) {
        int cn = tips[i], an;
        do {
            an = nodeLink[cn].first;
            colorEdge[nodeLink[cn].second] = color;
            cn = an;
        } while (an != rootNode);
    }
    return colorEdge;
}
