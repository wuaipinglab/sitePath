#include <cmath>
#include <string>
#include <vector>
#include <Rcpp.h>

#include "minEntropy.h"

float MinEntropy::shannonEntropy(
        const aaSummary &values,
        const unsigned int tipNum
) {
    float res = 0.0;
    for (
            aaSummary::const_iterator it = values.begin();
            it != values.end(); ++it
    ) {
        // Probability of drawing distinct AA
        float p = it->second / static_cast<float>(tipNum);
        res -= p * std::log(p);
    }
    return res;
}

Rcpp::ListOf<Rcpp::IntegerVector> MinEntropy::updatedSegmentation(
        const Rcpp::ListOf<Rcpp::IntegerVector> &nodeSummaries,
        const segment &final
) {
    std::vector<Rcpp::IntegerVector> res;
    // Keep track of previous AA to see if the segment can be combined
    std::string prevFixedAA = "";
    // Summarize the AA if combined
    aaSummary combNode;
    // Group the tips if combined
    std::vector<int> combTips;
    // Iterate segment points and grouping tips
    segIndex start = 0;
    // To get the ancestral nodes of the raw groups
    Rcpp::CharacterVector ancestralNodes = nodeSummaries.names();
    // Index of the ancestral node of the group
    segIndex aNodeIndex = 0;
    // To find the segmentation index of current group, there is a need to look
    // back at the last segmentation index that separate two groups with
    // different dominant AA
    unsigned int indexShift = 1;
    // Ancestral node of the group
    std::string aNode = Rcpp::as<std::string>(ancestralNodes.at(0));
    for (
            segment::const_iterator final_itr = final.begin();
            final_itr != final.end(); ++final_itr
    ) {
        std::vector<int> tips; // To store and grow the group
        aaSummary node; // To summarize the amino acids of the group
        for (unsigned int i = start; i < *final_itr; ++i) {
            Rcpp::IntegerVector nodeTips = nodeSummaries[i];
            // Group the tips with a segment
            tips.insert(tips.end(), nodeTips.begin(), nodeTips.end());
            // Summarize the AA within a segment
            Rcpp::IntegerVector summary = nodeTips.attr("aaSummary");
            Rcpp::CharacterVector aa = summary.names();
            for (int j = 0; j < aa.size(); ++j) {
                node[Rcpp::as<std::string>(aa.at(j))] += summary.at(j);
            }
        }
        // Find the most dominant AA in the current segment
        aaSummary::iterator node_itr = node.begin();
        std::string fixedAA = node_itr->first;
        int maxFreq = node_itr->second;
        for (++node_itr; node_itr != node.end(); ++node_itr) {
            if (node_itr->second > maxFreq) {
                fixedAA = node_itr->first;
                maxFreq = node_itr->second;
            }
        }
        // If the most dominant AA is the same as the previous segment.
        if (fixedAA == prevFixedAA) {
            // Combine with the previous segment
            for (node_itr = node.begin(); node_itr != node.end(); ++node_itr) {
                combNode[node_itr->first] += node_itr->second;
            }
            combTips.insert(combTips.end(), tips.begin(), tips.end());
            indexShift++;
        } else {
            // Add the previous combination as a new segment
            Rcpp::IntegerVector combined = Rcpp::wrap(combTips);
            combined.attr("aaSummary") = Rcpp::wrap(combNode);
            combined.attr("AA") = prevFixedAA;
            combined.attr("node") = aNode;
            res.push_back(combined);
            // Initiate a new segment to be combined
            combTips = tips;
            combNode = node;
            // Find the segmentation index for the next group
            aNodeIndex = (final_itr == final.begin()) ? 0 : *(final_itr - indexShift);
            aNode = Rcpp::as<std::string>(ancestralNodes.at(aNodeIndex));
            indexShift = 1;
        }
        // Update the starting segment point for the next segment
        start = *final_itr;
        // Update the previous fixed AA
        prevFixedAA = fixedAA;
    }
    // Add the last segment to the return list
    Rcpp::IntegerVector combined = Rcpp::wrap(combTips);
    combined.attr("aaSummary") = Rcpp::wrap(combNode);
    combined.attr("AA") = prevFixedAA;
    combined.attr("node") = aNode;
    res.push_back(combined);
    // The first combination of segments is actually empty
    res.erase(res.begin());
    Rcpp::ListOf<Rcpp::IntegerVector> res2 = Rcpp::wrap(res);
    res2.attr("final") = Rcpp::wrap(final);
    res2.attr("nodeSummaries") = nodeSummaries.names();
    return res2;
}
