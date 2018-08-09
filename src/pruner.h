#ifndef PRUNER_H
#define PRUNER_H

#include "util.h"

class Pruner: public TreeAlignmentMatch {
public:
  Pruner(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs,
    const float simThreshold
  );
  std::map< int, std::vector<int> > getTips();
  std::vector<IntegerVector> getPaths();
};

#endif