#ifndef RUNER_H
#define RUNER_H

#include "util.h"

class Pruner: public TreeAlignmentMatch {
public:
  Pruner(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs
  );
  std::map< int, std::vector<int> > groupTips ();
};

#endif