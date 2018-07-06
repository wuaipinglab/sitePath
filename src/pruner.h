#ifndef RUNER_H
#define RUNER_H

#include "util.h"

class Pruner: public TreeAlignmentMatch {
public:
  Pruner(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs
  );
  std::map< std::string, std::vector<int> > groupTips ();
};

#endif