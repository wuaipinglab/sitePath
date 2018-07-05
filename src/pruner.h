#ifndef RUNER_H
#define RUNER_H

#include "util.h"

class Pruner: public TreeAlignmentMatch {
public:
  Pruner(
    ListOf<IntegerVector> tipPaths, 
    ListOf<CharacterVector> alignedSeqs
  );
  map< string, vector<int> > groupTips ();
};

#endif