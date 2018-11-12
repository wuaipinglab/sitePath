#ifndef SITEPATH_UTIL_H
#define SITEPATH_UTIL_H

#include <string>
#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;

const float compare(const std::string &query, const std::string &subject);

class TipSeqLinker {
public:
  TipSeqLinker(CharacterVector sequence, IntegerVector tipPath);
  void proceed();
  const int nextClade() const;
  const int currentClade() const;
  const int getTip() const;
  const int getRoot() const;
  const int getSeqLen() const;
  IntegerVector getPath() const;
  std::string getSeq() const;
private:
  std::string seq;
  IntegerVector path;
  const int tipIndex;
  int cIndex;
};

#endif // SITEPATH_UTIL_H