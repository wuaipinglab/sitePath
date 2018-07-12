#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <map>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector timesTwo() {
  int *a = nullptr;
  std::vector<int> z;
  std::map<int, int> x;
  int y = 1;
  x[1] = y;
  x[3] = 4;
  x[5] = 6;
  return wrap(x);
}

/*** R
x <- timesTwo()
*/
