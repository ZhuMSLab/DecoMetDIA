#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix mergeEIC(NumericMatrix x, NumericMatrix y) {
  int xr = x.nrow()-1;
  int yr = y.nrow()-1;
  int is, ie;
  if (x(0, 0) < y(0, 0)) {
    is = x(0, 0);
  } else {
    is = y(0, 0);
  }

  if (x(xr, 0) > y(yr, 0)) {
    ie = x(xr, 0);
  } else {
    ie = y(yr, 0);
  }
  NumericMatrix dm(ie-is+1, 3);
  for (int i=0; i<=xr; i++) {
    int idx = x(0, 0) - is + i;
    dm(idx ,0) = x(i, 0);
    dm(idx, 1) = x(i, 1);
  }
  for (int i=0; i<=yr; i++) {
    int idx = y(0, 0) - is + i;

    dm(idx ,0) = y(i, 0);
    dm(idx, 2) = y(i, 1);
  }
  return dm;
}
