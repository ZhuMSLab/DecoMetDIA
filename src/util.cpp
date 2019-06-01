#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
bool checkContinuousPtsAboveThr(NumericVector v, int iStart, int num, double thr, int nSkipMax) {
  int cnt = 0;
  bool res = false;
  int nSkip = 0;
  int nSkipPre = 0;
  for (int i = iStart; i < v.length(); i++) {
    if (v[i] > thr) {
      cnt++;
      nSkip = 0;
      nSkipPre = 0;
    } else {
      if (cnt > 0) {
        nSkipPre = nSkip;
        nSkip++;
      } else {
        nSkip = nSkipMax + 1;
      }
      if (nSkip > nSkipMax) {
        cnt = 0;
        nSkip = 0;
        nSkipPre = 0;
      } else {
        if (nSkipPre < nSkipMax) {
          cnt++;
        } else {
          cnt = cnt - nSkip + 1;
          nSkipPre = 0;
        }
      }
    }
    if (cnt >= num) {
      return(true);
    }
  }

  return(res);
}

// [[Rcpp::export]]
NumericVector getContinuousPtsAboveThrIdx(NumericVector v, int iStart, int num, double thr, int nSkipMax) {
  int cnt = 0;
  int stidx = 0;
  int enidx = 0;
  int nv = v.length();
  int nSkip = 0;
  int nSkipPre = 0;
  NumericVector res(nv);
  for (int i = iStart; i < nv; i++) {
    res[i] = false;
    if (v[i] > thr) {
      cnt++;
      nSkip = 0;
      nSkipPre = 0;
      if (cnt == 1) {
        stidx = i;
      } else {
        enidx = i;
      }
    } else {
      if (cnt > 0) {
        nSkipPre = nSkip;
        nSkip++;
      } else {
        nSkip = nSkipMax + 1;
      }
      if (nSkip > nSkipMax) {
        cnt = 0;
        nSkip = 0;
        nSkipPre = 0;
      } else {
        if (nSkipPre < nSkipMax) {
          cnt++;
        } else {
          cnt = cnt - nSkip + 1;
          nSkipPre = 0;
        }
      }
    }

    if ((cnt == 0 || i == nv - 1) && (enidx - stidx + 1) >= num) {
      for (int j = stidx; j <= enidx; j++) {
        res[j] = true;
      }
      stidx = 0;
      enidx = 0;
    }
  }
  return(res);
}
