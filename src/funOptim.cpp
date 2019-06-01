#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
double funOptim(NumericVector x, NumericVector M, NumericMatrix S) {
  int Sr = S.nrow();
  int Sc = S.ncol();
  double  res = 0;
  NumericVector m(Sr);
  for (int i=0; i<Sr; i++) {
    m(i) = 0;
  }
  for (int i=0; i<Sr; i++) {
    for (int j=0; j<Sc; j++) {
      m(i) = m(i) + x(j) * S(i, j);
    }
  }
  for (int i=0; i<Sr; i++) {
    res = res + pow(M(i) - m(i), 2);
  }
  return res;
}
