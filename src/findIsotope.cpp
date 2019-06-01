#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector findIsotope(NumericMatrix spec, int numIsotope, double ppm, NumericVector prop) {
  const double C13 = 1.003355;

  NumericVector isIsotope(spec.nrow());

  for (int n = 0; n < numIsotope; n++) {
    for (int i = 0; i < spec.nrow(); i++) {
      if (isIsotope(i) == 1) {
        continue;
      }
      double mz = spec(i, 0);
      double intensity = spec(i, 1);
      double mzDiffThr = (ppm * mz) / 1e6;
      double mzMin = mz + C13 * (n + 1) - mzDiffThr;
      double mzMax = mz + C13 * (n + 1) + mzDiffThr;
      for (int j = i + 1; j < spec.nrow(); j++) {
        if (isIsotope(j) == 1) {
          continue;
        }
        double intRatio = spec(j, 1) / intensity;
        if (spec(j, 0) >= mzMin && spec(j, 0) <= mzMax && intRatio < prop(n)) {
          isIsotope(j) = 1;
        }
      }

    }
  }
  return isIsotope;
}
