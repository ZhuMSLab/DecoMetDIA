#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

int findMaxIdx2(double d[], int len) {
  int iMin = 0;
  for (int i = 1; i < len; i++){
    if (d[iMin] < d[i]) {
      iMin = i;
    }
  }
  return iMin;
}

// [[Rcpp::export]]
List extractEIC(List specExp, NumericMatrix mzRange, NumericVector mz) {
  int lenMZ = mz.length();
  int lenSpec = specExp.size();
  Rcpp::List result(lenMZ);

  for (int iMZ = 0; iMZ < lenMZ; iMZ++) {
    NumericVector mzSpec(lenSpec), intSpec(lenSpec);

    for (int iSpec = 0; iSpec < lenSpec; iSpec++) {
      NumericMatrix spec = as <NumericMatrix> (specExp[iSpec]);
      double mz1 = mz(iMZ);
      double mzStart = mzRange(iMZ, 0);
      double mzEnd = mzRange(iMZ, 1);
      int idxInRange[lenSpec];
      int nr = spec.nrow();
      int iInRange = 0;

      for (int iRec = 0; iRec < nr; iRec++) {
        if (spec(iRec, 0) >= mzStart && spec(iRec, 0) <= mzEnd) {
          idxInRange[iInRange] = iRec;
          iInRange++;
        }
      }

      if (iInRange > 1) {
        double intInRange[iInRange];
        for (int i = 0; i < iInRange; i++) {
          intInRange[i] = spec(idxInRange[i], 1);
        }
        int iMax = findMaxIdx2(intInRange, iInRange);
        mzSpec(iSpec) = spec(idxInRange[iMax], 0);
        intSpec(iSpec) = spec(idxInRange[iMax], 1);
      } else if (iInRange == 1){
        mzSpec(iSpec) = spec(idxInRange[0], 0);
        intSpec(iSpec) = spec(idxInRange[0], 1);
      } else {
        mzSpec(iSpec) = mz1;
        intSpec(iSpec) = 0;
      }
    }
    result(iMZ) = Rcpp::List::create(Rcpp::Named("mz") = mzSpec,
           Rcpp::Named("intensity") = intSpec);
  }
  return result;
}
