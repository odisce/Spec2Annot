#include <Rcpp.h>
using namespace Rcpp;
void combn_recursive(NumericVector m, long unsigned int r, int start, std::vector<double>& current, List& result) {
  if (current.size() == r) {
      result.push_back(NumericVector(current.begin(), current.end()));
      return;
  }
  for (int i = start; i < m.size(); ++i) {
      current.push_back(m[i]);
      combn_recursive(m, r, i + 1, current, result);
      current.pop_back();
  }
}

// [[Rcpp::export]]
NumericMatrix combn_cpp(NumericVector m, int r) {
  List result;
  std::vector<double> current;
  combn_recursive(m, r, 0, current, result);
  
  int n = result.size();
  double diffval = 0.0;
  NumericMatrix out(n, r+1);
  bool diffL = false;
  if (r >= 2) {
    diffL = true;
  }
  
  for (int i = 0; i < n; ++i) {
      NumericVector row = result[i];
      for (int j = 0; j < r; ++j) {
          out(i, j) = row[j];
      }
      if (diffL) {
      // Add difference
        diffval = std::abs(out(i,0) - out(i,1));
        out(i, r) = diffval;
      }
  }
  
  return out;
}