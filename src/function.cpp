#include <Rcpp.h>
using namespace Rcpp;

//' Group m/Z from a vector based on tolerance
//'
//' This function returns a matrix with m/Z groups
//' based on a tolerance.
//'
//' @param xx A numeric and sorted vector of m/Z
//' @param tt Tolerance in absolute m/Z to group peaks
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericVector mz_vec_aggregate(Rcpp::NumericVector xx, double tt) {
  int xsize = xx.length();
  // Sort the vector in ascending order
  std::sort(xx.begin(), xx.end());
  // Iterate through vector and aggregate if i +- in tolerance
  // works on sorted vector only
  int group_id = 1;
  Rcpp::NumericVector grp (xsize);

  for (int i = 0; i < (xsize - 1); i++) {
    grp[i] = group_id;
    if (i < xsize) {
      if ((xx[i+1] - xx[i]) > tt) {
        group_id = group_id + 1;
      } else {
        grp[i+1] = group_id;
      }
    }
  }
  return grp;
}
