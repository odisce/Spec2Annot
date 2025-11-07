#include <Rcpp.h>
#include "mz_ppm.h"
using namespace Rcpp;

//' Get mass deviation in ppm
//'
//' This function returns a numeric value
//' corresponding to the mass deviation
//' between `massa` and `massb` in ppm.
//'
//' @param massa First mass
//' @param massb Second mass
//' @export
//'
// [[Rcpp::export]]
double mz_ppm(
  double massa = 120.1253,
  double massb = 120.1263
) {  
  Rcpp::NumericVector in_data = {massa, massb};
  double mass_error = max(in_data) - min(in_data);
  return (double)((mass_error / mean(in_data)) * 1e6);
}

//' Get mass range with ppm
//'
//' @param mass First mass
//' @param ppm ppm tolerance mass +- ppm(mass)
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericVector mz_range(
  double mass = 120.1253,
  double ppm = 10
) {  
  double mzdiff = (ppm * mass / 1e6);
  Rcpp::NumericVector output = {mass - mzdiff, mass + mzdiff};
  return output;
}
