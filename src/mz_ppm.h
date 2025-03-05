#ifndef __MZPPM__
#define __MZPPM__

double mz_ppm(double massa, double massb);

#endif // __MZPPM__

#ifndef __MZRANGE__
#define __MZRANGE__

Rcpp::NumericVector mz_range(double mass, double ppm);

#endif // __MZRANGE__