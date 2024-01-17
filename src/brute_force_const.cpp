#include <Rcpp.h>
#include <Rcpp.h>
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

//' Group m/Z from a vector based on tolerance
//'
//' This function returns a matrix with m/Z groups
//' based on a tolerance.
//'
//' @param mass Targeted mass
//' @param ppm Mass tolerance in ppm to restrict results
//' @param maxiter_vc_ Vector of element maximum limits (optional) (size = n)
//' @param mass_vc Vector of elemnt masses (size = n)
//' @param name_vc Vector of elemnt names (size = n)
//' @param debugl Integer (0: no message, 1: short, 2: verbose)
//' @param debugit Integer for max iter to do
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix brute_force_const(
  double mass = 120,
  double ppm = 5,
  Rcpp::NumericVector mass_vc = 0,
  Nullable<Rcpp::NumericVector> maxiter_vc_ = R_NilValue,
  Rcpp::CharacterVector name_vc = 0,
  int debugl = 0,
  int debugit = 0
) {  
  // INPUT
  Rcpp::NumericVector cur_cnt = ceil(mass / mass_vc);
  Rcpp::NumericVector max_vc = ceil(mass / mass_vc);

  if (maxiter_vc_.isNotNull()) {
    max_vc = maxiter_vc_;

    for (int i=0; i<max_vc.length(); i++) {
      cur_cnt[i] = max_vc[i];
    }
  }
  
  double mass_diff = (ppm * mass / 1e6);
  double mass_max = mass + mass_diff;
  double mass_min = mass - mass_diff;
  Rcpp::NumericVector cur_cnt_prev = rep(0.0, name_vc.length());
  
  // OUTPUT
  //    Store results in a vector with n element at each iteration
  //    then parse to a matrix by creating a new line each n element
  //    Store int and double values in two different vector
  std::vector<int> out_intv = {};
  std::vector<double> out_doublev = {};
  int out_doubleN = 3;
  int out_intN = name_vc.length();
  
  // PROGRAM
  if (debugl != 0) {
    Rcout << "Element names:\n" << name_vc << "\n";
    Rcout << "Mass-range target: " << mass_min << "-" << mass_max << "\n";
  }

  // Initialize starting iterators
  for (int i = 0; i < name_vc.length(); i++) {
    if (i == 0) {
      cur_cnt_prev[i] = cur_cnt[i] = cur_cnt[i] + 1; 
    } else {
      cur_cnt_prev[i] = cur_cnt[i] = 0; 
    }
  }

  if (debugl != 0) {
    Rcout << "Starting: " << cur_cnt << "\n";
  }
  int cnt_pos = 0;
  bool stopcond = true;
  double cur_mass = 0;
  int iter_val = 0;
  
  while (stopcond) {
    if (debugl == 1) {
      Rcout << "\riter: " << iter_val << " cnt: " << cur_cnt;
    }

    if (debugit > 0) {
      if (iter_val >= debugit) {
        stopcond = false;
        break;
      }
    }

    // A: Calculate current mass
    //     If in range append to result, if not: substract current
    //     iterator and recalculate the others until max
    
    //   A1: Calculate current mass
    iter_val = iter_val+1;
    cur_mass = sum(cur_cnt * mass_vc);

    //   Check: 
    bool keepL = false;
    if ((cur_mass >= mass_min) & (cur_mass <= mass_max)) {
      keepL = true;
      if (debugl == 1) {
        Rcout << " OK\n";
      }
      // Save to output
      //    Element composition
      for (int i=0; i<out_intN; i++) {
        out_intv.push_back(cur_cnt[i]);
      }
      // ppm deviation
      double cur_mass_ppm = mz_ppm(cur_mass, mass);
      // double values
      out_doublev.push_back(cur_mass);
      out_doublev.push_back(cur_mass_ppm);
      out_doublev.push_back(iter_val);
    }

    if (debugl == 2) {
      Rcout << "Cur Iter: " << cur_cnt << " mass: " << cur_mass << " " << keepL << "\n";
    }
    // }

    //    A2: New iterators (decreas the right most > 0)
    //        decrease the right most > 0 (but the last) with stop
    
    // Check next position
    cnt_pos = 0;
    for (int i=0; i<cur_cnt.length()-1; i++) {
      if (cur_cnt[i] > 0) {
        cnt_pos = i;
      }
    }

    // If all iterators are zero, stop
    if (sum(cur_cnt) <= 0) {
      stopcond = false;
      break;
    }

    // Decrease current iterator and check if < 0
    cur_cnt[cnt_pos] = cur_cnt[cnt_pos]-1;
    if (cur_cnt[cnt_pos] < 0) {
      break;
    }
    // Reset last count
    for (int i=cur_cnt.length()-1; i<cur_cnt.length(); i++) {
      cur_cnt[i] = 0;
    }
    
    // Calculate next iterators by dividing differences with
    // element mass
    for (int i=cnt_pos+1; i<cur_cnt.length(); i++) {
      double mass_diff = mass - sum(cur_cnt * mass_vc);
      int divider = (int)round(mass_diff / mass_vc[i]);
      if (divider <= 0) {
        cur_cnt[i] = 0;
      } else if (divider > max_vc[i]) {
        cur_cnt[i] = max_vc[i];
      } else {
        cur_cnt[i] = divider;
      }
    }
  }

  if (debugl != 0) {
    Rcout << "\n\r\n";
  }

  // Convert output to Numeric matrix
  Rcpp::NumericMatrix ac(out_intN, out_intv.size()/out_intN, out_intv.begin());
  Rcpp::NumericMatrix a = transpose(ac);
  Rcpp::NumericMatrix bc(out_doubleN, out_doublev.size()/out_doubleN, out_doublev.begin());
  Rcpp::NumericMatrix b = transpose(bc);

  // This part comes from: 
  // https://stackoverflow.com/a/31921185
  int acoln = a.ncol();
  int bcoln = b.ncol();
  Rcpp::NumericMatrix out = Rcpp::no_init(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      out(_, j) = b(_, j - acoln);
    }
  }
  // --- //

  Rcpp::CharacterVector out_names = Rcpp::no_init(acoln + bcoln);

  for (int i=0; i<out_names.length(); i++) {
    if (i < acoln) {
      out_names[i] = name_vc[i];
    } else {
      if (i == acoln) {
        out_names[i] = "mass";
      } else if (i == acoln + 1) {
        out_names[i] = "ppm";
      } else if (i == acoln + 2) {
        out_names[i] = "iteration";
      }
    }
  }
  
  colnames(out) = out_names;

  if (debugl != 0) {
    Rcout << "\nFound: " << out.nrow() << " match in " << iter_val << " iterations." << "\n";
  }
  
  return out;
}
