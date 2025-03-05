#include <Rcpp.h>
#include "mz_ppm.h"
using namespace Rcpp;

//' Get index in range
//'
//' @export
//'
// [[Rcpp::export]]
std::vector<int> get_range(std::vector<double> input, double valA, double valB) {
  std::vector<int> output(2, 0);
  std::vector<double>::iterator mzlow, mzhigh;
  
  mzlow = std::lower_bound(input.begin(), input.end(), valA);
  mzhigh = std::upper_bound(input.begin(), input.end(), valB);
  
  if (mzlow <= input.begin()) {
    output[0] = 0;
  } else if (mzlow >= input.end()) {
    output[0] = input.size()-1;
  } else {
    output[0] = (mzlow - input.begin());
  }
  if (mzhigh <= input.begin()) {
    output[1] = 0;
  } else if (mzhigh >= input.end()) {
    output[1] = input.size()-1;
  } else {
    output[1] = (mzhigh - input.begin()) - 1;
  }
  return output;
}

template <typename T>
void sort_indices(std::vector<T> &data, std::vector<size_t> &indices){
    std::sort(indices.begin(), indices.end(), [&data](size_t a, size_t b){ return data[a] < data[b]; });
}

std::vector<size_t> sorted_index(Rcpp::NumericVector& indt) {
  std::vector<double> indt_v = as<std::vector<double> >(indt);
  std::vector<size_t> indt_v_index(indt_v.size());
  for (size_t it = 0 ; it != indt_v_index.size() ; it++)
  indt_v_index[it] = it;
  sort_indices(indt_v, indt_v_index);
  return indt_v_index;
}

//' Match two table
//'
//' This function match two table, searching every
//' entries in B which are contained in A mz +- ppmtol
//' and A.rt +- rttol
//'
//' @param db_dt a data.frame with at least mz and rt values (ref table)
//' @param exp_dt a data.frame with at least mz and rt values (exp table)
//' @param ppmtol a numeric value for the ppm tolerance
//' @param rttol a numeric value for the rt tolerance
//' @export
//'
// [[Rcpp::export]]
Rcpp::DataFrame match_tables(DataFrame db_dt, DataFrame exp_dt, double ppmtol, double rttol, bool debugL = false) {
  Rcpp::NumericVector mzrange(2), rtrange(2);
  Rcpp::NumericVector bigMZ, smallMZ, bigRT, smallRT;

  // Set shortest table as reference
  Rcpp::DataFrame smallDT, bigDT;
  if (exp_dt.nrows() < db_dt.nrows()) {
    smallDT = exp_dt;
    bigDT = db_dt;
  } else {
    smallDT = db_dt;
    bigDT = exp_dt;
  }
  
  bigMZ = bigDT["mz"];
  smallMZ = smallDT["mz"];
  bigRT = bigDT["rt"];
  smallRT = smallDT["rt"];

  // Get sorted m/Z index for smallDT and bigDT
  std::vector<size_t> smallDTmzindex, bigDTmzindex;
  smallDTmzindex = sorted_index(smallMZ);
  bigDTmzindex = sorted_index(bigMZ);

  // Loop over smallDT and search corresponding entries in bigDT
  //  Opt: may use binary search to get the starting point
  std::vector<int> matched_Small_index, matched_Big_index;
  size_t prevstart = 0;
  double rtsel = 0.0;
  for (size_t small_it : smallDTmzindex) {
    mzrange = mz_range(smallMZ[small_it], ppmtol);
    rtsel = smallRT[small_it];
    rtrange = {rtsel - rttol, rtsel + rttol};
    if (debugL) {
      std::cout << "Getting ref: " << small_it;
      std::cout << " rt: " << rtrange[0] << "-" << rtrange[1];
      std::cout << " mz: " << mzrange[0] << "-" << mzrange[1] << std::endl; 
    }
    bool firstL = false;
    for (size_t big_it = prevstart ; big_it != bigDTmzindex.size() ; big_it++) {
      if (debugL) {
        std::cout << "   exp it: " << big_it;
        std::cout << " mz: " << bigMZ[bigDTmzindex[big_it]];
        std::cout << " rt: " << bigRT[bigDTmzindex[big_it]];
      }
      if (bigMZ[bigDTmzindex[big_it]] >= mzrange[0]) {
        if (!firstL) {
          prevstart = big_it;
          firstL = true;
          if (debugL) {
            std::cout << " new start: " << prevstart;
          }
        }
        if (debugL) {
          std::cout << " >= mzmin ";
        }
        if (bigMZ[bigDTmzindex[big_it]] > mzrange[1]) {
          if (debugL) {
            std::cout << " > mzmax ";
            std::cout << std::endl;
          }
          break;
        } else {
          // Check if in rtime
          if (debugL) {
            std::cout << " x < mzmax ";
          }
          if (
            bigRT[bigDTmzindex[big_it]] >= rtrange[0] &&
              bigRT[bigDTmzindex[big_it]] <= rtrange[1]
          ) {
            if (debugL) {
              std::cout << " in rtrange ";
            }
            matched_Small_index.push_back(small_it + 1); // +1 to convert to R indexing
            matched_Big_index.push_back(bigDTmzindex[big_it] + 1); // +1 to convert to R indexing
          }
        }
      }
      if (debugL) {
        std::cout << std::endl;
      }
    }
  }

  // Create a DataFrame to store the results
  std::vector<int> exp_dt_index, db_dt_index;
  if (exp_dt.nrows() < db_dt.nrows()) {
    exp_dt_index = matched_Small_index;
    db_dt_index = matched_Big_index;
  } else {
    db_dt_index = matched_Small_index;
    exp_dt_index = matched_Big_index;
  }
  DataFrame result = DataFrame::create(
    Named("dbid") = db_dt_index,
    Named("expid") = exp_dt_index
  );
  
  return result;
      


  // // Extract columns from table A
  // if (ppmtol) {
  //   mz_A = db_dt["mz"];
  //   mz_B = exp_dt["mz"];
  // } else {
  //   mz_A(db_dt.nrows());
  //   mz_B(exp_dt.nrows());
  // }
  // if (rttol) {
  //   rt_A = db_dt["rt"];
  //   rt_B = exp_dt["rt"];
  // } else {
  //   rt_A(db_dt.nrows());
  //   rt_B(exp_dt.nrows());
  // }
  
  // std::vector<int> matched_B_index; // To store which row in B the match corresponds to
  // std::vector<int> matched_A_index; // To store which row in A the match corresponds to
  // // Loop through each row in table B
  // for (int i = 0; i < exp_dt.nrows(); ++i) {
  //   // Loop through each row in table A
  //   if (ppmtol) {
  //     mzrange = mz_range(mz_B[i], ppmtol);
  //   } else {
  //     mzrange = {0, 0};
  //   }
  //   if (rttol) {
  //     rtrange = {rt_B[i] - rttol, rt_B[i] + rttol};
  //   } else {
  //     rtrange = {0, 0};
  //   }
  //   for (int j = 0; j < db_dt.nrows(); ++j) {
  //     // Check if mz and rt from A are within the ranges specified in B
  //     if (
  //         (mz_A[j] >= mzrange[0] && mz_A[j] <= mzrange[1]) &&
  //         (rt_A[j] >= rtrange[0] && rt_A[j] <= rtrange[1])
  //     ) {
  //       // Store the matching rows from A
  //       matched_A_index.push_back(j + 1); // +1 to convert to R indexing
  //       matched_B_index.push_back(i + 1); // +1 to convert to R indexing
  //     }
  //   }
  // }
  

}