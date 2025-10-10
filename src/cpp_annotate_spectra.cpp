#include <Rcpp.h>
using namespace Rcpp;

struct c_unique {
  int current;
  c_unique() {current=0;}
  int operator()() {return ++current;}
} UniqueNumber;

void fun_print (int i) {
  Rcpp::Rcout << i << ' ';
};

//' Find all combination
//'
//' This function returns edges of networked
//' ions
//'
//' @param n number of pairs to create
//' @param r number of element in each pairs
//' @export
//'
// [[Rcpp::export]]
int find_comb(int n, int r) {
    // int n = 100;
    // int r = 2;
    std::vector<int> myints(r);
    std::vector<int>::iterator first = myints.begin(), last = myints.end();
    std::generate(first, last, UniqueNumber);
    std::for_each(first, last, fun_print);
    Rcpp::Rcout << std::endl;
    while((*first) != n-r+1){
        std::vector<int>::iterator mt = last;
        while (*(--mt) == n-(last-mt)+1);
        (*mt)++;
        while (++mt != last) *mt = *(mt-1)+1;
        std::for_each(first, last, fun_print);
        Rcpp::Rcout << std::endl;
    }
    return 0;
}

//' Annotate a matrix and return network edges
//'
//' This function returns edges of networked
//' ions
//'
//' @param mass A numeric vector
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericVector annotate_spectra(
  Rcpp::NumericVector mass
) {
  // Iterate over every unique combination of two elements
  int n = mass.length()-1;
  int r = 2;
  double diff_i = 0;
  std::vector<int> myints(r);
  std::vector<int>::iterator first = myints.begin(), last = myints.end();
  std::generate(first, last, UniqueNumber);
  std::for_each(first, last, fun_print);
  while((*first) != n-r+1) {
      std::vector<int>::iterator mt = last;
      while (*(--mt) == n-(last-mt)+1);
      (*mt)++;
      while (++mt != last) *mt = *(mt-1)+1;
      diff_i = mass[*last] - mass[*first];
      Rcout << *last << ":" << *first << "=" << diff_i << "\n";
      //Rcpp::Rcout << *last << ":" << *first;
  }
  return mass;
}
