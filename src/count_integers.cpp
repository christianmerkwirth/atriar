// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector count_integers(IntegerVector bins, long max_bin = -1) {
  // Count positive integers.
  if (max_bin < 0) {
    max_bin = *std::max_element(bins.begin(), bins.end());
  }
  if (max_bin < 1) {
    throw Rcpp::exception("Values in bins must be positive.");
  }
  IntegerVector count(max_bin);
  
  for (auto bin : bins) {
    if (bin < 1) {
      throw Rcpp::exception("Values in bins must be positive.");
    }
    if (bin > max_bin) {
      throw Rcpp::exception("Values in bins must be smaller or equal to max_bin.");
    }
    count[bin-1] += 1;   // Rcpp is 0-based, R is 1-based.
  }
  return count;
}
