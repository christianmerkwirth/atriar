#include "Rcpp.h"

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix henon(long length = 2000,
                    NumericVector params = NumericVector::create(),
                    long transient = 5000) {
  if (length < 1) {
    stop("Requested number of samples must be positive.");
  }
  if (params.size() == 0) {
    params = NumericVector::create(-1.4, 0.3, 0.2, 0.12);
  }
  if (params.size() != 4) {
    stop("Need parameter vector of length 4 (a,b,x0,y0)");
  }
  
  NumericMatrix out(length, 2);
  
  const double a = params(0);
  const double b = params(1);
  
  double xn = params(2);
  double yn = params(3);
  
  for (long i = 0; i < length + transient; i++) {
    const double xn1 = 1 + a * xn * xn + b * yn;
    const double yn1 = xn;
    
    if (fabs(xn1) > 1000000) {
      Rcpp::Rcerr << "Values exceed threshold of 1000000, limiting values"
                  << std::endl;
      xn = 1000000;
    } else {
      xn = xn1;
    }
    yn = yn1;
    if (i >= transient) {
      out(i - transient, 0) = xn;
      out(i - transient, 1) = yn;
    }
  }
  return out;
}
