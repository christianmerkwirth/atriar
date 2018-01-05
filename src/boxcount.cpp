// Fast box counting for a data set of points X (row vectors of integers of
// dimension D). We assume that each row of the input data set addresses single
// box in a D-dimensional space. By means of a fast ternary search tree
// algorithm, we count the number of times each box has been visited in all
// subspace from 1, 2, 3, ... up to D dimensions. A subpsace of dimension 2 e.g.
// is constructed from the two leftmost rows of the input matrix X. The
// algorithm return the returns boxcounting, information and correlation
// dimension (D0, D1 and D2) (using log2) for all prefix-subspaces 1, ..., D.
//
// A fast algorithm based on ternary search trees to store nonempty boxes is
// used. Based on code written originally by Christian Merkwirth and Joerg Wichard at
// the DPI Goettingen 1998. Adapted to R by Christian Merkwirth 2017/2018.

#include "Rcpp.h"
#include "ternary_search_tree.h"

using namespace std;
using namespace Rcpp;

inline double max(double a, double b) { return ((a >= b) ? a : b); }
inline double min(double a, double b) { return ((a <= b) ? a : b); }

class dim_estimator {
private:
  const long N_;
  const long dim_;
  long total_points_;

  vector<long> boxes_;   // scaling of number of boxes for dimensions from 1 to dim
  vector<double> info_;  // scaling of information for dimensions from 1 to dim
  vector<double> corr_;  // scaling of correlation for dimensions from 1 to dim
public:
  dim_estimator(const long n, const long dim)
      : N_(n), dim_(dim), total_points_(0), boxes_(dim), info_(dim), corr_(dim) {}
  ~dim_estimator() {};

  //  Mass is the absolute frequency of points falling into that box at level
  // "level" of the tree.
  void operator()(const long mass, const int level) {
    const double p_i = (((double)mass) / (double)N_); // relative frequency

    total_points_ += mass;                    // Check we got all points.
    if (mass)
      boxes_[level]++;                          // Count boxes with at least one point
    info_[level] += (p_i * log2(p_i));        // entropy
    corr_[level] += p_i * p_i;                // correlation
  }

  // Mainly used to check if all points were inserted correctly into the tree.
  long get_total_points() const {
    return total_points_;
  }
  // Functions to read out results below.
  double boxd(const long d) const { return -log2(boxes_[d]); }
  double infod(const long d) const { return info_[d]; }
  double corrd(const long d) const { return log2(corr_[d]); }
};

// [[Rcpp::export]]
List boxcount(IntegerMatrix x, bool verbose = false) {
  const long N = x.nrow();
  const long dim = x.ncol();

  NumericMatrix result(3, dim);

  dim_estimator estimator(N, dim);
  ternary_search_tree<long> tree(dim);

  // Fill up tree.
  if (verbose) {
    Rcout << "Filling ternary search tree." << std::endl;
  }
  for (long index = 0; index < N; index++) {
    long key[dim];
    for (long d = 0; d < dim; d++) {
      key[d] = x(index, d);
    }
    if (tree.insert(key)) {
      throw Rcpp::exception("Ran out of memory.");
    }
  }
  // While Traversing the tree all relevant quantities are calculated.
  if (verbose) {
    Rcout << "Reading out ternary search tree." << std::endl;
  }
  tree.traverse(estimator);
  if (verbose) {
    Rcout << "Input data matrix is " << N << " by " << dim << std::endl;
    Rcout << "Total nodes allocated   : " << tree.total_nodes()
          << std::endl;
    Rcout << "Total memory allocated   : "
          << tree.total_nodes() * sizeof(Tnode<int>) << std::endl;
  }
  if (estimator.get_total_points() != dim * N) {
    throw Rcpp::exception("Internal error, tree did not contain all points.");
  }
  NumericVector boxd(dim, 0.0);
  NumericVector infod(dim, 0.0);
  NumericVector corrd(dim, 0.0);

  // Copy results to output.
  for (long d = 0; d < dim; d++) {
    boxd(d) = estimator.boxd(d);
    infod(d) = estimator.infod(d);
    corrd(d) = estimator.corrd(d);
  }
  return List::create(Named("boxd") = boxd,
                      Named("infod") = infod,
                      Named("corrd") = corrd);
}
