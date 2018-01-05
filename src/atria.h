#ifndef ATRIA_H
#define ATRIA_H

// This includes the code for the nearest neighbor searcher.
#define VERBOSE
#define PARTIAL_SEARCH
#include "NNSearcher/metric.h"
#include "NNSearcher/nearneigh_search.h"
#include "NNSearcher/point_set.h"
#undef PARTIAL_SEARCH
#undef VERBOSE

#include <Rcpp.h>

// Lots of boilderplate code below. C++ does not support virtual member
// templates, so we have to flesh out the dispatch logic for every exported
// function here for all metrics.
class Searcher {
private:
  std::string metric_;
  ATRIA<rm_point_set<euclidian_distance>> *euclidian_;
  ATRIA<rm_point_set<manhattan_distance>> *manhattan_;
  ATRIA<rm_point_set<maximum_distance>> *maximum_;
  ATRIA<rm_point_set<hamming_distance>> *hamming_;

public:
  Searcher() = delete;
  Searcher(const Searcher &) = delete;

  Searcher(const Rcpp::NumericMatrix x, const std::string metric,
           const long excl = 0, const long minpts = 64, const uint32 seed = 9345356234)
      : metric_("euclidian"), euclidian_(nullptr), manhattan_(nullptr),
        maximum_(nullptr), hamming_(nullptr) {
    // Sanitize input metric.
    if (metric.compare("euclidian") == 0) {
      metric_ = "euclidian";
    } else if (metric.compare("manhattan") == 0) {
      metric_ = "manhattan";
    } else if (metric.compare("maximum") == 0) {
      metric_ = "maximum";
    } else if (metric.compare("hamming") == 0) {
      metric_ = "hamming";
    } else {
      std::string exception_string = "Unknown metric " + metric + " specified.";
      throw Rcpp::exception(exception_string.c_str());
    }
    Rcpp::Rcout << "Using " << metric_ << " metric." << endl;
    if (metric_ == "euclidian") {
      rm_point_set<euclidian_distance> points(x);
      euclidian_ = new ATRIA<rm_point_set<euclidian_distance>>(
          std::move(points), excl, minpts, seed);
    } else if (metric_ == "manhattan") {
      rm_point_set<manhattan_distance> points(x);
      manhattan_ = new ATRIA<rm_point_set<manhattan_distance>>(
          std::move(points), excl, minpts, seed);
    } else if (metric_ == "maximum") {
      rm_point_set<maximum_distance> points(x);
      maximum_ = new ATRIA<rm_point_set<maximum_distance>>(std::move(points),
                                                           excl, minpts, seed);
    } else if (metric_ == "hamming") {
      rm_point_set<hamming_distance> points(x);
      hamming_ = new ATRIA<rm_point_set<hamming_distance>>(std::move(points),
                                                           excl, minpts, seed);
    }
  }
  ~Searcher() {
    delete euclidian_;
    delete maximum_;
    delete manhattan_;
    delete hamming_;
  }

  // Search for k nearest neighbors of the point query_point, excluding
  // points with indices between first and last from the search. Returns a
  // sorted vector of neighbors (by reference).
  template <class ForwardIterator>
  long search_k_neighbors(vector<neighbor> &v, const long k,
                          ForwardIterator query_point, const long first = -1,
                          const long last = -1, const double epsilon = 0) {

    if (metric_ == "euclidian") {
      return euclidian_->search_k_neighbors(v, k, query_point, first, last,
                                            epsilon);
    } else if (metric_ == "manhattan") {
      return manhattan_->search_k_neighbors(v, k, query_point, first, last,
                                            epsilon);
    } else if (metric_ == "maximum") {
      return maximum_->search_k_neighbors(v, k, query_point, first, last,
                                          epsilon);
    } else if (metric_ == "hamming") {
      return hamming_->search_k_neighbors(v, k, query_point, first, last,
                                          epsilon);
    }
    return 0;
  };

  // Count the number of points within distance 'radius' from the query point,
  // excluding points with indices between first and last from the search.
  template <class ForwardIterator>
  long count_range(const double radius, ForwardIterator query_point,
                   const long first = -1, const long last = -1) {
    if (metric_ == "euclidian") {
      return euclidian_->count_range(radius, query_point, first, last);
    } else if (metric_ == "manhattan") {
      return manhattan_->count_range(radius, query_point, first, last);
    } else if (metric_ == "maximum") {
      return maximum_->count_range(radius, query_point, first, last);
    } else if (metric_ == "hamming") {
      return hamming_->count_range(radius, query_point, first, last);
    }
    return 0;
  };

  // Search points within distance 'radius' from the query point,  excluding
  // points with indices between first and last  Returns an unsorted vector v of
  // neigbors by reference.
  template <class ForwardIterator>
  long search_range(vector<neighbor> &v, const double radius,
                    ForwardIterator query_point, const long first = -1,
                    const long last = -1) {
    if (metric_ == "euclidian") {
      return euclidian_->search_range(v, radius, query_point, first, last);
    } else if (metric_ == "manhattan") {
      return manhattan_->search_range(v, radius, query_point, first, last);
    } else if (metric_ == "maximum") {
      return maximum_->search_range(v, radius, query_point, first, last);
    } else if (metric_ == "hamming") {
      return hamming_->search_range(v, radius, query_point, first, last);
    }
    return 0;
  };

  // Returns an approximation of the data set radius such that any pairwise
  // distance in the data set is smaller than twice this radius. This bound is
  // not necessarily tight.
  double data_set_radius() const {
    if (metric_ == "euclidian") {
      return euclidian_->data_set_radius();
    } else if (metric_ == "manhattan") {
      return manhattan_->data_set_radius();
    } else if (metric_ == "maximum") {
      return maximum_->data_set_radius();
    } else if (metric_ == "hamming") {
      return hamming_->data_set_radius();
    }
    return 0.0;
  };

  long total_tree_nodes() const {
    if (metric_ == "euclidian") {
      return euclidian_->total_tree_nodes();
    } else if (metric_ == "manhattan") {
      return manhattan_->total_tree_nodes();
    } else if (metric_ == "maximum") {
      return maximum_->total_tree_nodes();
    } else if (metric_ == "hamming") {
      return hamming_->total_tree_nodes();
    }
    return 0;
  };

  long number_of_points() const {
    if (metric_ == "euclidian") {
      return euclidian_->number_of_points();
    } else if (metric_ == "manhattan") {
      return manhattan_->number_of_points();
    } else if (metric_ == "maximum") {
      return maximum_->number_of_points();
    } else if (metric_ == "hamming") {
      return hamming_->number_of_points();
    }
    return 0;
  };
};

#endif
