// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>

using namespace Rcpp;

#include "atria.h"

//' Create ATRIA searcher (ball-tree).
//' @title FUNCTION_TITLE
//' @description FUNCTION_DESCRIPTION
//' @param x PARAM_DESCRIPTION
//' @param metric PARAM_DESCRIPTION, Default: 'euclidian'
//' @param exclude_samples PARAM_DESCRIPTION, Default: 0
//' @param cluster_max_points PARAM_DESCRIPTION, Default: 64
//' @param seed PARAM_DESCRIPTION, Default: 93453562
//' @return OUTPUT_DESCRIPTION
//' @details DETAILS
//' @examples
//' \dontrun{
//' if(interactive()){
//'  #EXAMPLE1
//'  }
//' }
//' @rdname create_searcher
//' @export
// [[Rcpp::export]]
XPtr<Searcher> create_searcher(NumericMatrix x,
                               const string metric = "euclidian",
                               const long exclude_samples = 0,
                               const long cluster_max_points = 64,
                               const uint32 seed=93453562L) {
  Searcher *s = new Searcher(x, metric, exclude_samples, cluster_max_points, seed);
  XPtr<Searcher> searcher(s);
  Rcout << "Approx. dataset radius: " << s->data_set_radius() << std::endl;
  return searcher;
}

//' Release searcher
//'
//' Destroy ATRIA searcher object and free all allocated memory.
//' @title FUNCTION_TITLE
//' @description FUNCTION_DESCRIPTION
//' @param searcher PARAM_DESCRIPTION
//' @return OUTPUT_DESCRIPTION
//' @details DETAILS
//' @examples
//' \dontrun{
//' if(interactive()){
//'  #EXAMPLE1
//'  }
//' }
//' @rdname release_searcher
//' @export
//[[Rcpp::export]]
bool release_searcher(XPtr<Searcher> searcher) {
  searcher.release();
  return !searcher;
}

//' Number of points,
//'
//' Return number of points used to create the searcher.
//'
//' @param searcher An exernal pointer to an ATRIA searcher.
//' @title FUNCTION_TITLE
//' @description FUNCTION_DESCRIPTION
//' @param searcher PARAM_DESCRIPTION
//' @return OUTPUT_DESCRIPTION
//' @details DETAILS
//' @examples
//' \dontrun{
//' if(interactive()){
//'  #EXAMPLE1
//'  }
//' }
//' @rdname number_of_points
//' @export
//[[Rcpp::export]]
long number_of_points(XPtr<Searcher> searcher) {
  return searcher->number_of_points();
}

//' Data set radius
//'
//' This function returns a logical vector identifying if
//' there are leading NA, marking the leadings NA as TRUE and
//' everything else as FALSE.
//'
//' @param x An integer vector
//' @export
//' @title FUNCTION_TITLE
//' @description FUNCTION_DESCRIPTION
//' @param searcher PARAM_DESCRIPTION
//' @return OUTPUT_DESCRIPTION
//' @details DETAILS
//' @examples
//' \dontrun{
//' if(interactive()){
//'  #EXAMPLE1
//'  }
//' }
//' @rdname data_set_radius
//' @export
//[[Rcpp::export]]
double data_set_radius(XPtr<Searcher> searcher) {
  return searcher->data_set_radius();
}


//' @title FUNCTION_TITLE
//' @description FUNCTION_DESCRIPTION
//' @param searcher PARAM_DESCRIPTION
//' @param k PARAM_DESCRIPTION
//' @param query_points PARAM_DESCRIPTION
//' @param exclude PARAM_DESCRIPTION, Default: matrix()
//' @param epsilon PARAM_DESCRIPTION, Default: 0
//' @return OUTPUT_DESCRIPTION
//' @details DETAILS
//' @examples
//' \dontrun{
//' if(interactive()){
//'  #EXAMPLE1
//'  }
//' }
//' @rdname search_k_neighbors
//' @export
//[[Rcpp::export]]
List search_k_neighbors(XPtr<Searcher> searcher, const long k,
                        NumericMatrix query_points,
                        IntegerMatrix exclude = IntegerMatrix(),
                        const double epsilon = 0) {
  if (k <= 0) {
    throw Rcpp::exception("Number of neighbors must be positive.");
  }
  bool use_exclude = false;
  if ((exclude.nrow() > 1) || (exclude.ncol() > 1)) {
    if ((exclude.nrow() != query_points.nrow()) || (exclude.ncol() != 2)) {
      std::string exception_string =
          "Wrong dimensions for input argument exclude, expected " +
          std::to_string(query_points.nrow()) + " by 2";
      throw Rcpp::exception(exception_string.c_str());
    }
    use_exclude = true;
  }
  IntegerMatrix index(query_points.nrow(), k);
  NumericMatrix dist(query_points.nrow(), k);

  for (long n = 0; n < query_points.nrow(); n++) {
    vector<neighbor> v;
    const auto query_point = query_points(n, Rcpp::_);
    // Search for neighbors.
    long first = -1;
    long last = -1;
    if (use_exclude) {
      // Convert exclude from one-based to zero-based indexing.
      // Include neighbor if index < first || index > last
      first = exclude(n, 0) - 1;
      last = exclude(n, 1) - 1;
    }
    searcher->search_k_neighbors(v, k, query_point.begin(), first, last,
                                 epsilon);
    for (long d = 0; d < k; d++) {
      index(n, d) = v[d].index() + 1; // Convert back to one-based indexing.
      dist(n, d) = v[d].dist();
    }
  }
  // Returns an IntegerMatrix and a NumericMatrix
  return List::create(Named("index") = index, Named("dist") = dist);
}

//' @title FUNCTION_TITLE
//' @description FUNCTION_DESCRIPTION
//' @param searcher PARAM_DESCRIPTION
//' @param radius PARAM_DESCRIPTION
//' @param query_points PARAM_DESCRIPTION
//' @param exclude PARAM_DESCRIPTION, Default: matrix()
//' @return OUTPUT_DESCRIPTION
//' @details DETAILS
//' @examples
//' \dontrun{
//' if(interactive()){
//'  #EXAMPLE1
//'  }
//' }
//' @rdname search_range
//' @export
//[[Rcpp::export]]
List search_range(XPtr<Searcher> searcher, const double radius,
                  NumericMatrix query_points,
                  IntegerMatrix exclude = IntegerMatrix()) {
  if (radius < 0) {
    throw Rcpp::exception("Radius can not be negative.");
  }
  bool use_exclude = false;
  if ((exclude.nrow() > 1) || (exclude.ncol() > 1)) {
    if ((exclude.nrow() != query_points.nrow()) || (exclude.ncol() != 2)) {
      std::string exception_string =
          "Wrong dimensions for input argument exclude, expected " +
          std::to_string(query_points.nrow()) + " by 2";
      throw Rcpp::exception(exception_string.c_str());
    }
    use_exclude = true;
  }

  // Returns an List of lists
  IntegerVector count(query_points.nrow());
  List nn(query_points.nrow());

  for (long n = 0; n < query_points.nrow(); n++) {
    vector<neighbor> v;
    const auto query_point = query_points(n, Rcpp::_);
    long first = -1;
    long last = -1;
    if (use_exclude) {
      //  Convert exclude from one-based to zero-based indexing.
      first = exclude(n, 0) - 1;
      last = exclude(n, 1) - 1;
    }
    // Search for neighbors.
    searcher->search_range(v, radius, query_point.begin(), first, last);
    count(n) = v.size();
    IntegerVector index(v.size());
    NumericVector dist(v.size());
    for (long d = 0; d < v.size(); d++) {
      index(d) = v[d].index() + 1; // Convert back to one-based indexing.
      dist(d) = v[d].dist();
    }
    nn(n) = List::create(Named("index") = index, Named("dist") = dist);
  }
  // Returns an IntegerMatrix and a NumericMatrix
  return List::create(Named("count") = count, Named("nn") = nn);
}
