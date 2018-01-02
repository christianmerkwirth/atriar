#ifndef POINT_SET_H
#define POINT_SET_H

#include <Rcpp.h>
#include <algorithm>

// This file gives an example class for the implementation of a point_set which
// can be used by the nearest neighbor algorithm This particular implementation
// can be used for use with Rcpp, where a point set is given as an R numeric
// matrix. The coordinates of one point is given by one row of this matrix.
template <class METRIC> class point_set_base {
protected:
  const long N; // number of points
  point_set_base<METRIC>(const long n) : N(n){};
public:
  ~point_set_base<METRIC>(){};
  inline long size() const { return N; };
  typedef METRIC Metric;
};

// Define a row major point_set, i.e. data belonging to the same point
// are stored consecutively in memory. An object of this class copies(!)
// the data from the input Rcpp:NumericMatrix and keeps the copied data
// in memory allocated on the heap until the object is deleted. Though this
// may not be memory-efficient, we considered the hassle of having the user
// care for providing the point set over and over again at every call to the
// nearest neighbor searcher worse than storing the data. Points can be accessed
// by an index, ranging from zero to N-1. The implementation of the points is
// not important as long as the class provides the possibility to calculate
// distances between two points of the point set and the distance between a
// point from the data set and an externally given point. Class point_set is
// parametrized by the METRIC that is used top compute distances. METRIC must be
// a class having an operator(). For possible implementations of a METRIC, see
// file "metric.h" in this directory.
template <class METRIC>
class rm_point_set : public point_set_base<METRIC> { 
  
protected:
  const long D; // dimension
  float* matrix_ptr; // points are stored row-major in a C style array
  const METRIC Distance; // a function object that calculates distances
public:
  rm_point_set() = delete;
  rm_point_set(const rm_point_set& from) = delete;
  rm_point_set(const Rcpp::NumericMatrix& m)
    : point_set_base<METRIC>(m.nrow()), D(m.ncol()), matrix_ptr(new float[m.nrow() * m.ncol()]), Distance(){
      for (long n=0; n < point_set_base<METRIC>::N; n++) {
        const auto v = m(n, Rcpp::_);
        std::copy(v.begin(), v.end(), matrix_ptr + n*D);
      }
#ifdef DEBUG
      Rcpp::Rcout << "Point set constructor called." << std::endl;
#endif
    };
  // Move constructor here.
  rm_point_set(rm_point_set&& from)
    : point_set_base<METRIC>(from.N), D(from.D), Distance(){
      matrix_ptr = from.matrix_ptr;
      from.matrix_ptr = nullptr;
#ifdef DEBUG
      Rcpp::Rcout << "Point set move constructor called." << std::endl;
#endif
    };
  ~rm_point_set(){
#ifdef DEBUG
    if (matrix_ptr == nullptr) {
      Rcpp::Rcout << "Point set destructor on nullptr called." << std::endl;
    } else {
      Rcpp::Rcout << "Point set destructor called." << std::endl;
    }
#endif
    delete[] matrix_ptr;
  };
  inline long dimension() const { return D; };
  
  typedef const float* row_iterator; // pointer that iterates over the elements
  // of one point in the rm_point_set (points are row vectors)
  
  row_iterator point_begin(const long n) const { return matrix_ptr + n*D; }
  row_iterator point_end(const long n) const {
    return matrix_ptr + (n + 1) * D; // past-the-end
  }
  
  template <class ForwardIterator>
  inline double distance(const long index1, ForwardIterator vec2) const {
    return Distance(point_begin(index1), point_end(index1), vec2);
  }
  
#ifdef PARTIAL_SEARCH
  template <class ForwardIterator>
  inline double distance(const long index1, ForwardIterator vec2,
                         const double thresh) const {
    return Distance(point_begin(index1), point_end(index1), vec2, thresh);
  }
#endif
  
  inline double distance(const long index1, const long index2) const {
    return Distance(point_begin(index1), point_end(index1),
                    point_begin(index2));
  }
};

#endif
