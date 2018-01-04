// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include<RcppArmadillo.h>

using namespace Rcpp;

#include "NNSearcher/metric.h"
#include "NNSearcher/point_set.h"
#include "portable_intrinsics.h"

float L2SqrSIMD(const float* pVect1, const float* pVect2, size_t qty) {
  size_t qty4  = qty/4;
  size_t qty16 = qty/16;

  const float* pEnd1 = pVect1 + 16 * qty16;
  const float* pEnd2 = pVect1 + 4  * qty4;
  const float* pEnd3 = pVect1 + qty;

  __m128  diff, v1, v2;
  __m128  sum = _mm_set1_ps(0);

  while (pVect1 < pEnd1) {
    v1   = _mm_loadu_ps(pVect1); pVect1 += 4;
    v2   = _mm_loadu_ps(pVect2); pVect2 += 4;
    diff = _mm_sub_ps(v1, v2);
    sum  = _mm_add_ps(sum, _mm_mul_ps(diff, diff));

    v1   = _mm_loadu_ps(pVect1); pVect1 += 4;
    v2   = _mm_loadu_ps(pVect2); pVect2 += 4;
    diff = _mm_sub_ps(v1, v2);
    sum  = _mm_add_ps(sum, _mm_mul_ps(diff, diff));

    v1   = _mm_loadu_ps(pVect1); pVect1 += 4;
    v2   = _mm_loadu_ps(pVect2); pVect2 += 4;
    diff = _mm_sub_ps(v1, v2);
    sum  = _mm_add_ps(sum, _mm_mul_ps(diff, diff));

    v1   = _mm_loadu_ps(pVect1); pVect1 += 4;
    v2   = _mm_loadu_ps(pVect2); pVect2 += 4;
    diff = _mm_sub_ps(v1, v2);
    sum  = _mm_add_ps(sum, _mm_mul_ps(diff, diff));
  }

  while (pVect1 < pEnd2) {
    v1   = _mm_loadu_ps(pVect1); pVect1 += 4;
    v2   = _mm_loadu_ps(pVect2); pVect2 += 4;
    diff = _mm_sub_ps(v1, v2);
    sum  = _mm_add_ps(sum, _mm_mul_ps(diff, diff));
  }

  float PORTABLE_ALIGN16 TmpRes[4];

  _mm_store_ps(TmpRes, sum);
  float res= TmpRes[0] + TmpRes[1] + TmpRes[2] + TmpRes[3];

  while (pVect1 < pEnd3) {
    float diff = *pVect1++ - *pVect2++;
    res += diff * diff;
  }

  return res;
}

// [[Rcpp::export]]
NumericMatrix all_distances(NumericMatrix x) {
  rm_point_set<euclidian_distance> points(x);
  NumericMatrix out(x.nrow(), x.nrow());
  for (long i=0; i < x.nrow(); i++) {
    for (long j=0; j < x.nrow(); j++) {
      out(i,j) = points.distance(i, j);
    }
  }
  return out;
}

// [[Rcpp::export]]
NumericMatrix all_distances2(NumericMatrix x) {
  rm_point_set<euclidian_distance> points(x);
  NumericMatrix out(x.nrow(), x.nrow());
  for (long i=0; i < x.nrow(); i++) {
    for (long j=0; j < x.nrow(); j++) {
      out(i,j) = sqrt(L2SqrSIMD(points.point_begin(i),
                           points.point_begin(j),
                           points.dimension()));
    }
  }
  return out;
}
