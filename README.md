<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

The atriar package contains two sets of functions for computing fractal or intrinsic dimensions of data sets.

*Nearest neighbor* based methods allow the estimation of the [correlation sum](https://en.wikipedia.org/wiki/Correlation_integral). and related properties. The correlation dimension is one of the intrinsic dimensions of a data set and can be viewed as a measure of its complexity. It is expected that the performance of many machine learning algorithms depends on the intrinsic dimension of the input data.

This package contains an implementation of the [ATRIA nearest neighbor algorithm](https://www.researchgate.net/publication/12238515_Fast_nearest-neighbor_searching_for_nonlinear_signal_processing). which is a variant of the ball-tree algorithm family. It supports both exact and approximate k-nearest neighbor and range searches.

*Boxcounting* methods bin the points of a data set into a D-dimensional grid and compute various statistics based on the distribution of points into grid boxes. Counting the number of non-empty boxes at different grid resolutions allow the estimation of the capacity dimension of the data set. The underlying workhorse for the boxcounting methods is a [ternary search tree](https://en.wikipedia.org/wiki/Ternary_search_tree) implemented in C++ which is both time and space efficient.

Part of the code was released before within the OpenTSTOOL Matlab toolbox by the DPI Goettingen, Germany.

Installation
------------

The package is not (yet) on CRAN. You can also use the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("christianmerkwirth/atriar", ref = "develop")
```

Nearest Neighbors Searching
---------------------------

### k-Nearest Neighbor Queries

### Range queries

Boxcounting
-----------

Fast box counting for a data set of points X (row vectors of integers of dimension D). We assume that each row of the input data set addresses single box in a D-dimensional space. By means of a fast ternary search tree algorithm, we count the number of times each box has been visited in all subspaces from 1, 2, 3, ... up to D dimensions. A subpsace of dimension 2 e.g. is constructed from the two leftmost rows of the input matrix X. The algorithm returns the boxcounting, information and correlation measures for all prefix-subspaces 1, ..., D.

``` r
# Create a million data points in 10-d.
D <- 10
X <- floor(10 * matrix(runif(1e6), ncol=D))
mode(X) <- "integer"

# Let's just count non-empty bins in in D-dimensional space.
y.unique <- unique(X);
bc <- boxcount(X)

boxes.count <- 2^(-bc$boxd[D])
unique.count <- nrow(y.unique)

round(boxes.count) == unique.count

library(microbenchmark)
microbenchmark(boxcount(X), unique(X), times=10)
```

Dimension Estimation
--------------------

References
----------

[Nearest-neighbor based methods for nonlinear time-series analysis, Christian Merkwirth, 2001](http://hdl.handle.net/11858/00-1735-0000-0006-B40F-A)
