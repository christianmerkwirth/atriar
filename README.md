<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

The atriar package contains two sets of functions for computing fractal or intrinsic dimensions of data sets.

*Nearest neighbor* based methods allow the estimation of the [correlation sum](https://en.wikipedia.org/wiki/Correlation_integral). and related properties. The correlation dimension is one of the intrinsic dimensions of a data set and can be viewed as a measure of its complexity. It is expected that the performance of many machine learning algorithms depends on the intrinsic dimension of the input data.

This package contains an implementation of the [ATRIA nearest neighbor algorithm](https://www.researchgate.net/publication/12238515_Fast_nearest-neighbor_searching_for_nonlinear_signal_processing). which is a variant of the ball-tree algorithm family. It supports both exact and approximate k-nearest neighbor and range searches. During preprocessing, a search tree is constructed by dividing the set of points in two (sub)clusters. Each cluster is than subdivided until a minimum number of points is reached When performing a query, the triangle inequality is used to exclude cluster from further searching. ATRIA might be a good choice for unevenly distributed points in very high dimensional spaces.

*Boxcounting* methods bin the points of a data set into a D-dimensional grid and compute various statistics based on the distribution of points into grid boxes. Counting the number of non-empty boxes at different grid resolutions allow the estimation of the capacity dimension of the data set. The underlying workhorse for the boxcounting methods is a [ternary search tree](https://en.wikipedia.org/wiki/Ternary_search_tree) implemented in C++ which is both time and space efficient.

Part of the code was released before within the OpenTSTOOL Matlab toolbox by the DPI Goettingen, Germany.

Installation
------------

The package is not (yet) on CRAN. You can also use the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("christianmerkwirth/atriar")
```

The package is tested to build successfully on OS X and on Linux.

Nearest Neighbors Searching
---------------------------

Fast nearest neighbor searching usually starts from preprocessing a data set of points into an index structure (e.g. ball tree, random projections) that is then used to accelerate subsequent neighbor queries:

``` r
# Create a set of points in 4-d.
D <- 4
points <- matrix(runif(1e6), ncol=D)
# Creating the ATRIA nearest neighbor searcher object.
searcher = create_searcher(points, metric="euclidian")
```

Function create\_searcher does the preprocessing for a given set of points. The returned searcher object contains preprocessing information that is a mandatory input argument for k-NN and range queries. Preprocessing and querying is divided into separate functions to give the user the possibility to re-use the searcher object when doing multiple searches on the same point set. Note that as soon as the underlying point set is changed or modified, one has to recompute the searcher for the updated set of points.

``` r
# Deleting the searcher object, freeing the allocated memory.
release_searcher(searcher)
```

### k-Nearest Neighbor Queries

``` r
library(atriar)

# Create a set of points in 4-d.
D <- 4
points <- matrix(runif(1e6), ncol=D)
# Creating the ATRIA nearest neighbor searcher object.
searcher = create_searcher(points, metric="euclidian")
```

    ## Using euclidian metric.
    ## ATRIA Constructor
    ## Size of point set : 250000  points of dimension 4
    ## Number of points used : 250000
    ## MINPOINTS : 64
    ## Root center : 225456
    ## Root starting index  : 1
    ## Root length : 249999
    ## Root Rmax : 1.42061
    ## Created tree structure for ATRIA searcher
    ## Approx. dataset radius: 1.42061

``` r
k.max <- 8
rand.sample <- sample.int(nrow(points), size = 1000)
nn <- search_k_neighbors(
    searcher = searcher,
    k = k.max,
    query_points = points[rand.sample, ],
    exclude = cbind(rand.sample, rand.sample)
  )
  
str(nn)
```

    ## List of 2
    ##  $ index: int [1:1000, 1:8] 97552 109147 109538 230959 242810 110531 227740 79973 179699 160963 ...
    ##  $ dist : num [1:1000, 1:8] 0.0356 0.0285 0.0146 0.0301 0.0251 ...

``` r
# Cleanup, delete the searcher object.
release_searcher(searcher)
```

    ## ATRIA Destructor
    ## Total_clusters : 11525
    ## Total number of points in terminal nodes : 238475
    ## Average number of points in a terminal node : 41.3804
    ## Average percentage of points searched 0.450853% (1128)
    ## Average number of terminal nodes visited : 34.476

    ## [1] TRUE

The output list returned by search\_k\_neighbors contains fields *index* and *dist* which both are matrices with as many rows as there were query points. Indices are 1-based and can be used to index R matrices without rebasing.

In general, reference points query.points can be arbitrarily located, but it is also possible that query.points are taken from the preprocessed points set. To allow the user to avoid self-matches, search\_k\_neighbors accepts the optional input argument *exclude* which is a nrow(query.points) by 2 integer matrix. *exclude* specifies a range of indices (i.e. first, last) that is not eligible as nearest neighbors for the given query point.

### Approximate k-Nearest Neighbor Queries

Approximate nearest neighbors algorithms report neighbors to the query point q with distances possibly greater than the true nearest neighbors distances. The maximal allowed relative error, named *epsilon*, is given as a parameter to the algorithm. For epsilon =0, the approximate search returns the true (exact) nearest neighbor(s). Computing exact nearest neighbors for data set with intrinsic dimension much higher than 6 seems to be a very time-consuming task. Few algorithms seem to perform significantly better than a brute-force computation of all distances. However, it has been shown that by computing nearest neighbors approximately, it is possible to achieve significantly faster execution times with relatively small actual errors in the reported distances. See (<https://github.com/erikbern/ann-benchmarks>) for a more thorough comparison of k-NN implementations.

``` r
# Create a set of points in 12-d.
D <- 12
points <- matrix(runif(1e6), ncol=D)
```

    ## Warning in matrix(runif(1e+06), ncol = D): data length [1000000] is not a
    ## sub-multiple or multiple of the number of rows [83334]

``` r
# Creating the ATRIA nearest neighbor searcher object.
searcher = create_searcher(points, metric="euclidian")
```

    ## Using euclidian metric.
    ## ATRIA Constructor
    ## Size of point set : 83334  points of dimension 12
    ## Number of points used : 83334
    ## MINPOINTS : 64
    ## Root center : 75152
    ## Root starting index  : 1
    ## Root length : 83333
    ## Root Rmax : 2.17909
    ## Created tree structure for ATRIA searcher
    ## Approx. dataset radius: 2.17909

``` r
k.max <- 8
rand.sample <- sample.int(nrow(points), size = 1000)

library(microbenchmark)
microbenchmark(
  search_k_neighbors(searcher, k.max, points[rand.sample, ]),
  search_k_neighbors(searcher, k.max, points[rand.sample, ], epsilon=3.0),
  times=10)
```

    ## Unit: milliseconds
    ##                                                                     expr
    ##               search_k_neighbors(searcher, k.max, points[rand.sample, ])
    ##  search_k_neighbors(searcher, k.max, points[rand.sample, ], epsilon = 3)
    ##      min      lq    mean  median      uq    max neval
    ##  1797.29 1808.94 1862.28 1846.89 1898.60 1998.0    10
    ##   608.54  628.76  636.64  633.13  658.89  660.7    10

``` r
# Cleanup, delete the searcher object.
release_searcher(searcher)
```

    ## ATRIA Destructor
    ## Total_clusters : 3873
    ## Total number of points in terminal nodes : 79461
    ## Average number of points in a terminal node : 41.0227
    ## Average percentage of points searched 31.524% (26271)
    ## Average number of terminal nodes visited : 955.31

    ## [1] TRUE

### Range queries

In the task of range searching , we ask for all points of data set P that have distance r or less from the query point q. Sometimes range searching is called a fixed size approach, while k nearest neighbors searching is called a fixed mass approach.

``` r
# Create a set of points in 4-d.
D <- 4
points <- matrix(runif(1e6), ncol=D)
# Creating the ATRIA nearest neighbor searcher object.
searcher = create_searcher(points, metric="euclidian")
```

    ## Using euclidian metric.
    ## ATRIA Constructor
    ## Size of point set : 250000  points of dimension 4
    ## Number of points used : 250000
    ## MINPOINTS : 64
    ## Root center : 225456
    ## Root starting index  : 1
    ## Root length : 249999
    ## Root Rmax : 1.66016
    ## Created tree structure for ATRIA searcher
    ## Approx. dataset radius: 1.66016

``` r
radius <- 0.2
rand.sample <- sample.int(nrow(points), size = 4)
nn <- search_range(
  searcher = searcher,
  radius = radius,
  query_points = points[rand.sample, ],
  # Ignore samples with index smaller than the query point
  # so we avoid counting the same pairwise distance twice.
  exclude = cbind(rep(-1, 4), rand.sample)
)
  
str(nn)
```

    ## List of 2
    ##  $ count: int [1:4] 1003 217 53 413
    ##  $ nn   :List of 4
    ##   ..$ :List of 2
    ##   .. ..$ index: int [1:1003] 121651 211106 169694 102288 190528 203998 137243 144967 173606 170611 ...
    ##   .. ..$ dist : num [1:1003] 0.177 0.165 0.181 0.183 0.161 ...
    ##   ..$ :List of 2
    ##   .. ..$ index: int [1:217] 197855 200214 208185 205750 187958 217101 202395 229350 235010 195533 ...
    ##   .. ..$ dist : num [1:217] 0.178 0.179 0.123 0.18 0.19 ...
    ##   ..$ :List of 2
    ##   .. ..$ index: int [1:53] 230722 247656 247223 245359 231327 237901 241614 228430 240998 231456 ...
    ##   .. ..$ dist : num [1:53] 0.132 0.197 0.163 0.188 0.161 ...
    ##   ..$ :List of 2
    ##   .. ..$ index: int [1:413] 246711 228607 202026 249380 173811 172821 249298 192497 182745 192394 ...
    ##   .. ..$ dist : num [1:413] 0.192 0.197 0.199 0.199 0.178 ...

``` r
# Cleanup, delete the searcher object.
release_searcher(searcher)
```

    ## ATRIA Destructor
    ## Total_clusters : 11485
    ## Total number of points in terminal nodes : 238515
    ## Average number of points in a terminal node : 41.5314
    ## Average percentage of points searched 1.1931% (2983)
    ## Average number of terminal nodes visited : 172.5

    ## [1] TRUE

The index and distance vectors for a single query point have the length that is given in count. Both vectors are not sorted by distance.

Boxcounting
-----------

Fast box counting for a data set of points X (row vectors of integers of dimension D). We assume that each row of the input data set addresses single box in a D-dimensional space. By means of a fast ternary search tree algorithm, we count the number of times each box has been visited in all subspaces from 1, 2, 3, ... up to D dimensions. A subpsace of dimension 2 e.g. is constructed from the two leftmost rows of the input matrix X. The algorithm returns the boxcounting, information and correlation measures for all prefix-subspaces 1, ..., D.

``` r
# Create a million data points in 10-d.
D <- 10
X <- floor(10 * matrix(runif(1e6), ncol=D))
mode(X) <- "integer"

# Let's just count non-empty bins in a D-dimensional space.
y.unique <- unique(X);
bc <- boxcount(X)

boxes.count <- 2^(-bc$boxd[D])
unique.count <- nrow(y.unique)

round(boxes.count) == unique.count
```

    ## [1] TRUE

``` r
library(microbenchmark)
microbenchmark(boxcount(X), unique(X), times=10)
```

    ## Unit: milliseconds
    ##         expr     min      lq    mean  median      uq     max neval
    ##  boxcount(X)  46.909  48.542  52.298  51.771  53.461  62.634    10
    ##    unique(X) 564.321 616.289 646.607 630.731 691.876 760.775    10

Dimension Estimation
--------------------

Below a typical workflow for estimating the correlation dimension by analying the slope of the correlation sum versus the distance in log scale:

``` r
# Create a downsampled data set of the terated Henon map.
data <- henon(2e6, params = c(-1.4, 0.3, 0.1 * runif(2)))
data <- data[sample.int(n = nrow(data), size = 2e5), ]

# Here we create the ATRIA nearest neighbor searcher object.
searcher = create_searcher(data, metric="euclidian", cluster_max_points = 64)
```

    ## Using euclidian metric.
    ## ATRIA Constructor
    ## Size of point set : 200000  points of dimension 2
    ## Number of points used : 200000
    ## MINPOINTS : 64
    ## Root center : 180365
    ## Root starting index  : 1
    ## Root length : 199999
    ## Root Rmax : 2.38589
    ## Created tree structure for ATRIA searcher
    ## Approx. dataset radius: 2.38589

``` r
# We need to get an idea about typical small and large distances
# in the data set in order to generate distance bins.
dist.limits <- distlimits(searcher, data)
dist.breaks <- logspace(dist.limits[1], dist.limits[2], 32)

# Invoke the correlation sum computation.
res <- corrsum(searcher=searcher,
               data=data,
               dist.breaks = dist.breaks,
               min.actual.pairs = 2000,
               min.nr.samples.at.scale = 256,
               max.nr.samples.at.scale = 1024,
               batch.size = 128)
```

    ## Radius: 0.287709145408075 count: 1080356 1080356
    ## Radius: 0.287709145408075 count: 1268003 1268003
    ## Radius: 0.214861366891889 count: 864592 864592
    ## Radius: 0.214861366891889 count: 907664 907664
    ## Radius: 0.160458601054102 count: 589487 589487
    ## Radius: 0.160458601054102 count: 625632 625632
    ## Radius: 0.119830582038485 count: 454115 454115
    ## Radius: 0.119830582038485 count: 453189 453189
    ## Radius: 0.089489552428794 count: 293961 293961
    ## Radius: 0.089489552428794 count: 355114 355114
    ## Radius: 0.0668308528396692 count: 239696 239696
    ## Radius: 0.0668308528396692 count: 253811 253811
    ## Radius: 0.0499093220388085 count: 155397 155397
    ## Radius: 0.0499093220388085 count: 155779 155779
    ## Radius: 0.0372723124205731 count: 113577 113577
    ## Radius: 0.0372723124205731 count: 108775 108775
    ## Radius: 0.0278349858588857 count: 77965 77965
    ## Radius: 0.0278349858588857 count: 69564 69564
    ## Radius: 0.020787184573413 count: 52247 52247
    ## Radius: 0.020787184573413 count: 52447 52447
    ## Radius: 0.0155238822351045 count: 34883 34883
    ## Radius: 0.0155238822351045 count: 31614 31614
    ## Radius: 0.0115932448089975 count: 27154 27154
    ## Radius: 0.0115932448089975 count: 31416 31416
    ## Radius: 0.00865784235965265 count: 18410 18410
    ## Radius: 0.00865784235965265 count: 17378 17378
    ## Radius: 0.00646568200357685 count: 12237 12237
    ## Radius: 0.00646568200357685 count: 13580 13580
    ## Radius: 0.00482857530026162 count: 9982 9982
    ## Radius: 0.00482857530026162 count: 9491 9491
    ## Radius: 0.00360598300649468 count: 6986 6986
    ## Radius: 0.00360598300649468 count: 6229 6229
    ## Radius: 0.0026929503289353 count: 5018 5018
    ## Radius: 0.0026929503289353 count: 4463 4463
    ## Radius: 0.00201109696331106 count: 2783 2783
    ## Radius: 0.00201109696331106 count: 3922 3922
    ## Radius: 0.00150188845014384 count: 1996 1996
    ## Radius: 0.00150188845014384 count: 3190 3190
    ## Radius: 0.00112161121906411 count: 1439 1439
    ## Radius: 0.00112161121906411 count: 1442 1442
    ## Radius: 0.000837619948811782 count: 1273 1273
    ## Radius: 0.000837619948811782 count: 1165 1165
    ## Radius: 0.000625535093374763 count: 727 727
    ## Radius: 0.000625535093374763 count: 697 697
    ## Radius: 0.000467149992784257 count: 480 480
    ## Radius: 0.000467149992784257 count: 537 537
    ## Radius: 0.000348867902168342 count: 324 324
    ## Radius: 0.000348867902168342 count: 325 325
    ## Radius: 0.000260534764087106 count: 271 271
    ## Radius: 0.000260534764087106 count: 214 214
    ## Radius: 0.000194567522194031 count: 228 228
    ## Radius: 0.000194567522194031 count: 166 166
    ## Radius: 0.000145303145341741 count: 120 120
    ## Radius: 0.000145303145341741 count: 127 127
    ## Radius: 0.000108512478383459 count: 96 96
    ## Radius: 0.000108512478383459 count: 89 89
    ## Radius: 0.000108512478383459 count: 86 86
    ## Radius: 0.000108512478383459 count: 102 102
    ## Radius: 0.000108512478383459 count: 93 93
    ## Radius: 0.000108512478383459 count: 113 113
    ## Radius: 0.000108512478383459 count: 119 119
    ## Radius: 0.000108512478383459 count: 111 111
    ## Radius: 8.10371856522924e-05 count: 59 59
    ## Radius: 8.10371856522924e-05 count: 73 73
    ## Radius: 8.10371856522924e-05 count: 65 65
    ## Radius: 8.10371856522924e-05 count: 55 55
    ## Radius: 8.10371856522924e-05 count: 62 62
    ## Radius: 8.10371856522924e-05 count: 79 79
    ## Radius: 8.10371856522924e-05 count: 64 64
    ## Radius: 8.10371856522924e-05 count: 64 64
    ## Radius: 6.05186201280713e-05 count: 31 31
    ## Radius: 6.05186201280713e-05 count: 34 34
    ## Radius: 6.05186201280713e-05 count: 35 35
    ## Radius: 6.05186201280713e-05 count: 41 41
    ## Radius: 6.05186201280713e-05 count: 36 36
    ## Radius: 6.05186201280713e-05 count: 51 51
    ## Radius: 6.05186201280713e-05 count: 34 34
    ## Radius: 6.05186201280713e-05 count: 46 46
    ## Radius: 4.51953427642534e-05 count: 38 38
    ## Radius: 4.51953427642534e-05 count: 42 42
    ## Radius: 4.51953427642534e-05 count: 27 27
    ## Radius: 4.51953427642534e-05 count: 25 25
    ## Radius: 4.51953427642534e-05 count: 28 28
    ## Radius: 4.51953427642534e-05 count: 30 30
    ## Radius: 4.51953427642534e-05 count: 29 29
    ## Radius: 4.51953427642534e-05 count: 31 31
    ## Radius: 3.37519097966164e-05 count: 33 33
    ## Radius: 3.37519097966164e-05 count: 29 29
    ## Radius: 3.37519097966164e-05 count: 23 23
    ## Radius: 3.37519097966164e-05 count: 29 29

``` r
# Print and visualize results.
print(res)
```

    ## $dists
    ##  [1] 0.0000e+00 3.3752e-05 4.5195e-05 6.0519e-05 8.1037e-05 1.0851e-04
    ##  [7] 1.4530e-04 1.9457e-04 2.6053e-04 3.4887e-04 4.6715e-04 6.2554e-04
    ## [13] 8.3762e-04 1.1216e-03 1.5019e-03 2.0111e-03 2.6930e-03 3.6060e-03
    ## [19] 4.8286e-03 6.4657e-03 8.6578e-03 1.1593e-02 1.5524e-02 2.0787e-02
    ## [25] 2.7835e-02 3.7272e-02 4.9909e-02 6.6831e-02 8.9490e-02 1.1983e-01
    ## [31] 1.6046e-01 2.1486e-01 2.8771e-01
    ## 
    ## $correlation.sum
    ##  [1] 0.0000e+00 1.7401e-06 2.4720e-06 3.4731e-06 4.9460e-06 7.1875e-06
    ##  [7] 1.0081e-05 1.4630e-05 2.0964e-05 2.9714e-05 4.2766e-05 6.2158e-05
    ## [13] 8.9825e-05 1.2805e-04 1.8189e-04 2.5833e-04 3.6571e-04 5.2019e-04
    ## [19] 7.3129e-04 1.0267e-03 1.4305e-03 1.9917e-03 2.7820e-03 4.0057e-03
    ## [25] 5.7792e-03 8.3173e-03 1.2012e-02 1.7060e-02 2.4390e-02 3.4854e-02
    ## [31] 4.8469e-02 6.6869e-02 9.4461e-02
    ## 
    ## $actual.pairs.count
    ##  [1]       0    2009     808    1004    1325    1785    2002    3032
    ##  [9]    4062    5388    7698   10963   14914   19640   26341   35425
    ## [17]   47034   63763   81656  106722  136233  174881  225219  318743
    ## [25]  414682  528493  672655  785556  945740 1075913 1048946  950782
    ## [33]  710226
    ## 
    ## $potential.pairs.count
    ##  [1] 1154513229 1154513229 1104088917 1002863257  899584342  796331606
    ##  [7]  691875883  666545028  641337805  615733287  589789061  565336142
    ## [13]  539050284  513781917  489296367  463413741  438018566  412749574
    ## [19]  386817963  361319264  337299672  311635099  284981120  260480062
    ## [25]  233818000  208221584  182061601  155628592  129022408  102821665
    ## [31]   77038870   51674766   25739778
    ## 
    ## $samples.used
    ##  [1]    0  512 1024 1024 1024 1024  256  256  256  256  256  256  256  256
    ## [15]  256  256  256  256  256  256  256  256  256  256  256  256  256  256
    ## [29]  256  256  256  256  256

``` r
x <- log2(res$dists)
y <-log2(res$correlation.sum)
plot(x, y)
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# Fit a linear model to the data in log-log scale. The slope should give
# us an estimate of the correlation dimension.
print(lm(y ~ x, data = data.frame(x = x, y = y)[2:18,]))
```

    ## 
    ## Call:
    ## lm(formula = y ~ x, data = data.frame(x = x, y = y)[2:18, ])
    ## 
    ## Coefficients:
    ## (Intercept)            x  
    ##      -0.931        1.227

``` r
# Cleanup, delete the searcher object.
release_searcher(searcher)
```

    ## ATRIA Destructor
    ## Total_clusters : 9313
    ## Total number of points in terminal nodes : 190687
    ## Average number of points in a terminal node : 40.9463
    ## Average percentage of points searched 0.467145% (935)
    ## Average number of terminal nodes visited : 38.941

    ## [1] TRUE

We also can use boxcouting to get a rough estimate of the capacity dimension:

``` r
# Next let's use the boxcounting approach to estimate the capacity dimension.
res <- boxcounting(data, dist.breaks)

x <- log2(res$dists)
y <- res$boxd[, 2]

plot(x, y)
```

![](README_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
print(lm(y ~ x, data = data.frame(x = x, y = y)[10:24,]))
```

    ## 
    ## Call:
    ## lm(formula = y ~ x, data = data.frame(x = x, y = y)[10:24, ])
    ## 
    ## Coefficients:
    ## (Intercept)            x  
    ##       -3.39         1.18

Author
------

Christian Merkwirth

License
-------

GPL (&gt;= 2)

References
----------

[Nearest-neighbor based methods for nonlinear time-series analysis, Christian Merkwirth, 2001](http://hdl.handle.net/11858/00-1735-0000-0006-B40F-A)
