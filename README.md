# Overview of the Rcorrsum R package


Christian Merkwirth

2017-12-30

Fast code for computing statistics of pairwise distances in the a data set, which can be used to estimate the 
correlation dimension and related properties. The correlation dimension is one of the intrinsic dimensions of a data set and can be viewed as a measure of its complexity. It is expected that the performance of many machine learning algorithms is dependent on the intrinsic dimension of the input data.

For more info on the correlation sum, see (https://en.wikipedia.org/wiki/Correlation_integral).

This package contains an implementation of the ATRIA nearest neighbor search algorithm which is a variant of the ball-tree algorithm family. It supports both exact and approximate k-nearest neighbor and range searches. For
details, see (https://www.researchgate.net/publication/12238515_Fast_nearest-neighbor_searching_for_nonlinear_signal_processing).

Part of the code was released before within the OpenTSTOOL Matlab toolbox by the DPI Goettingen, Germany.


Chapter 4
Nearest Neighbors Searching
An integral part of a majority of methods for nonlinear time series analysis is searching for nearest
neighbors. The perfomance of these methods depends strongly of the perfomance of the employed
nearest neighbor algorithm. Thus, choosing an efficient nearest neighbor algorithm should be done
very carefully.
4.1  Definition
Definition : A set
P
of data points in D-dimensional space is given.  Then we define the nearest
neighbor to some reference point
q
(also called
query point
) to be the point of data set
P
that has
the smallest distance to
q
(we don’t issue the problem of ambiguity at this point).The more general
task of finding more than one nearest neighbor is called
k nearest neighbors problem
. In general, the
reference point
q
is an arbitrarily located point, but it is also possible that
q
is itself a member of
data set
P
(as illustrated in the figure, where five neighbors to
q
(excluding self match) are found).
4.2  Approximate nearest neighbors searching
Approximate nearest neighbors algorithms report neighbors to the query point
q
with distances pos-
sibly greater than the true nearest neighbors distances. The maximal allowed relative error
epsilon
is
given as a parameter to the algorithm. For
epsilon
=0, the approximate search returns the true (exact)
nearest neighbor(s).
Computing exact nearest neighbors for data set with
fractal dimension
much higher than 6 seems to
be a very time-consuming task. Few algorithms seem to perform significantly better than a brute-
force computation of all distances. However, it has been shown that by computing nearest neighbors
approximately, it is possible to achieve significantly faster execution times with relatively small actual
errors in the reported distances.
19
4.3  Range searching
In the task of
range searching
, we ask for all points of data set
P
that have distance
r
or less from the
query point
q
. Sometimes
range searching
is called a
fixed size approach
, while
k nearest neighbors
searching
is called a
fixed mass approach
.
4.4  Matlab mex-functions
4.4.1  nn
prepare
nn
prepare
does the preprocessing for a given data set
pointset
. The returned data structure
atria
contains preprocessing information that is necessary to use
nn
search
or
range
search
.
Preprocessing and searching is divided into different mex-files to give the user the possibility to
re-
use
the preprocessing data (contained in
atria
) when doing multiple searches on the same point set.
However, as soon as the underlying point set is changed or modified, one has to recompute
atria
for
the changed point set.
Syntax:
•
atria = nn_prepare(pointset)
•
atria = nn_prepare(pointset, metric)
•
atria = nn_prepare(pointset, metric, clustersize)
Input arguments:
•
pointset
- a
N
by
D
double matrix containing the coordinates of the point set, organized as
N
points of dimension
D
•
metric
- (optional) either ’euclidian’ or ’maximum’ (default is ’euclidian’)
•
clustersize
- (optional) threshold for clustering algorithm, defaults to 64
4.4.2  nn
search
nn
search
does exact or approximate k-nearest neighbor queries to one or more query points. These
query points can be given explicitly or taken from the data set of points (see below).
Before one can use
nn
search
, one has to call
nn
prepare
to compute the preprocessing information.
However, as long as the input point set isn’t modified, the preprocessing information is valid and can
be re-used for multiple calls to
nn
search
or
range
search
.
Syntax:
20
•
[index, distance] = nn_search(pointset, atria, query_points, k)
•
[index, distance] = nn_search(pointset, atria, query_points, k, epsilon)
•
[index, distance] = nn_search(pointset, atria, query_indices, k, exclude)
•
[index, distance] = nn_search(pointset, atria, query_indices, k,
exclude, epsilon)
Input arguments:
•
pointset
- a
N
by
D
double matrix containing the coordinates of the point set, organized as
N
points of dimension
D
•
atria
- output of nn
prepare for pointset
•
query
points
- a
R
by
D
double matrix containing the coordinates of the query points, orga-
nized as
R
points of dimension
D
•
query
indices
- query points are taken out of the pointset, query
indices is a vector of length
R
which contains the indices of the query points (indices may vary from 1 to N)
•
k
- number of nearest neighbors to compute
•
epsilon
- (optional) relative error for approximate nearest neighbors queries, defaults to 0 (=
exact search)
•
exclude
- in case the query points are taken out of the pointset, exclude specifies a range of
indices which are omitted from search. For example if the index of the query point is 124 and
exclude is set to 3, points with indices 121 to 127 are omitted from search. Using exclude = 0
means: exclude self-matches
Output arguments:
•
index
- a matrix of size
R
by
k
which contains the indices of the nearest neighbors. Each row
of
index
contains
k
indices of the nearest neighbors to the corresponding query point.
•
distance
- a matrix of size
R
by
k
which contains the distances of the nearest neighbors to the
corresponding query points, sorted in increasing order.
4.4.3  range
search
The routine
range
search
does a range search to one or more query points. These query points can
be given explicitly or taken from the data set of points (see below).
Before one can use
range
search
, one has to call
nn
prepare
to compute the preprocessing informa-
tion. However, as long as the input point set isn’t modified, the preprocessing information is valid
and can be re-used for multiple calls to
nn
search
or
range
search
.
Syntax:
•
[count, neighbors] = range_search(pointset, atria, query_points, r)
•
[count, neighbors] = range_search(pointset, atria, query_indices, r, exclude)
Input arguments:
•
pointset
- a
N
by
D
double matrix containing the coordinates of the point set, organized as
N
points of dimension
D
21
•
query
points
- a
R
by
D
double matrix containing the coordinates of the query points, orga-
nized as
R
points of dimension
D
•
query
indices
- query points are taken out of the pointset, query
indices is a vector of length
R
which contains the indices of the query points
•
r
- range or search radius (r
>
0)
•
exclude
- in case the query points are taken out of the pointset, exclude specifies a range of
indices which are omitted from search. For example if the index of the query point is 124 and
exclude is set to 3, points with indices 121 to 127 are omitted from search. Using exclude = 0
means: exclude self-matches
Output arguments:
•
count
- a vector of length
R
contains the number of points within distance r to the corresponding
query point
•
neighbors
- a Matlab cell structure of size
R
by
2
which contains vectors of indices and vectors
of distances to the neighbors for each given query point. This output argument can not be stored
in a standard Matlab matrix because the number of neighbors within distance r is not the same
for all query points. The vectors if indices and distances for one query point have exactly the
length that is given in count. The values in the distances vectors are
not
sorted.
4.5  Example session
% create a 3-dimensional data set with 100000 points
pointset = rand(100000, 3);
% do the preprocessing for this point set
atria = nn_prepare(pointset, ’euclidian’);
% now search for 2 (exact) nearest neighbors, using points 1 to
% 10 as query points, excluding self-matches
[index, distance] = nn_search(pointset, atria, 1:10, 2, 0)
index =
5618       96574
38209       84549
54991       60397
38429       59732
4114       76991
72121         452
13678       59332
26022       16718
86042       38436
24830       44434
distance =
22
0.0101    0.0175
0.0078    0.0134
0.0132    0.0167
0.0050    0.0223
0.0087    0.0097
0.0124    0.0189
0.0129    0.0168
0.0046    0.0110
0.0101    0.0103
0.0156    0.0177
% now do a range search for radius 0.0224, using points 1 to 10 as
% query points, excluding self-matches
[count, neighbors] = range_search(pointset, atria, 1:10, 0.0224, 0)
count =
4
10
7
2
5
6
2
4
7
5
neighbors =
[1x4  double]    [1x4  double]
[1x10 double]    [1x10 double]
[1x7  double]    [1x7  double]
[1x2  double]    [1x2  double]
[1x5  double]    [1x5  double]
[1x6  double]    [1x6  double]
[1x2  double]    [1x2  double]
[1x4  double]    [1x4  double]
[1x7  double]    [1x7  double]
[1x5  double]    [1x5  double]
% let’s see the indices of the points that are within range to the first query point
neighbors{1,1}
ans =
56921       97100       96574        5618
% let’s see the corresponding distances of the points that are
% within range to the first query point
neighbors{1,2}
23



# References

Nearest-neighbor based methods for nonlinear time-series analysis, Christian Merkwirth, 2001 (http://hdl.handle.net/11858/00-1735-0000-0006-B40F-A)







