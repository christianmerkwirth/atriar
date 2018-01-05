library(RcppAnnoy)
library(RANN)
library(microbenchmark)
library(h5)
library(atriar)

file <- h5file('/Users/cmerk/Downloads/mnist-784-euclidean.hdf5', mode='r')
train <- file['train'][]
test <- file['test'][]
expected.distances <- file['distances']
expected.neighbors <- file['neighbors']

k <- 10

work <- function(max.points, eps=0) {
  s <- create_searcher(train,
                       metric = 'manhattan',
                       cluster_max_points = max.points,
                       seed=ceiling(runif(1)*1e9))
  nn.atria <- search_k_neighbors(searcher = s,
                                 k = k,
                                 query_points = test[1:2000, ],
                                 epsilon = eps)
  release_searcher(s)
  return(nn.atria)
}


mb <- microbenchmark(
  work(8, 0),
  work(64, 0),
  work(256, 0),
  work(512, 0),
  work(1024, 0),
  work(8, 2.0),
  work(64, 2.0),
  work(256, 2.0),
  work(512, 2.0),
  work(1024, 2.0),
  times = 1
)

print(mb)




if (0) {
a <- new(AnnoyEuclidean, ncol(train))
for (i in seq(nrow(train))) {
  a$addItem(i-1, train[i, ])
}
a$build(50)                           	# 50 trees

k <- 10

nn.annoy <- apply(test, 1, function(x) { a$getNNsByVector(x, k) } )
}
