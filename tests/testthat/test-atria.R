context("atria nearest neighbors")


search_knn <-
  function(train,
           test,
           k,
           metric = 'euclidian',
           epsilon = 0.0) {
    searcher <- create_searcher(train, metric = metric)
    nn.atria <- search_k_neighbors(
      searcher = searcher,
      k = k,
      query_points = test,
      epsilon = epsilon
    )
    release_searcher(searcher)
    return(nn.atria)
  }

search_radius <- function(train, test, radius, metric = 'euclidian') {
  searcher <- create_searcher(train, metric = metric)
  nn.atria <- search_range(searcher = searcher,
                           radius = radius,
                           query_points = test)
  release_searcher(searcher)
  return(nn.atria)
}

eucl.dist <- function(x, y) {
  return(sqrt(sum((x - y) ^ 2)))
}

man.dist <- function(x, y) {
  return(sum((abs(x - y))))
}

max.dist <- function(x, y) {
  return(max((abs(x - y))))
}

check.distances <- function(nn, train, test, dist.func) {
  for (i in 1:nrow(nn$index)) {
    for (j in 1:ncol(nn$index)) {
      index <- nn$index[i, j]
      dist <- nn$dist[i, j]
      x <- train[index,]
      y <- test[i,]
      expected.dist <- dist.func(x, y)
      expect_equal(dist, expected.dist, tolerance = 1e-3)
    }
  }
}

test_that('create and destruct searcher', {
  d <- 5
  train <- matrix(rnorm(1000 * d), ncol = d)
  searcher <- create_searcher(train, metric = 'manhattan')
  expect_equal(class(searcher), 'externalptr')
  expect_equal(number_of_points(searcher), nrow(train))
  expect_true(release_searcher(searcher))
})

test_that('search range', {
  d <- 4
  radius <- 1.2
  train <- matrix(rnorm(1000 * d), ncol = d)
  test <- matrix(rnorm(20 * d), ncol = d)
  nn.atria <- search_radius(train, test, radius, metric = 'manhattan')

  for (i in 1:length(nn.atria$nn)) {
    nn <- nn.atria$nn[[i]]
    # Check that count agrees with length of index and dist
    expect_equal(length(nn$index), nn.atria$count[i])
    expect_equal(length(nn$dist), nn.atria$count[i])
    if (length(nn$index)) {
      for (j in 1:length(nn$index)) {
        index <- nn$index[j]
        dist <- nn$dist[j]
        x <- train[index,]
        y <- test[i,]
        expected.dist <- man.dist(x, y)
        expect_equal(dist, expected.dist, tolerance = 1e-3)
      }
    }
  }
})

test_that('distances are correct in different metrics', {
  k <- 5
  for (d in c(4, 6, 8)) {
    train <- matrix(rnorm(1000 * d), ncol = d)
    test <- matrix(rnorm(20 * d), ncol = d)
    nn.eucl <- search_knn(train, test, k, metric = 'euclidian')
    nn.man <- search_knn(train, test, k, metric = 'manhattan')
    nn.max <- search_knn(train, test, k, metric = 'maximum')
    check.distances(nn.eucl, train, test, eucl.dist)
    check.distances(nn.man, train, test, man.dist)
    check.distances(nn.max, train, test, max.dist)
  }
})

test_that('atria and RANN:nn2 agree', {
  if (require('RANN')) {
    for (d in c(4, 6, 8)) {
      k <- 8
      train <- matrix(rnorm(1000 * d), ncol = d)
      test <- matrix(rnorm(100 * d), ncol = d)
      nn.rann <- nn2(
        data = train,
        query = test,
        k = k,
        eps = 0.0
      )
      nn.atria <- search_knn(train, test, k = k, epsilon = 0.0)
      expect_equal(nn.atria$index, nn.rann$nn.idx)
      expect_equal(nn.atria$dist, nn.rann$nn.dists, tolerance = 1e-4)
    }
  }
})
