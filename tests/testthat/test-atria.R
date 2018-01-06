context("atria nearest neighbors")


test_that('atria gives appropriate result structure',{
  # externalpointer
})

# release searcher returns true

# search k-NN

# search range

# count agrees with length of index and dist

# distances are correct in different metrics

test_that('atria and RANN:nn2 agree',{
  if(require('RANN')){
    for (d in c(4, 6, 8)) {
      k <- 8
      train <- matrix(rnorm(1000*d), ncol=d)
      test <- matrix(rnorm(100*d), ncol=d)
      nn.rann <- nn2(data=train, query=test, k=k, eps=0.0)
      searcher <- create_searcher(train, metric='euclidian')
      nn.atria <- search_k_neighbors(
        searcher = searcher,
        k = k,
        query_points = test,
        epsilon = 0.0)
      release_searcher(searcher)
      expect_equal(nn.atria$index, nn.rann$nn.idx)
      expect_equal(nn.atria$dist, nn.rann$nn.dists, tolerance=1e-4)
    }
  }
})

