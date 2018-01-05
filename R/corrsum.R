#' @useDynLib atriar
#' @importFrom Rcpp sourceCpp
NULL

distlimits <- function(searcher, data) {
  N <- nrow(data)
  k.max <- ceiling(0.10 * N)

  # First get a rough idea about distances in the data set.
  rand.sample <- sample.int(N, size = 64)
  nn <- search_k_neighbors(
    searcher = searcher,
    k = k.max,
    query_points = data[rand.sample, ],
    exclude = cbind(rand.sample, rand.sample)
  )

  min.dist <- 0
  d <- 1
  while ((min.dist == 0) && (d < ncol(nn$dist))) {
    min.dist <- quantile(nn$dist[, d], probs = 0.20)
    d <- d + 1
  }
  if (min.dist == 0) {
    stop('Can not determine a minimum lenght scale.')
  }
  # This translates to a dist that should cover 10% of all pairs.
  max.dist <- quantile(nn$dist[, k.max], probs = 0.50)
  return(c(min.dist, max.dist))
}

logspace <- function(min.dist, max.dist, n.bins) {
  exp(seq(log(min.dist), log(max.dist), length.out=n.bins))
}

corrsum <-
  function(searcher,
           data,
           dist.breaks,
           min.actual.pairs = 2000,
           min.nr.samples.at.scale = 128,
           max.nr.samples.at.scale = 1024,
           batch.size = 32,
           verbose = FALSE) {
    N <- nrow(data)
    dists <- c(0, dist.breaks)
    actual.pairs.count <- numeric(length(dists))
    potential.pairs.count <- numeric(length(dists))
    samples.used <- numeric(length(dists))

    # We iterate from higher to lower distance scale. We switch to smaller length scale as
    # soon as the statistics is good enough. Good enough means enough count in both numerator
    # and denominator.
    for (pos in rev(1:length(dists))) {
      radius <- dists[pos]
      if (radius == 0) {
        break
      }
      # Perform at least one range search at each distance bin.
      while ((actual.pairs.count[pos] < min.actual.pairs) ||
             (samples.used[pos] < min.nr.samples.at.scale)) {
        if (samples.used[pos] >= max.nr.samples.at.scale) {
          break
        }
        rand.sample <- sample.int(N, size = batch.size)
        nn <- search_range(
          searcher = searcher,
          radius = radius,
          query_points = data[rand.sample, ],
          # Ignore samples with index smaller than the query point
          # so we avoid counting the same pairwise distance twice.
          exclude = cbind(rep(-1, batch.size), rand.sample)
        )
        # Calculate the potential number of pairs.
        potential.pairs <- sum(N - rand.sample)
        potential.pairs.count[dists <= radius] <-
          potential.pairs.count[dists <= radius] + potential.pairs
        new.pair.counts <- sapply(nn$nn,
                                  function(x) {
                                    q <- .bincode(x$dist, breaks = c(0, dists), TRUE, TRUE)
                                    #ścat(q)
                                    count_integers(q, length(dists))
                                  })
        actual.pairs.count <-
          actual.pairs.count + rowSums(new.pair.counts)
        samples.used[pos] <- samples.used[pos] + length(nn$count)
        if (verbose) {
          cat(paste0(
            'Radius: ',
            radius,
            ' count: ',
            sum(nn$count),
            ' ',
            sum(rowSums(new.pair.counts)),
            '\n'
          ))
        }
      }
    }
    return(
      list(
        dists = dists,
        correlation.sum = cumsum(actual.pairs.count / potential.pairs.count),
        actual.pairs.count = actual.pairs.count,
        potential.pairs.count = potential.pairs.count,
        samples.used = samples.used
      )
    )
  }
