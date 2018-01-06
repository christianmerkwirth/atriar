context("boxcounting")


test_that('boxcount computes sample cases right',{
  # All points in one box.
  d <- 4
  x <- matrix(0.0, nrow = 1000, ncol=d)
  mode(x) <- "integer"
  bc <- boxcount(x)
  expect_equal(bc$boxes, rep(1.0, d))
  expect_equal(bc$entropy, rep(0.0, d))
  expect_equal(bc$correlation, rep(1.0, d))

  # All boxes filled with one points.
  s <- 1:4
  x <- as.matrix(expand.grid(s, s, s, s, s))
  mode(x) <- "integer"
  bc <- boxcount(x)
  expect_equal(bc$boxes, 4^(1:5))
  expect_equal(bc$entropy, seq(2, 10, by=2))
  expect_equal(bc$correlation, 4^(-(1:5)), tolerance=1e-5)
})

test_that('boxcount gives right count of boxes',{
  for (d in c(4, 6, 8)) {
    x <- round(10 * matrix(runif(1000 * d), ncol=d))
    mode(x) <- "integer"

    # Let's just count non-empty bins in a D-dimensional space.
    y.unique <- unique(x);
    bc <- boxcount(x)

    atria.boxes.count <-bc$boxes[d]
    unique.count <- nrow(y.unique)

    expect_equal(round(atria.boxes.count), unique.count)
  }
})
