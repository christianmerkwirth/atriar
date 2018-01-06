context("boxcounting")


test_that('boxcount gives appropriate output structure',{

})


test_that('boxcount gives right count of boxes',{
  for (d in c(4, 6, 8)) {
    x <- round(10 * matrix(runif(1000 * d), ncol=d))
    mode(x) <- "integer"

    # Let's just count non-empty bins in a D-dimensional space.
    y.unique <- unique(x);
    bc <- boxcount(x)

    atria.boxes.count <- 2^(-bc$boxd[d])
    unique.count <- nrow(y.unique)

    expect_equal(round(atria.boxes.count), unique.count)
  }
})
