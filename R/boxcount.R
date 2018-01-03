boxcounting <- function(data, dist.breaks) {
  boxd <- numeric()
  infod <- numeric()
  corrd <- numeric()
  for (dist in dist.breaks) {
    x <- data / dist
    mode(x) <- 'integer'
    bc <- boxcount(x)
    if (length(boxd) == 0) {
      boxd <- bc$boxd
      infod <- bc$infod
      corrd <- bc$corrd
    } else {
      boxd <- rbind(boxd, bc$boxd)
      infod <- rbind(infod, bc$infod)
      corrd <- rbind(corrd, bc$corrd)
    }
  }
  return(
    list(
      dists = dist.breaks,
      boxd = boxd,
      infod = infod,
      corrd = corrd
    ))
}
