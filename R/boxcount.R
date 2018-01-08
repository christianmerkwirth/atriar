#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param dist.breaks PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname boxcounting
boxcounting <- function(data, dist.breaks) {
  boxes <- numeric()
  entropy <- numeric()
  correlation <- numeric()
  for (dist in dist.breaks) {
    x <- data / dist
    mode(x) <- 'integer'
    bc <- boxcount(x)
    if (length(boxes) == 0) {
      boxes <- bc$boxes
      entropy <- bc$entropy
      correlation <- bc$correlation
    } else {
      boxes <- rbind(boxes, bc$boxes)
      entropy <- rbind(entropy, bc$entropy)
      correlation <- rbind(correlation, bc$correlation)
    }
  }
  return(
    list(
      dists = dist.breaks,
      boxes = boxes,
      entropy = entropy,
      correlation = correlation
    ))
}
