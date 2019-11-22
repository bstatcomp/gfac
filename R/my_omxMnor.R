#' Calculate the integrated density
#'
#' This function calculates the MVN density by integrating along
#' corresponding bounds.
#'
#' @export
#' @param my_vec numeric vector. First \code{l}^2 values represent the
#' estimated residual covariance, the next \code{l} values represent
#' the predictions (means), the next \code{l} values the lower bounds for
#' integration, and the last \code{l} values the upper bounds for integration.
#' @param l numeric. Number of dimensions.
#' @param logl logical. Calculate log density? Defaults to true.
#' @return numeric vector. The (log) density.
my_omxMnor <- function (my_vec, l, logl) {
  out <- omxMnor(covariance = matrix(my_vec[1:l^2], nrow = l, byrow = T),
                 means      = my_vec[(l^2+1):(l^2+l)],
                 lbound     = my_vec[(l^2+l+1):(l^2+2*l)],
                 ubound     = my_vec[(l^2+2*l+1):(l^2+3*l)])
  if (logl) out <- log(out)
  return (out)
}
