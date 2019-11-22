#' Calculate the density
#'
#' This function calculates the MVN density.
#'
#' @export
#' @param my_vec numeric vector. First \code{l} values represent the
#' observations, the next \code{l} values represent the predictions (means),
#' and
#' the last \code{l}^2 values represent the estimated residual covariance
#' matrix.
#' @param l numeric. Number of dimensions.
#' @param logl logical. Calculate log density? Defaults to true.
#' @return numeric vector. The (log) density.
my_dmvn    <- function (my_vec, l, logl) {
  return (dtmvnorm(x     = my_vec[1:l],
                   mean  = my_vec[(l + 1):(2*l)],
                   sigma = matrix(my_vec[(2*l + 1):(2*l+l^2)], nrow = l,
                                  byrow = T),
                   log   = logl))
}
