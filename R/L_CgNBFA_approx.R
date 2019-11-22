#' Compute approximate copula likelihood with negative binomial marginals
#'
#' This function computes the approximated likelihood of a Gaussian copula with
#' NB marginals. Integration in the calculation of the exact
#' likelihood can be computationally expensive, especially in high dimensions,
#' therefore this function provides us with an approximation of the integral.
#'
#' @export
#' @param obs vector. The data point.
#' @param size vector. The estimated size parameter.
#' @param mu vector. The estimated mean parameter.
#' @param Gamma matrix. The estimated correlation matrix.
#' @param logl logical. Should the output be log-likelihood?
#' @return The (log-)likelihood.
L_CgNBFA_approx <- function (obs, size, mu, Gamma, logl = T) {
  obs    <- as.numeric(obs)
  a      <- pnbinom(obs - 1, size = size, mu = mu)
  b      <- pnbinom(obs, size = size, mu = mu)
  int_fun <- function (x, Sigma) {
    val <- dmvnorm(qnorm(x), sigma = Sigma)
    jac <- prod(1 / (dnorm(qnorm(x))))
    return(val * jac)
  }
  if (logl) {
    out <- sum(log(abs(b - a))) + log(int_fun((b + a) / 2, Gamma))
  } else {
    out <- prod(abs(b - a)) * int_fun((b + a) / 2, Gamma)
  }
  return (out)
}
