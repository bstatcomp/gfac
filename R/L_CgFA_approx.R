#' Compute approximate copula likelihood with normal marginals
#'
#' This function computes the approximated likelihood of a Gaussian copula with
#' normal marginals. Integration in the calculation of the exact
#' likelihood can be computationally expensive, especially in high dimensions,
#' therefore this function provides us with an approximation of the integral.
#'
#' @export
#' @param obs vector. The data point.
#' @param mu vector. The estimated mean parameter.
#' @param sd vector. The estimated standard deviation.
#' @param Gamma matrix. The estimated correlation matrix.
#' @param logl logical. Should the output be log-likelihood?
#' @param sd_train vector. The standard deviation of the original
#' (unstandardized) data.
#' @return The (log-)likelihood.
L_CgFA_approx <- function (obs, mu, sd, Gamma, logl = T, sd_train) {
  obs    <- as.numeric(obs)
  a      <- pnorm(obs - 0.5 / sd_train, mean = mu, sd = sd)
  b      <- pnorm(obs + 0.5 / sd_train, mean = mu, sd = sd)
  int_fun <- function (x, Sigma) {
    val <- dmvnorm(qnorm(x), sigma = Sigma)
    jac <- prod(1 / (dnorm(qnorm(x))))
    return(val * jac)
  }
  out <- sum(log(abs(b - a))) + log(int_fun((b + a) / 2, Gamma))
  return (out)
}
