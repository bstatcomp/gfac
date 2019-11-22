#' Calculate the likelihood for the FA and gFA models
#'
#' This function calculates the MVN likelihood for factor analysis methods.
#' If \code{L_intg = FALSE} then an approximation with point-density is used
#' instead of integrating over suitable bounds.
#'
#' @export
#' @param test data frame. Test data set for evaluation.
#' @param pred data frame or matrix. The provided predictions.
#' @param tr_fa Test dataset (for evaluation).
#' @param logl logical. Should the output be log-likelihood?
#' @param L_intg logical. Should the likelihood be integrated over
#' suitable borders? Defaults to FALSE, which uses an approximation with
#' point-MVN density.
#' Special care needs to be taken when integrating, as the likelihood is
#' likely to be 0 in higher dimensions due to machine precision, resulting in
#' non-finite log-likelihood.
#' @return numeric vector. The (log-)likelihood for each test observation.
#' @seealso \code{\link{my_dmvn}} \code{\link{my_omxMnor}}
likelihood_cl <- function (test, pred, tr_fa, logl = TRUE, integrate = FALSE,
                           truncate = FALSE) {
  out     <- vector(mode = "numeric", length = nrow(test))
  sigma   <- as.vector(t(tr_fa$noise))
  l       <- length(sigma)
  n       <- nrow(pred)
  x       <- as.matrix(test)
  sigma_m <- matrix(rep(sigma, times = n), nrow = n, byrow = TRUE)
  if (integrate) {
    my_data <- cbind(sigma_m, pred, x - 0.5, x + 0.5)
    out     <- apply(my_data, MARGIN = 1, FUN = my_omxMnor, l = sqrt(l),
                     logl = logl)
  } else {
    my_data <- cbind(x, pred, sigma_m)
    out <- apply(my_data, MARGIN = 1, FUN = my_dmvn, l = sqrt(l), logl = logl)
  }
  return (out)
}
