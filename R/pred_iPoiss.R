#' Predict from an independent Poisson model
#'
#' This function predicts based on the estimation of the independent Poisson
#' model \code{train_iPoiss}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_iPoiss()}.
#' @param gp_test vector. Groups for each observation in the test set.
#' @param test data frame. Test data set for evaluation.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{log_lik}{A numeric vector of log-likelihoods for all observations
#'   in the test set.}
#' @seealso \code{\link{train_iPoiss}}
pred_iPoiss <- function (tr_fa, gp_test, test) {
  uniq_tr     <- unique(tr_fa$gp_train)
  predictions <- test
  log_lik     <- vector(mode = "numeric", length = nrow(test))
  for (i in 1L:length(gp_test)) {
    my_g             <- gp_test[i]
    my_vec           <- as.numeric(test[i, ])
    wh_g             <- which(my_g == uniq_tr)
    if (length(wh_g) > 0) {
      predictions[i, ] <- round(tr_fa$params[wh_g, ])
      params           <- as.numeric(tr_fa$params[wh_g, ])
    } else {
      params           <- apply(tr_fa$params, 2, mean)
      predictions[i, ] <- round(params)
      params           <- as.numeric(params)
    }
    log_lik[i] <- sum(dpois(my_vec, params, log = TRUE))
  }
  return (list(predictions = predictions,
               log_lik     = log_lik))
}
