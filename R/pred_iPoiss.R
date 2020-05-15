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
#'   \item{lpd}{A numeric vector of log predictive densities for
#'   all observations
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
    if (length(wh_g) > 0) { # if my_g in train data set
      predictions[i, ] <- round(tr_fa$params[wh_g, ])
      params           <- as.numeric(tr_fa$params[wh_g, ])
      alpha_tmp        <- tr_fa$alpha[wh_g, ]
      beta_tmp         <- tr_fa$beta[wh_g, ]
    } else {
      params           <- apply(tr_fa$params, 2, mean)
      predictions[i, ] <- round(params)
      params           <- as.numeric(params)
      mean_params <- matrix(params, ncol = length(params),
                            nrow = nrow(tr_fa$params), byrow = T)
      diff        <- tr_fa$params - mean_params
      diff        <- apply(diff, 1, sum)
      min_ind     <- which.min(diff)
      alpha_tmp   <- tr_fa$alpha[min_ind, ]
      beta_tmp    <- tr_fa$beta[min_ind, ]
    }
    tmp <- 0
    for (j in 1:1000) {
      lambda_tmp <- rgamma(ncol(test), alpha_tmp, beta_tmp)
      tmp        <- tmp + prod(dpois(my_vec, lambda_tmp, log = F))
    }
    tmp <- tmp / 1000
    log_lik[i] <- log(tmp)
  }
  return (list(predictions = predictions,
               lpd         = log_lik))
}
