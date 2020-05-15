#' Predict from gPFA model
#'
#' This function predicts based on the estimation of the gPFA
#' model \code{train_gPFA}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_gPFA}.
#' @param gp_test vector. Groups for each observation.
#' @param test data frame. Test data set for evaluation.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{lpd}{A numeric vector of log predictive densities for all
#'   observations in the test set.}
#' @seealso \code{\link{train_gPFA}}
pred_gPFA <- function (tr_fa, gp_test, test) {
  predictions <- matrix(data = NA, nrow=length(gp_test), ncol=ncol(tr_fa$train))
  gp_df <- unique(data.frame(orig = tr_fa$gp_train,
                             id   = as.integer(as.factor(as.character(
                               tr_fa$gp_train)))))
  gp_df_t <- left_join(data.frame(orig = gp_test), gp_df)
  ext     <- rstan::extract(tr_fa$stan_mod)
  nit     <- dim(ext$Theta)[1]
  log_lik <- vector(mode = "numeric", length = nrow(test))
  for (i in 1:nrow(test)) {
    print(paste0("Test: ", i))
    tmp    <- 0
    tmp_gp <- gp_test[i]
    my_mus <- matrix(data = NA, nrow = nit, ncol = ncol(test))
    for (j in 1:nit) {
      if (is.na(gp_df_t[i,2])) {
        theta_tmp <- apply(ext$Theta[j, , ], 1, mean)
      } else {
        theta_tmp <- ext$Theta[j, ,gp_df_t[i,2]]
      }
      Psi_tmp <- ext$Psi[j, , ]
      mu_tmp  <- Psi_tmp %*% theta_tmp
      tmp     <- tmp + prod(dpois(as.numeric(test[i, ]),
                                  lambda = as.numeric(mu_tmp),
                                  log = F))
      my_mus[j, ] <- mu_tmp
    }
    tmp        <- tmp / nit
    log_lik[i] <- log(tmp)
    predictions[i, ] <- as.numeric(apply(my_mus, 2, mean))
  }

  return (list(predictions = predictions,
               lpd         = log_lik))
}
