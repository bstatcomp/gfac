#' Predict from FA model
#'
#' This function predicts based on the estimation of the FA
#' model \code{train_FA}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_FA}.
#' @param gp_test vector. Groups for each observation.
#' @param test data frame. Test data set for evaluation.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{log_lik}{A numeric vector of (log-)likelihoods for all observations
#'   in the test set.}
#' @seealso \code{\link{train_FA}}
pred_FA <- function (tr_fa, gp_test, test) {
  predictions <- matrix(data = NA, nrow=length(gp_test), ncol=ncol(tr_fa$train))
  gp_df <- unique(data.frame(orig = tr_fa$gp_train,
                             id   = as.integer(as.factor(
                               tr_fa$gp_train))))
  train_ids <- as.integer(as.factor(tr_fa$gp_train))
  gp_df_t <- left_join(data.frame(orig = gp_test), gp_df)
  ext     <- rstan::extract(tr_fa$stan_mod)
  nit     <- dim(ext$L_tri)[1]
  log_lik <- vector(mode = "numeric", length = nrow(test))
  means   <- apply(tr_fa$train, MARGIN = 2, FUN = mean)
  sds     <- apply(tr_fa$train, MARGIN = 2, FUN = sd)
  for (i in 1:nrow(test)) {
    tmp    <- 0
    tmp_gp <- gp_test[i]
    my_mus <- matrix(data = NA, nrow = nit, ncol = ncol(test))
    for (j in 1:nit) {
      if (!is.null(ext$iFactor)) {
        if (is.na(gp_df_t[i,2])) {
          fac_tmp <- apply(ext$iFactor[j, , ], 1, mean)
        } else {
          fac_tmp <- ext$iFactor[j, ,gp_df_t[i,2]]
        }
      } else {
        if (is.na(gp_df_t[i,2])) {
          fac_tmp <- apply(ext$Factor[j, , ], 1, mean)
        } else {
          factor_gp_id <- gp_df_t[i,2]
          factor_gp    <- ext$Factor[j, ,train_ids == factor_gp_id]
          fac_tmp      <- apply(factor_gp, 1, mean)
        }
      }
      L_tri_tmp <- ext$L_tri[j, , ]
      Psi       <- diag(ext$Psi[j, , ])
      n         <- 100
      tmp_t     <- 0
      mu_tmp    <- as.numeric(L_tri_tmp %*% fac_tmp)

      # de-standardize
      mu_tmp <- mu_tmp * sds + means
      Psi    <- Psi * sds^2

      set.seed(1)
      U <- runif(n * ncol(test), as.numeric(test[i, ]) - 0.5, as.numeric(test[i, ]) + 0.5)
      U <- matrix(U, ncol = ncol(test), byrow = T)
      tmp_t <- dnorm(t(U), mean = mu_tmp, sd = sqrt(Psi), log = F)
      tmp_t <- t(tmp_t)
      tmp_t <- apply(tmp_t, 1, prod)
      tmp_t <- mean(tmp_t)
      tmp   <- tmp + tmp_t
      my_mus[j, ] <- mu_tmp
    }
    tmp        <- tmp / nit
    log_lik[i] <- log(tmp)
    predictions[i, ] <- as.numeric(apply(my_mus, 2, mean))
  }
  return (list(predictions  = predictions,
               lpd          = log_lik))
}
