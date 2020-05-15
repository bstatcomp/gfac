#' Predict from CgNBFA model
#'
#' This function predicts based on the estimation of the CgNBFA
#' model \code{train_CgNBFA}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_CgNBFA}.
#' @param gp_test vector. Groups for each observation.
#' @param test data frame. Test data set for evaluation.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{lpd}{A numeric vector of log predictive densities for all
#'   observations in the test set.}
#' @seealso \code{\link{train_CgNBFA}}
pred_CgNBFA <- function (tr_fa, gp_test, test) {
  predictions <- matrix(data = NA, nrow=length(gp_test), ncol=ncol(tr_fa$train))
  gp_df <- unique(data.frame(orig = tr_fa$gp_train,
                             id   = as.integer(as.factor(as.character(
                               tr_fa$gp_train)))))
  gp_df_t <- left_join(data.frame(orig = gp_test), gp_df)
  ext     <- rstan::extract(tr_fa$stan_mod$NBFA)
  nit     <- dim(ext$mu)[1]
  log_lik <- vector(mode = "numeric", length = nrow(test))
  for (i in 1:nrow(test)) {
    tmp    <- 0
    tmp_gp <- gp_test[i]
    my_mus <- matrix(data = NA, nrow = nit, ncol = ncol(test))
    for (j in 1:nit) {
      if (is.na(gp_df_t[i,2])) {
        mu_tmp <- apply(ext$mu[j, , ], 1, mean)
      } else {
        mu_tmp  <- ext$mu[j, ,gp_df_t[i,2]]
      }
      phi_tmp <- ext$phi[j, ]
      int_fun <- function (x, Sigma) {
        val <- dmvnorm(qnorm(x), sigma = Sigma)
        jac <- prod(1 / (dnorm(qnorm(x))))
        return(val * jac)
      }
      mu_tmp  <- as.numeric(mu_tmp)
      phi_tmp <- as.numeric(phi_tmp)
      n   <- 100
      set.seed(1)
      obs <- as.numeric(test[i, ])
      a   <- pnbinom(obs - 1, size = phi_tmp, mu = mu_tmp)
      b   <- pnbinom(obs, size = phi_tmp, mu = mu_tmp)
      U   <- runif(n * length(obs), a, b)
      U   <- matrix(U, ncol = ncol(test), byrow = T)
      qU  <- qnorm(U)
      mU  <- dmvnorm(qU, sigma = tr_fa$Sigma)
      jac <- 1 / dnorm(qnorm(U))
      jac <- apply(jac, 1, prod)
      tmp_t <- mU * jac
      fst   <- prod(abs(b - a))
      tmp_t <- fst * tmp_t
      tmp_t <- mean(tmp_t)
      if (is.nan(tmp_t)) {
        tmp_t <- 0
      }
      tmp         <- tmp + tmp_t
      my_mus[j, ] <- mu_tmp
    }
    tmp        <- tmp / nit
    log_lik[i] <- log(tmp)
    predictions[i, ] <- as.numeric(apply(my_mus, 2, mean))
  }
  return (list(predictions  = predictions,
               lpd          = log_lik))
}
