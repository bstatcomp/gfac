#' Predict from CgPFA model
#'
#' This function predicts based on the estimation of the CgPFA
#' model \code{train_CgPFA}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_CgPFA}.
#' @param gp_test vector. Groups for each observation.
#' @param test data frame. Test data set for evaluation.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{lpd}{A numeric vector of log predictive densities for all
#'   observations in the test set.}
#' @seealso \code{\link{train_CgPFA}}
pred_CgPFA <- function (tr_fa, gp_test, test) {
  predictions <- matrix(data = NA, nrow=length(gp_test), ncol=ncol(tr_fa$train))
  gp_df <- unique(data.frame(orig = tr_fa$gp_train,
                             id   = as.integer(as.factor(as.character(
                               tr_fa$gp_train)))))
  gp_df_t <- left_join(data.frame(orig = gp_test), gp_df)
  ext     <- rstan::extract(tr_fa$stan_mod$PFA)
  nit     <- dim(ext$Theta)[1]
  log_lik <- vector(mode = "numeric", length = nrow(test))
  for (i in 1:nrow(test)) {
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
      mu_tmp  <- as.numeric(Psi_tmp %*% theta_tmp)
      int_fun <- function (x, Sigma) {
        val <- dmvnorm(qnorm(x), sigma = Sigma)
        jac <- prod(1 / (dnorm(qnorm(x))))
        return(val * jac)
      }
      mu_tmp  <- as.numeric(mu_tmp)
      n   <- 100
      set.seed(1)
      obs <- as.numeric(test[i, ])
      a   <- ppois(obs - 1, lambda = mu_tmp)
      b   <- ppois(obs, lambda = mu_tmp)
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
