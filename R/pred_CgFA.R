#' Predict from CgFA model
#'
#' This function predicts based on the estimation of the CgFA
#' model \code{train_CgFA}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_CgFA}.
#' @param gp_test vector. Groups for each observation in the test set.
#' @param test data frame. Test data set for evaluation.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{lpd}{A numeric vector of log predictive densities for all
#'   observations in the test set.}
#' @seealso \code{\link{train_CgFA}}
pred_CgFA <- function (tr_fa, gp_test, test) {
  if (is.na(tr_fa$Sigma)) {
    predictions <- NA
    log_lik     <- NA
  } else {
    predictions <- matrix(data = NA, nrow=length(gp_test), ncol=ncol(tr_fa$train))
    gp_df <- unique(data.frame(orig = tr_fa$gp_train,
                               id   = as.integer(as.factor(
                                 tr_fa$gp_train))))
    gp_df_t <- left_join(data.frame(orig = gp_test), gp_df)
    ext     <- rstan::extract(tr_fa$stan_mod$gFA)
    nit     <- dim(ext$L_tri)[1]
    log_lik <- vector(mode = "numeric", length = nrow(test))
    means   <- apply(tr_fa$train, MARGIN = 2, FUN = mean)
    sds     <- apply(tr_fa$train, MARGIN = 2, FUN = sd)
    test_st <- scale(test, center = means, scale = sds)
    for (i in 1:nrow(test)) {
      tmp    <- 0
      tmp_gp <- gp_test[i]
      my_mus <- matrix(data = NA, nrow = nit, ncol = ncol(test))
      for (j in 1:nit) {
        if (is.na(gp_df_t[i,2])) {
          fac_tmp <- apply(ext$iFactor[j, , ], 1, mean)
        } else {
          fac_tmp <- ext$iFactor[j, ,gp_df_t[i,2]]
        }
        L_tri_tmp <- ext$L_tri[j, , ]
        Psi       <- diag(ext$Psi[j, , ])
        n         <- 100
        set.seed(1)
        tmp_t     <- 0
        mu_tmp    <- as.numeric(L_tri_tmp %*% fac_tmp)
        obs <- as.numeric(test_st[i, ])
        a   <- pnorm(obs - 0.5 / sds, mean = mu_tmp, sd = sqrt(Psi))
        b   <- pnorm(obs + 0.5 / sds, mean = mu_tmp, sd = sqrt(Psi))
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
        tmp    <- tmp + tmp_t
        mu_tmp <- mu_tmp * sds + means
        my_mus[j, ] <- mu_tmp
      }
      tmp        <- tmp / nit
      log_lik[i] <- log(tmp)
      predictions[i, ] <- as.numeric(apply(my_mus, 2, mean))
    }
  }
  return (list(predictions = predictions,
               lpd         = log_lik))
}
