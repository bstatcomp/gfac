#' Train a copula FA model
#'
#' This function trains a copula factor analysis model in Stan.
#'
#' The transformations are made separately, between the modeling of
#' marginals and the correlation.
#'
#' @export
#' @param train data frame. The training data set, consisting of counts.
#'   Columns represent variables, rows represent observations.
#' @param gp_train vector. Groups for each observation.
#' @param nfac numeric. The number of factors.
#' @param ... Arguments passed to \code{rstan::sampling} (e.g. iter,
#' chains).
#' @return A list.
#' \item{train}{Training data set.}
#' \item{gp_train}{A vector of groups for each observation.}
#' \item{loadings}{Aggregated loadings.}
#' \item{scores}{Aggregated factor scores.}
#' \item{noise}{Aggregated residual covariance matrix.}
#' \item{noise_org}{Aggregated residual covariance matrix, scaled by the
#' variance of the original data.}
#' \item{Sigma}{Aggregated correlation matrix.}
#' \item{stan_mod}{A list with two elements.
#' Both are objects of S4 class \code{stanfit}. First is the marginal
#' distribution estimation, second is the copula.}
train_CgFA <- function (train, gp_train, nfac = 2, ...) {
  X <- scale(train)
  X <- t(X)
  n <- ncol(X)
  m <- nrow(X)
  p <- nfac
  g <- as.integer(as.factor(gp_train))
  stan_data <- list(n  = n,
                    m  = m,
                    p  = p,
                    g  = g,
                    ng = max(g),
                    X  = X)
  blr      <- rstan::sampling(stanmodels$gFA, data = stan_data, ...)
  ext      <- extract(blr)
  lambda   <- ext$L_tri
  lambda   <- apply(lambda, MARGIN = c(2,3), FUN = median)
  factors  <- ext$Factor
  factors  <- apply(factors, MARGIN = c(2,3), FUN = median)
  psi      <- ext$Psi
  psi      <- apply(psi, MARGIN = c(2,3), FUN = median)
  params   <- list(lambda, factors, psi)
  sds      <- apply(train, MARGIN = 2, FUN = sd)
  sds      <- matrix(sds, ncol = 1)
  sds_mat  <- sds %*% t(sds)
  psi_org  <- psi
  psi      <- psi * sds_mat
  mns      <- lambda %*% factors
  for (i in 1L:nrow(X)) {
    for (j in 1L:ncol(X)) {
      u      <- pnorm(X[i,j], mean = mns[i,j], sd = psi_org[i,i])
      X[i,j] <- qnorm(u, mean = 0, sd = 1)
    }
  }
  X <- t(X)
  stan_data <- list(
    n      = nrow(X),
    m      = ncol(X),
    X      = as.matrix(X)
  )
  if (anyNA(X)) {
    sig_est <- NA
  } else {
    blr2    <- rstan::sampling(stanmodels$corr_est, data = stan_data)
    ext     <- extract(blr2)
    sig_est <- apply(ext$Sigma, c(2,3), median)
  }
  return (list(train     = train,
               gp_train  = gp_train,
               loadings  = lambda,
               scores    = factors,
               noise     = psi,
               noise_org = psi_org,
               Sigma     = sig_est,
               stan_mod  = list("gFA" = blr, "Copula" = blr2)))
}
