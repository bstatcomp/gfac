#' Train a factor analysis model
#'
#' This function trains a factor analysis model in Stan.
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
#' \item{stan_mod}{An object of S4 class \code{stanfit}.}
train_FA <- function (train, gp_train, nfac = 2, ...) {
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
  blr      <- rstan::sampling(stanmodels$FA, data = stan_data, ...)
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
  psi      <- psi * sds_mat
  return (list(train           = train,
               gp_train        = gp_train,
               loadings        = lambda,
               scores          = factors,
               noise           = psi,
               stan_mod        = blr))
}
