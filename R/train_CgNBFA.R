#' Train a copula NBFA2 model
#'
#' This function trains a copula factor analysis model with
#' negative binomial marginals in Stan. Each variable is assigned its own
#' scaling parameter.
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
#' \item{mu}{Aggregated mean parameters.}
#' \item{phi}{Aggregated size parameters.}
#' \item{Sigma}{Aggregated correlation matrix.}
#' \item{stan_mod}{A list with two elements.
#' Both are objects of S4 class \code{stanfit}. First is the marginal
#' distribution estimation, second is the copula.}
train_CgNBFA <- function (train, gp_train, nfac = 2, ...) {
  X <- t(train)
  n <- ncol(X)
  m <- nrow(X)
  p <- nfac
  g <- as.integer(as.factor(as.character(gp_train)))
  stan_data <- list(n  = n,
                    m  = m,
                    p  = p,
                    g  = g,
                    ng = max(g),
                    X  = X)
  blr      <- rstan::sampling(stanmodels$gNBFA, data = stan_data, ...)
  ext      <- extract(blr)
  lambda   <- apply(ext$Psi, MARGIN = c(2,3), FUN = median)
  factors  <- apply(ext$Theta, MARGIN = c(2,3), FUN = median)
  mu       <- apply(ext$mu, MARGIN = c(2,3), FUN = median)
  size     <- apply(ext$phi, MARGIN = c(2), FUN = median)
  cop_size <- size
  cop_mu   <- mu
  X        <- train
  for (i in 1L:nrow(X)) {
    for (j in 1L:ncol(X)) {
      u <- runif(1,
                 pnbinom(train[i,j]-1,
                         size = cop_size[j],
                         mu   = cop_mu[j, g[i]]),
                 pnbinom(train[i,j],
                         size = cop_size[j],
                         mu   = cop_mu[j, g[i]]))
      X[i,j] <- qnorm(u, mean = 0, sd = 1)
    }
  }
  stan_data <- list(
    n      = nrow(X),
    m      = ncol(X),
    X      = as.matrix(X)
  )
  blr2    <- rstan::sampling(stanmodels$corr_est, data = stan_data, ...)
  ext     <- extract(blr2)
  sig_est <- apply(ext$Sigma, c(2,3), median)
  return (list(train    = train,
               gp_train = gp_train,
               loadings = lambda,
               scores   = factors,
               mu       = mu,
               phi      = size,
               Sigma    = sig_est,
               stan_mod = list("PFA" = blr, "Copula" = blr2)))
}
