#' Train a copula PFA model
#'
#' This function trains a copula factor analysis model with
#' Poisson marginals in Stan.
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
#' \item{Sigma}{Aggregated correlation matrix.}
#' \item{stan_mod}{A list with two elements.
#' Both are objects of S4 class \code{stanfit}. First is the marginal
#' distribution estimation, second is the copula.}
train_CgPFA <- function (train, gp_train, nfac = 2, ...) {
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
  blr      <- rstan::sampling(stanmodels$gPFA, data = stan_data, ...)
  ext      <- extract(blr)
  lambda   <- apply(ext$Psi, MARGIN = c(2,3), FUN = median)
  factors  <- apply(ext$Theta, MARGIN = c(2,3), FUN = median)

  gp_score    <- data.frame(group = unique(gp_train), scores = t(factors))
  gp_uniq     <- unique(gp_score)
  colnames(gp_uniq)[1] <- "group"
  pr_factors  <- join(data.frame(group = gp_train), gp_uniq, by = "group")
  thetas      <- pr_factors[ ,-1]
  thetas      <- t(lambda %*% t(thetas))
  X_lat       <- thetas
  X           <- t(X)
  for (i in 1L:nrow(X_lat)) {
    for (j in 1L:ncol(X_lat)) {
      u <- runif(1, ppois(X[i,j] - 1, thetas[i,j]), ppois(X[i,j], thetas[i,j]))
      X_lat[i,j] <- qnorm(u, mean = 0, sd = 1)
    }
  }
  stan_data <- list(
    n      = nrow(X_lat),
    m      = ncol(X_lat),
    X      = as.matrix(X_lat)
  )
  blr2    <- rstan::sampling(stanmodels$corr_est,
                             data   = stan_data,
                             ...)
  ext     <- extract(blr2)
  sig_est <- apply(ext$Sigma, c(2,3), median)
  return (list(train    = train,
               gp_train = gp_train,
               loadings = lambda,
               scores   = factors,
               Sigma    = sig_est,
               stan_mod = list("PFA" = blr, "Copula" = blr2)))
}
