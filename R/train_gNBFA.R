#' Train a negative binomial factor analysis model
#'
#' This function trains a negative binomial factor analysis model in Stan.
#' Each variable has its own scaling parameter.
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
#' \item{probs}{Aggregated probability parameters.}
#' \item{phi}{Aggregated size parameters.}
#' \item{stan_mod}{An object of S4 class \code{stanfit}.}
train_gNBFA <- function (train, gp_train, nfac = 2, ...) {
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
  return (list(train    = train,
               gp_train = gp_train,
               loadings = lambda,
               scores   = factors,
               mu       = mu,
               phi      = size,
               stan_mod = blr))
}
