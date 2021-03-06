#' Train a Gaussian process factor analysis model with full covariance
#'
#' This function trains a Gaussian process factor analysis model with full covariance and a locally periodic kernel in Stan. Uses a normal likelihood and standardized data.
#'
#' @export
#' @param X data frame. The training data set, consisting of counts. Columns represent variables, rows represent observations.
#' @param ts vector. Time-points of observations in X.
#' @param gp vector. Groups for each observation.
#' @param nfac numeric. Number of factors.
#' @param nit numeric. Number of iterations.
#' @param nchain numeric. Number of chains.
#' @param period_length numeric. Length of the period.
#' @param prior_r2 vector. Parameters (2) for the inverse-gamma prior distribution for the lengthscale of the squared exponential
#' part of the locally periodic kernel.
#' @param prior_r3 vector. Parameters (2) for the inverse-gamma prior distribution for the lengthscale of the periodic
#' part of the locally periodic kernel.
#' @param transform_t function. Function for the transformation of the time-series (for example if we want it on a specific interval).
#' @param tseed numeric. Seed for sampling Defaults to 1.
#' @return A list.
#' \item{data}{Data that goes into sampling.}
#' \item{samps}{An object of S4 class \code{stanfit}. Samples.}
#' \item{gmap}{Mapping of groups to indexes.}
#' \item{orig_X}{Original X.}
train_gMVNGPFA_LP <- function (X, ts, gp, nfac, nit, nchain, period_length,
                            prior_r2, prior_r3,
                            transform_t, tseed = NULL,...) {
  tmp <- bind_cols(as.data.frame(X), ts = ts, group = gp) %>%
    arrange(group)

  n <- nrow(X)
  m <- ncol(X)
  p <- nfac

  X_ar <- dplyr::select(tmp, - c(ts, group)) %>%
    as.matrix()

  g <- count(tmp, group) %>%
    dplyr::select(n) %>%
    unlist() %>%
    as.vector()
  gmap <- count(tmp, group) %>%
    bind_cols(ind = 1:length(g))

  # rescale time-series
  t_ar       <- transform_t(tmp$ts)


  # standardize
  orig_X <- X_ar
  X_ar   <- scale(X_ar)


  stan_data  <- list(
    N      = n,
    NG     = length(g),
    M      = m,
    P      = p,
    t      = t_ar,
    G      = as.array(g),
    X      = t(X_ar),

    p_dist          = c(1, 1, 1),
    p_val           = c(prior_r2, prior_r3, prior_phi),
    period_length   = period_length
  )
  # stan_mod <- get_stan_compiled("gGPFA_LP")
  if (!is.null(tseed)) {
    samps    <- rstan::sampling(stanmodels$gMVNGPFA_LP,
                                stan_data,
                                seed = tseed,
                                chains = nchain,
                                iter = nit,
                                ...)
  } else {
    samps    <- rstan::sampling(stanmodels$gMVNGPFA_LP,
                                stan_data,
                                seed = 1,
                                chains = nchain,
                                iter = nit,
                                ...)
  }

  out_list <- list(data   = stan_data,
                   samps  = samps,
                   gmap   = gmap,
                   orig_X = orig_X)
  return(out_list)
}
