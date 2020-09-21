gNBGPFAC_train <- function (X, ts, gp, nfac, stan_pars, init_list = NULL, other_pars = NULL, ...) {
  tmp <- bind_cols(as.data.frame(X), ts = ts, group = gp) %>%
    arrange(group)
  # browser()
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
  t_ar       <- other_pars$transform_t(tmp$ts)
  
  # get prior parameters
  prior_pars <- other_pars$prior_pars
  
  # seed
  if (is.null(other_pars$tseed)) {
    tseed <- 1
  } else {
    tseed <- other_pars$tseed
  }
  if (is.null(other_pars$adelta)) {
    tmp_adelta <- 0.8
  } else {
    tmp_adelta <- other_pars$adelta
  }
  
  
  stan_data  <- list(
    N      = n,
    NG     = length(g),
    M      = m,
    P      = p,
    t      = t_ar,
    G      = as.array(g),
    X      = t(X_ar),
    
    p_dist          = prior_pars$p_dist,
    p_val           = prior_pars$p_val,
    period_length   = other_pars$period_length
  )
  
  stan_mod <- get_stan_compiled("gNBGPFAC_LP")
  
  if (!is.null(init_list)) {
    il <- list()
    for (i in 1:stan_pars$nchain) {
      il[[length(il) + 1]] <- init_list
    }
    samps    <- rstan::sampling(stan_mod, stan_data, seed = tseed,
                                chains = stan_pars$nchain,
                                iter = stan_pars$nit,
                                init = il,
                                control = list(adapt_delta = tmp_adelta),
                                ...)
  } else {
    samps    <- rstan::sampling(stan_mod, stan_data, seed = tseed,
                                chains = stan_pars$nchain,
                                iter = stan_pars$nit,
                                control = list(adapt_delta = tmp_adelta),
                                ...)
  }
  
  out_list <- list(data   = stan_data,
                   samps  = samps,
                   gmap   = gmap)
  return(out_list)
}