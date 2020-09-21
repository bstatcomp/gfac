gNBGPFAC_LP_predict <- function (mod, X_test, ts_test, gp_test, ts_train, gp_train, other_pars = NULL) {
  # Generate quantities
  orig_ts_test  <- ts_test
  orig_ts_train <- ts_train

  ts_test  <- other_pars$transform_t(ts_test)
  ts_train <- other_pars$transform_t(ts_train)
  
  test_map  <- data.frame(orig = orig_ts_test, trans = ts_test)
  train_map <- data.frame(orig = orig_ts_train, trans = ts_train)
  
  lpk_kernel <- function (x1, x2, a2, r2, r3) {
    x1_mat <- matrix(data = x1, ncol = length(x2), nrow = length(x1))
    x2_mat <- matrix(data = x2, ncol = length(x2), nrow = length(x1), byrow = T)
    k      <- a2 * 
      exp(-2 * (sin(pi * abs(x1_mat - x2_mat) / other_pars$period_length))^2 / r2^2) *
      exp(-(x1_mat - x2_mat)^2/ (2 * r3^2))
    return (k)
  }
  
  get_pred <- function (x_star, x, y, sK, a1, a2, r1, r2, r3) {
    K_star_lpk    <- lpk_kernel(x_star, x, a2, r2, r3)
    K_star        <- K_star_lpk
    K_2star_lpk   <- lpk_kernel(x_star, x_star, a2, r2, r3)
    K_2star       <- K_2star_lpk
    f2_mu         <- K_star %*% sK %*% y
    cov_f2        <- K_2star - K_star %*% sK %*% t(K_star)
    tmp           <- f2_mu
    return(tmp)
  }
  
  tmp  <- mod
  tmod <- tmp$samps
  df   <- tmp$data
  
  ## extract model and get dimension
  ext   <- rstan::extract(tmod)
  nit   <- dim(ext$simp)[1]
  gmaps <- tmp$gmap
  
  ## predicted values
  f2     <- array(data = NA, dim = c(nit, df$P, length(ts_train) + length(ts_test)))
  theta2 <- array(data = NA, dim = c(nit, df$P, length(ts_train) + length(ts_test)))
  mu2    <- array(data = NA, dim = c(nit, df$M, length(ts_train) + length(ts_test)))
  
  gmaps <- gmaps[order(gmaps$ind), ] # ordered the same as in train
  for (i in 1:nit) {
    print(paste0("Iteration: ", i))
    Psi   <- as.matrix(ext$Psi[i, , ])
    a2    <- matrix(ext$a2[i, , ], nrow = mod$data$NG)
    r2    <- matrix(ext$r2[i, , ], nrow = mod$data$NG)
    r3    <- matrix(ext$r3[i, , ], nrow = mod$data$NG)
    sigma <- 0.00000001
    cs_g  <- cumsum(df$G)
    for (g in gmaps$group) {
      ts_new   <- c(ts_train[gp_train == g], ts_test[gp_test == g])
      ind_map  <- gmaps$ind[gmaps$group == g]
      if (ind_map == 1) {
        ind_train <- 1:cs_g[ind_map]
      } else {
        ind_train <- (cs_g[(ind_map - 1)] + 1):cs_g[ind_map]
      }
      ts_tr    <- df$t[ind_train]
      g_ind    <- c((gp_train == g), (gp_test == g))
      for (p in 1:df$P) {
        
        fn_train  <- ext$fn[i, p, ind_train]
        x2        <- ts_new
        x1        <- ts_tr
        y1        <- fn_train
        a2_tmp    <- a2[ind_map, p]
        r2_tmp    <- r2[ind_map, p]
        r3_tmp    <- r3[ind_map, p]
        sigma_tmp <- sigma
        K2        <- lpk_kernel(x1, x1, a2_tmp, r2_tmp, r3_tmp)
        K         <- K2
        diag(K)   <- diag(K) + sigma_tmp
        sK        <- solve(K, tol = 1e-40)
        f2_tmp <- sapply(x2, 
                         get_pred, 
                         x     = x1,
                         y     = y1,
                         sK    = sK,
                         a2    = a2_tmp,
                         r2    = r2_tmp,
                         r3    = r3_tmp)
        f2[i,p,g_ind]     <- f2_tmp
        theta2[i,p,g_ind] <- f2[i,p,g_ind]
        
      }
      base     <- ext$base
      base_tmp <- base[i,ind_map, ]
      tmp_mu   <- Psi %*% theta2[i, ,g_ind]
      tmp_o    <- sweep(tmp_mu, 1, base_tmp, FUN = "+")
      mu2[i, ,g_ind] <- exp(tmp_o)
    }
  }
  tmaps   <- rbind(train_map, test_map)
  gq_pars <- list(f2     = f2,
                  theta2 = theta2,
                  mu2    = mu2,
                  ts     = c(orig_ts_train, orig_ts_test),
                  gp     = c(as.character(gp_train), as.character(gp_test)),
                  t_maps = tmaps)
  mod <- list(trained_model = tmp, generated_quantities = gq_pars)
  
  
  # Evaluate
  ghk_cop <- function (u, L, y, mu, phi) {
    K       <- length(u)
    my_l    <- vector(mode = "numeric", length = K)
    z       <- vector(mode = "numeric", length = K)
    mbounds <- matrix(data = NA, nrow = K, ncol = 2)
    for (k in 1:K) {
      mbounds[k,1] <- qnorm(pnbinom(y[k] - 1, mu = mu[k], size = phi[k]))
      mbounds[k,2] <- qnorm(pnbinom(y[k], mu = mu[k], size = phi[k]))
      if (k == 1) {
        tmp_prod <- 0
      } else {
        L_tmp    <- L[k, 1:(k-1)]
        z_tmp    <- z[1:(k-1)]
        tmp_prod <- L_tmp %*% z_tmp
      }
      z_star_lower <- (mbounds[k,1] - tmp_prod) / L[k,k]
      z_star_upper <- (mbounds[k,2] - tmp_prod) / L[k,k]
      u_star_lower <- pnorm(z_star_lower)
      u_star_upper <- pnorm(z_star_upper)
      v       <- u_star_lower + (u_star_upper - u_star_lower) * u[k]
      z[k]    <- qnorm(v)
      my_l[k] <- u_star_upper - u_star_lower
    }
    out <- prod(my_l)
    if (is.nan(out)) out <- 0
    return (out)
  }
  
  test_df <- data.frame(ts = ts_test, gp = gp_test)
  gqs <- mod$generated_quantities
  ext <- rstan::extract(mod$trained_model$samps)
  
  L1   <- ext$L
  phis <- ext$phi
  mus1 <- gqs$mu2
  ts   <- gqs$ts
  gp   <- gqs$gp
  
  n   <- dim(mus1)[3]
  nit <- dim(mus1)[1]
  LS_mat  <- matrix(data = NA, nrow = nrow(X_test), ncol = nit)
  for (i in 1:nit) {
    print(paste0("Iteration: ", i))
    L   <- L1[i, , ]
    mus <- t(mus1[i, , ])
    pred_df    <- data.frame(ts = ts, gp = gp, pred = mus)
    pred_df$gp <- as.character(pred_df$gp)
    test_df$gp <- as.character(test_df$gp)
    spreds     <- left_join(test_df, pred_df)
    X_pred     <- spreds[ ,-which(colnames(spreds) %in% c("ts", "gp"))]
    sc <- vector(mode = "numeric", length = nrow(X_test))
    for (j in 1:length(sc)) {
      nit2 <- 100
      tmp_ll <- vector(mode = "numeric", length = nit2)
      for (k in 1:nit2) {
        u         <- runif(nrow(L), 0, 1)
        tmp_ll[k] <- ghk_cop(u, 
                             L, 
                             as.numeric(X_test[j, ]), 
                             as.numeric(X_pred[j, ]), 
                             phis[i, ])
        if (is.nan(mean(tmp_ll[k]))) browser()
        a <- ghk_cop(u, 
                     L, 
                     as.numeric(X_test[j, ]), 
                     as.numeric(X_pred[j, ]), 
                     phis[i, ])
      }
      sc[j] <- mean(tmp_ll)
    }
    LS_mat[ ,i]  <- sc
  }
  LS     <- log(apply(LS_mat, 1, mean))
  
  return (list(trained_model = tmp$trained_model,
               generated_quantities = gq_pars,
               log_scores = LS))
}