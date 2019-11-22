#' Predict from CgNBFA model
#'
#' This function predicts based on the estimation of the CgNBFA
#' model \code{train_CgNBFA}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_CgNBFA}.
#' @param gp_test vector. Groups for each observation.
#' @param test data frame. Test data set for evaluation.
#' @param logl logical. Return log-likelihood? Defaults to TRUE.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{likelihoods}{A numeric vector of (log) likelihoods for all
#'   observations in the test set.}
#' @seealso \code{\link{train_CgNBFA}}
pred_CgNBFA <- function (tr_fa, gp_test, test, logl = TRUE) {
  predictions <- matrix(data = NA, nrow=length(gp_test), ncol=ncol(tr_fa$train))
  mus         <- t(as.matrix(tr_fa$mu))
  size        <- t(as.matrix(tr_fa$phi))
  gpt         <- unique(tr_fa$gp_train)
  gpt         <- sort(as.character(gpt))
  gp_mu       <- data.frame(group = gpt, mus = mus)
  gp_size     <- data.frame(group = gpt, size = size)
  colnames(gp_mu)[1]   <- "group"
  colnames(gp_size)[1] <- "group"
  pr_mu         <- join(data.frame(group = gp_test), gp_mu, by = "group")
  pr_size     <- join(data.frame(group = gp_test), gp_size, by = "group")
  if (anyNA(pr_mu)) {
    ind                 <- apply(is.na(pr_mu), 1, sum)
    ind[ind != 0]       <- 1
    ind                 <- as.logical(ind)
    ap_df               <- apply(gp_mu[ ,-1], 1, mean)
    md                  <- median(ap_df, na.rm = T)
    wh_md               <- which.min(abs(ap_df - md))
    tmp_mu              <- pr_mu[wh_md, -1]
    tmp_size            <- pr_size[wh_md, -1]
    pr_mu[ind, -1]      <- tmp_mu
    pr_size[ind, -1]    <- tmp_size
  }
  log_lik <- vector(mode = "numeric", length = nrow(test))
  for (i in 1L:length(gp_test)) {
    sizes1   <- as.numeric(pr_size[i, -1])
    mus1     <- as.numeric(pr_mu[i, -1])
    predictions[i, ] <- mus1
    log_lik[i]       <- L_CgNBFA_approx(test[i,],
                                      sizes1,
                                      mus1,
                                      tr_fa$Sigma,
                                      logl = logl)
  }
  return (list(predictions = predictions,
               likelihoods = log_lik))
}
