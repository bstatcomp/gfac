#' Predict from CgPFA model
#'
#' This function predicts based on the estimation of the CgPFA
#' model \code{train_CgPFA}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_CgPFA}.
#' @param gp_test vector. Groups for each observation.
#' @param test data frame. Test data set for evaluation.
#' @param logl logical. Return log-likelihood? Defaults to TRUE.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{likelihoods}{A numeric vector of (log-)likelihoods for all
#'   observations in the test set.}
#' @seealso \code{\link{train_CgPFA}}
pred_CgPFA <- function (tr_fa, gp_test, test, logl = T) {
  predictions <- matrix(data = NA, nrow=length(gp_test), ncol=ncol(tr_fa$train))
  loadings    <- matrix(tr_fa$loadings, ncol = ncol(tr_fa$loadings))
  scores      <- t(as.matrix(tr_fa$scores))
  gp_score    <- data.frame(group = unique(tr_fa$gp_train), scores = scores)
  gp_uniq     <- unique(gp_score)
  colnames(gp_uniq)[1] <- "group"
  pr_factors  <- join(data.frame(group = gp_test), gp_uniq, by = "group")
  if (anyNA(pr_factors)) {
    ind                 <- apply(is.na(pr_factors), 1, sum)
    ind[ind != 0]       <- 1
    ind                 <- as.logical(ind)
    pr_factors[ind, -1] <- matrix(rep(apply(as.matrix(gp_score[ ,-1]), 2, mean),
                                      times = sum(ind)),
                                  nrow = sum(ind), byrow = T)
  }
  predictions <- as.matrix(pr_factors[ ,-1]) %*% t(loadings)
  tmp_ts      <- as.matrix(test)
  log_lik     <- vector(mode = "numeric", length = nrow(test))
  for (i in 1L:length(gp_test)) {
    log_lik[i] <- L_CgPFA_approx(test[i,],
                                predictions[i,],
                                tr_fa$Sigma,
                                logl = logl)
  }

  return (list(predictions = predictions,
               likelihoods = log_lik))
}
