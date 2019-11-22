#' Predict from FA model
#'
#' This function predicts based on the estimation of the FA
#' model \code{train_FA}.
#'
#' @export
#' @param tr_fa An object resulting from \code{train_FA}.
#' @param gp_test vector. Groups for each observation.
#' @param test data frame. Test data set for evaluation.
#' @param logl logical. Return log-likelihood? Defaults to TRUE.
#' @param L_intg logical. Should the likelihood be integrated across
#' suitable borders? Defaults to FALSE, which uses an approximation with
#' MVN density.
#' Special care needs to be taken when integrating, as the likelihood is
#' likely to be 0 in higher dimensions due to machine precision.
#' @return A list.
#'   \item{predictions}{A numeric vector of predictions for all observations
#'   in the test set.}
#'   \item{log_lik}{A numeric vector of (log-)likelihoods for all observations
#'   in the test set.}
#' @seealso \code{\link{train_FA}}
pred_FA <- function (tr_fa, gp_test, test, logl = T, L_intg = F) {
  predictions <- matrix(data = NA, nrow=length(gp_test), ncol=ncol(tr_fa$train))
  loadings    <- matrix(tr_fa$loadings, ncol = ncol(tr_fa$loadings))
  scores      <- t(as.matrix(tr_fa$scores))
  gp_score    <- data.frame(group = tr_fa$gp_train, scores = scores)
  gp_score <- aggregate(gp_score[ ,-1], by = list(gp_score$group), FUN = mean)
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
  means       <- apply(tr_fa$train, MARGIN = 2, FUN = mean)
  sds         <- apply(tr_fa$train, MARGIN = 2, FUN = sd)
  predictions <- sweep(predictions, 2, sds, '*')
  predictions <- sweep(predictions, 2, means, '+')
  predictions <- round(predictions)
  log_lik     <- likelihood_cl(test, predictions, tr_fa, logl,
                               integrate = L_intg)
  return (list(predictions = predictions,
               log_lik     = log_lik))
}
