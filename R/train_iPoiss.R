#' Train an independent Poisson model
#'
#' This function trains an independent Poisson model for multivariate
#' count estimation.
#'
#' @export
#' @param train data frame. The training data set, consisting of counts.
#'   Columns represent variables, rows represent observations.
#' @param gp_train vector. Groups for each observation.
#' @return A list.
#'   \item{params}{A matrix of Poisson parameters for each group-variable
#'   combination. Rows represent groups, columns represent variables.}
#'   \item{train}{Train data set.}
#'   \item{gp_train}{Groups of the train data set.}
#' @seealso \code{\link{pred_iPoiss}}
train_iPoiss <- function (train, gp_train) {
  params <- matrix(data = NA, nrow = length(unique(gp_train)), ncol=ncol(train))
  ind_g  <- 1
  for (g in unique(gp_train)) {
    for (j in 1L:ncol(train)) {
      g_data          <- train[gp_train == g,j]
      gamma_0         <- c(0.01, 0.01)
      alpha_p         <- sum(g_data) + gamma_0[1]
      beta_p          <- length(g_data) + gamma_0[2]
      lamb_mean       <- alpha_p / beta_p
      params[ind_g,j] <- lamb_mean
    }
    ind_g <- ind_g + 1
  }
  return (list(params   = params,
               train    = train,
               gp_train = gp_train))
}
