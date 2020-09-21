#' toy A data
#'
#' A generated data set with low overdispersed counts (Poisson), temporal latent 3-factor structure, and residual covariance not due to the latent structure.
#'
#' @name toy_A
#' @format A list with four elements.
#'   \describe{
#'     \item{df}{A data frame with 400 rows and 10 columns.}
#'     \item{ts}{A vector of length 400 of time-points corresponding to each observation.}
#'     \item{gp}{A vector of length 400 of groups for each observation.}
#'     \item{splits}{A data frame with 400 rows and two columns. Indexing of a 2-fold CV split of the data.}
#'   }
NULL
