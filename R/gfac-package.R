#' The 'gfac' package.
#'
#' @description The package contains several factor analysis methods, with
#' emphasis to count data. The main methods rely on copulas to estimate the
#' residual covariance, not explained by the latent structure. Additionally, the
#' package contains several data sets for empirical evaluation. Methods use RStan for inference.
#'
#' @docType package
#' @name gfac-package
#' @aliases gfac
#' @useDynLib gfac, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#'
NULL
