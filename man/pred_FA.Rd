% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred_FA.R
\name{pred_FA}
\alias{pred_FA}
\title{Predict from FA model}
\usage{
pred_FA(tr_fa, gp_test, test)
}
\arguments{
\item{tr_fa}{An object resulting from \code{train_FA}.}

\item{gp_test}{vector. Groups for each observation.}

\item{test}{data frame. Test data set for evaluation.}
}
\value{
A list.
  \item{predictions}{A numeric vector of predictions for all observations
  in the test set.}
  \item{log_lik}{A numeric vector of (log-)likelihoods for all observations
  in the test set.}
}
\description{
This function predicts based on the estimation of the FA
model \code{train_FA}.
}
\seealso{
\code{\link{train_FA}}
}
