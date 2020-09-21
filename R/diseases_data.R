#' Diseases data
#'
#' A data set containing monthly counts of 7 major diseases in Slovenia from 2000-2017.
#'
#' \itemize{
#'     \item{A02}{salmonella}
#'     \item{A04}{e.coli infection}
#'     \item{A08}{rotaviral enteritis}
#'     \item{A09}{astroenteritis and colitis of unspecified origin}
#'     \item{A38}{scarlet fever}
#'     \item{A46}{erysipelas}
#'     \item{A69}{lyme disease}
#' }
#'
#' @name diseases_data
#' @format A list with four elements.
#'   \describe{
#'     \item{fa_data}{A data frame with 216 rows and 7 columns.}
#'     \item{ts}{A vector of length 216 of dates.}
#'     \item{gp}{A vector of length 216. For this data set we only have one group.}
#'     \item{splits}{A data frame with 216 rows and 12 columns. Indexing for fone-year-ahead forecasting.}
#'   }
#' @source \url{https://www.nijz.si}
NULL
