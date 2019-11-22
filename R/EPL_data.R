#' EPL data
#'
#' A data set containing the EPL statistics for several seasons.
#'
#' #' \itemize{
#'     \item{FTHG}{home goals}
#'     \item{FTAG}{away goals}
#'     \item{HS}{home shots}
#'     \item{AS}{away shots}
#'     \item{HST}{home shots on target}
#'     \item{AST}{away shots on target}
#'     \item{HF}{home fouls}
#'     \item{AF}{away fouls}
#'     \item{HC}{home corners}
#'     \item{AC}{away corners}
#'     \item{HY}{home yellow cards}
#'     \item{AY}{away yellow cards}
#'     \item{HR}{home red cards}
#'     \item{AR}{away red cards}
#' }
#'
#' @name EPL_data
#' @format A list with four elements.
#'   \describe{
#'     \item{fa_data}{A data frame with 1900 rows and 14 columns.}
#'     \item{group}{A vector of team names for each observation.}
#'     \item{splits}{A data frame with three columns which represent
#'     time-respecting data splits. 2 stands for the train
#'     observations and 1 stands for test observations.}
#'     \item{sup_data}{A data frame with two columns which represent the date
#'     of the game and the away team.}
#'   }
#' @source \url{http://www.football-data.co.uk/data.php}
NULL
