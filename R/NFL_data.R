#' NFL data
#'
#' A data set containing the NFL statistics for several seasons.
#'
#' \itemize{
#'     \item{APAS}{attempted passes}
#'     \item{CPAS}{completed passes}
#'     \item{PASY}{yards gained by passing}
#'     \item{PASTD}{touchdowns by passing}
#'     \item{INT}{intercepted passes}
#'     \item{ARSH}{attempted rushes}
#'     \item{RSHY}{yards gained by rushing}
#'     \item{RSHTD}{touchdowns by rushing}
#'     \item{TAC}{tackles}
#'     \item{SAC}{sacks}
#' }
#' @name NFL_data
#' @format A list with four elements.
#'   \describe{
#'     \item{fa_data}{A data frame with 1536 rows and 10 columns.}
#'     \item{group}{A vector of team names for each observation.}
#'     \item{splits}{A data frame with three columns which represent
#'     time-respecting data splits. 2 stands for the train
#'     observations and 1 stands for test observations.}
#'     \item{sup_data}{A data frame with one column which represents the date
#'     of the game.}
#'   }
#' @references Horowitz, M., Yurko, R., Ventura, S. (2018) nflscrapR: Compiling
#' the NFL play-by-play API for easy use in R.
#' (\href{https://github.com/maksimhorowitz/nflscrapR})
NULL
