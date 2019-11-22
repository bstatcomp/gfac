#' NBA data
#'
#' A data set containing the NBA statistics for several seasons.
#'
#' \itemize{
#'     \item{M2FG}{2-point field goals made}
#'     \item{A2FG}{2-point field goals attempted}
#'     \item{M3FG}{3-point field goals made}
#'     \item{A3FG}{3-point field goals attempted}
#'     \item{MFT}{free-throws made}
#'     \item{AFT}{free-throws attempted}
#'     \item{O}{offensive rebounds}
#'     \item{D}{defensive rebounds}
#'     \item{As}{assists}
#'     \item{St}{steals}
#'     \item{To}{lost balls}
#'     \item{Fv}{blocks}
#'     \item{Ag}{blocks against}
#'     \item{Cm}{personal fouls}
#'     \item{Rv}{personal fouls against}
#' }
#'
#' @name NBA_data
#' @format A list with four elements.
#'   \describe{
#'     \item{fa_data}{A data frame with 43814 rows and 15 columns.}
#'     \item{group}{A vector of team names for each observation.}
#'     \item{splits}{A data frame with three columns which represent
#'     time-respecting data splits. 2 stands for the train
#'     observations and 1 stands for test observations.}
#'     \item{sup_data}{A data frame with three columns which represent the date
#'     of the game, season, and where the team is playing.}
#'   }
#' @source \url{https://www.basketball-reference.com/}
NULL
