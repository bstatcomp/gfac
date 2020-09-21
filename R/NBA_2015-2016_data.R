#' NBA 2015-2016 data
#'
#' A data set containing the NBA statistics for 2015-2016 season.
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
#'     \item{fa_data}{A data frame with 2632 rows and 15 columns.}
#'     \item{ts}{A vector of length 2632 of time-points corresponding to each observation.}
#'     \item{gp}{A vector of length 2632 of groups for each observation.}
#'     \item{splits}{A data frame with 2632 rows and 4 columns. Indexing of a 4-fold CV split of the data.}
#'   }
#' @source \url{https://www.basketball-reference.com/}
NULL
