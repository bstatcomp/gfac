#' collisions data
#'
#' A data set containing NYC collision counts for different contributing
#' factors.
#'
#' \itemize{
#'     \item{Alcohol involvement}
#'     \item{Aggressive driving/Road rage}
#'     \item{Backing unsafely}
#'     \item{Brakes defective}
#'     \item{Driver inexperience}
#'     \item{Failure to yield right-of-way}
#'     \item{Other vehicular}{Vehicular problems, which need to be specified by the officer case by case.}
#'     \item{Turning improperly}
#'     \item{Passenger distraction}
#'     \item{View obstructed/limited}
#'     \item{Obstruction/Debris}
#'     \item{Traffic control disregarded}
#'     \item{Pavement slippery}
#' }
#'
#' @name collisions_data
#' @format A list with four elements.
#'   \describe{
#'     \item{fa_data}{A data frame with 1365 rows and 14 columns.}
#'     \item{group}{A vector of street names for each observation.}
#'     \item{splits}{A data frame with four columns which represent
#'     data splits. 2 stands for the train
#'     observations and 1 stands for test observations.}
#'     \item{sup_data}{A data frame with one column which represents the
#'     month.}
#'   }
#' @source \url{https://data.cityofnewyork.us/Public-Safety/Motor-Vehicle-Collisions-Crashes/h9gi-nx95}
NULL
