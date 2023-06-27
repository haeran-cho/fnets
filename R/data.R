#data

#' Simulated data from the unrestricted factor-adjusted vector autoregression model
#'
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.unrestricted(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#'
#' @format ## `unrestricted`
#' A data frame with 500 rows (observations) and 50 columns (series)
"unrestricted"


#' Simulated data from the restricted factor-adjusted vector autoregression model
#'
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.restricted(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#'
#' @format ## `restricted`
#' A data frame with 500 rows (observations) and 50 columns (series)
"restricted"
