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
#' @format ## `data.unrestricted`
#' A data frame with 500 rows (observations) and 50 columns (series)
"data.unrestricted"


#' Simulated data from the restricted factor-adjusted vector autoregression model
#'
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.restricted(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#'
#' @format ## `data.restricted`
#' A data frame with 500 rows (observations) and 50 columns (series)
"data.restricted"
