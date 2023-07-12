#data

#' Simulated data from the unrestricted factor-adjusted vector autoregression model
#'
#'
#' \code{set.seed(123)} \cr
#' \code{n <- 500} \cr
#' \code{p <- 50} \cr
#' \code{common <- sim.unrestricted(n, p)} \cr
#' \code{idio <- sim.var(n, p)} \cr
#' \code{x <- common$data + idio$data}
#'
#' @format ## `data.unrestricted`
#' A ts object with 500 rows (observations) and 50 columns (series)
"data.unrestricted"


#' Simulated data from the restricted factor-adjusted vector autoregression model
#'
#' \code{set.seed(123)} \cr
#' \code{n <- 500} \cr
#' \code{p <- 50} \cr
#' \code{common <- sim.restricted(n, p)} \cr
#' \code{idio <- sim.var(n, p)} \cr
#' \code{x <- common$data + idio$data}
#'
#' @format ## `data.restricted`
#' A ts object with 500 rows (observations) and 50 columns (series)
"data.restricted"
