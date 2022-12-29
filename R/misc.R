#' @title Bartlett weights
#' @keywords internal
Bartlett.weights <- function(x)
  1 - abs(x)

#' @keywords internal
make.symmetric <- function(DD, symmetric) {
  symmetric <- match.arg(symmetric, c("min", "max", "avg", "none"))
  if(symmetric != "none") {
    if(symmetric == "min") {
      DD[abs(DD) > t(abs(DD))] <- 0
      DD <- DD + t(DD)
      diag(DD) <- diag(DD) / 2
    }
    if(symmetric == "max") {
      DD[abs(DD) < t(abs(DD))] <- 0
      DD <- DD + t(DD)
      diag(DD) <- diag(DD) / 2
    }
    if(symmetric == "avg")
      DD <- (DD + t(DD)) / 2
  }
  DD
}


#' @keywords internal
var.to.vma <- function(A, trunc.lags) {
  d <- dim(A)[1]
  s <- dim(A)[2]
  l <- s / d
  B <- array(0, dim = c(d, d, trunc.lags + 1))
  B[, , 1] <- diag(1, d)
  for (ii in 1:trunc.lags) {
    for (jj in 1:min(ii, l))
      B[, , ii + 1] <-
        B[, , ii + 1] + B[, , ii - jj + 1] %*% A[, (jj - 1) * d + 1:d]
  }
  B
}

#' @keywords internal
check.list.arg <- function(arg) {
  arg.name <- deparse(substitute(arg))
  if(arg.name == "var.args")
    default <-
      list(
        n.iter = 100,
        tol = 0,
        n.cores = min(parallel::detectCores() - 1, 3)
      )
  if(arg.name == "tuning.args")
    default <-
      list(
        tuning = c("cv", "bic"),
        n.folds = 1,
        penalty = NULL,
        path.length = 10,
        do.plot = FALSE
      )
  if(arg.name == "common.args")
    default <-
      list(
        var.order = NULL,
        max.var.order = NULL,
        trunc.lags = 20,
        n.perm = 10
      )
  out <- check.list.match(arg, default)
  return(out)
}

#' @keywords internal
check.list.match <- function(arg, default) {
  for (ii in names(default)) {
    if(is.null(arg[[ii]]))
      arg[[ii]] <- default[[ii]]
  }
  return(arg)
}
