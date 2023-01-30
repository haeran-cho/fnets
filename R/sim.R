#' @title Simulate data from an unrestricted factor model
#' @description Simulate the common component following an unrestricted factor model that does not admit a restricted representation;
#' see the model (C1) in the reference.
#' @param n sample size
#' @param p dimension
#' @param q number of unrestricted factors
#' @param heavy if \code{heavy = FALSE}, common shocks are generated from \code{rnorm} whereas if \code{heavy = TRUE}, from \code{rt} with \code{df = 5} and then scaled by \code{sqrt(3 / 5)}
#' @return a list containing
#' \item{data}{ generated series}
#' \item{q}{ number of factors}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. arXiv preprint arXiv:2301.11675.
#' @examples
#' common <- sim.unrestricted(500, 50)
#' @importFrom stats rnorm runif rt
#' @export
sim.unrestricted <- function(n, p, q = 2, heavy = FALSE) {
  trunc.lags <- min(20, round(n / log(n)))
  chi <- matrix(0, p, n)
  if(!heavy) {
    uu <- matrix(rnorm((n + trunc.lags) * q), ncol = q)
  } else {
    uu <-
      matrix(rt((n + trunc.lags) * q, df = 5), ncol = q) * sqrt(3 / 5)
  }
  a <- matrix(runif(p * q,-1, 1), ncol = q)
  alpha <- matrix(runif(p * q,-.8, .8), ncol = q)
  for (ii in 1:p) {
    for (jj in 1:q) {
      coeffs <-
        alpha[ii, jj] * as.numeric(var.to.vma(as.matrix(a[ii, jj]), trunc.lags))
      for (tt in 1:n)
        chi[ii, tt] <-
          chi[ii, tt] + coeffs %*% uu[(tt + trunc.lags):tt, jj]
    }
  }
  return(list(data = chi, q = q))
}

#' @title Simulate data from a restricted factor model
#' @description Simulate the common component following an unrestricted factor model that admits a restricted representation;
#' see the model (C2) in the reference.
#' @param n sample size
#' @param p dimension
#' @param q number of unrestricted factors; number of restricted factors is given by \code{2 * q}
#' @param heavy if \code{heavy = FALSE}, common shocks are generated from \code{rnorm} whereas if \code{heavy = TRUE}, from \code{rt} with \code{df = 5} and then scaled by \code{sqrt(3 / 5)}
#' @return a list containing
#' \item{data}{ generated series}
#' \item{q}{ number of factors}
#' \item{r}{ number of restricted factors}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling
#' @examples
#' common <- sim.restricted(500, 50)
#' @importFrom stats rnorm runif rt
#' @export
sim.restricted <- function(n, p, q = 2, heavy = FALSE) {
  lags <- 1
  r <- q * (lags + 1)
  burnin <- 100
  if(!heavy) {
    uu <- matrix(rnorm((n + burnin) * q), nrow = q)
  } else {
    uu <- matrix(rt((n + burnin) * q, df = 5), nrow = q) * sqrt(3 / 5)
  }
  D0 <- matrix(runif(q^2, 0, .3), nrow = q)
  diag(D0) <- runif(q, .5, .8)
  D <- 0.7 * D0 / norm(D0, type = "2")

  f <- matrix(0, nrow = q, ncol = n + burnin)
  f[, 1] <- uu[, 1]
  for (tt in 2:(n + burnin))
    f[, tt] <- D %*% f[, tt - 1] + uu[, tt]
  f <- f[,-(1:(burnin - lags))]

  loadings <- matrix(rnorm(p * r, 0, 1), nrow = p)
  chi <- matrix(0, p, n)
  for (ii in 0:lags)
    chi <- chi + loadings[, ii * q + 1:q] %*% f[, 1:n + lags - ii]
  return(list(data = chi, q = q, r = r))
}

#' @title Simulate a VAR(1) process
#' @description Simulate a VAR(1) process; see the reference for the generation of the transition matrix.
#' @param n sample size
#' @param p dimension
#' @param Gamma innovation covariance matrix; ignored if \code{heavy = TRUE}
#' @param heavy if \code{heavy = FALSE}, common shocks are generated from \code{rnorm} whereas if \code{heavy = TRUE}, from \code{rt} with \code{df = 5} and then scaled by \code{sqrt(3 / 5)}
#' @return a list containing
#' \item{data}{ generated series}
#' \item{A}{ transition matrix}
#' \item{Gamma}{ innovation covariance matrix}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling
#' @examples
#' idio <- sim.var(500, 50)
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rt
#' @export
sim.var <- function(n,
                    p,
                    Gamma = diag(1, p),
                    heavy = FALSE) {
  burnin <- 100
  prob <- 1 / p

  if(!heavy) {
    if(identical(Gamma, diag(1, p))) {
      xi <- matrix(rnorm((n + burnin) * p), nrow = p)
    } else {
      xi <- t(MASS::mvrnorm(n + burnin, mu = rep(0, p), Sigma = Gamma))
    }
  } else {
    xi <- matrix(rt((n + burnin) * p, df = 5), nrow = p) * sqrt(3 / 5)
  }

  A <- matrix(0, p, p)
  index <- sample(c(0, 1), p^2, TRUE, prob = c(1 - prob, prob))
  A[which(index == 1)] <- .275
  A <- A / norm(A, "2")

  for (tt in 2:(n + burnin))
    xi[, tt] <- xi[, tt] + A %*% xi[, tt - 1]
  xi <- xi[,-(1:burnin)]

  return(list(data = xi, A = A, Gamma = Gamma))
}
