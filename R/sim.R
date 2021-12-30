#' @title Simulate data from a dynamic factor model
#' @description Simulate the common component following a dynamic factor model that does not admit a static representation;
#' see the model (C1) in the reference.
#' @param n sample size
#' @param p dimension
#' @param q number of dynamic factors
#' @return a list containing
#' \itemize{
#' \item{data}{ generated series}
#' \item{q}{ number of factors}
#' }
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @examples
#' common <- sim.factor.M1(500, 50)
#' @export
sim.factor.M1 <- function(n, p, q = 2){
  
  trunc.lags <- min(20, round(n/log(n)))
  chi <- matrix(0, p, n)
  uu <- matrix(rnorm((n + trunc.lags) * q), ncol = q)
  a <- matrix(runif(p * q, -1, 1), ncol = q)
  alpha <- matrix(runif(p * q, -.8, .8), ncol = q)
  for(ii in 1:p){
    for(jj in 1:q){
      coeffs <- alpha[ii, jj] * c(1, as.numeric(var.to.vma(as.matrix(a[ii, jj]), trunc.lags)))
      for(tt in 1:n) chi[ii, tt] <- chi[ii, tt] + coeffs %*% uu[(tt + trunc.lags + 1):tt, jj]
    }
  }
  return(list(data = chi, q = q))
  
}

#' @title Simulate data from a static factor model
#' @description Simulate the common component following a dynamic factor model that admits a static representation;
#' see the model (C2) in the reference.
#' @param n sample size
#' @param p dimension
#' @param q number of dynamic factors; number of static factors is given by \code{2 * q}
#' @return a list containing
#' \itemize{
#' \item{data}{ generated series}
#' \item{q}{ number of factors}
#' \item{r}{ number of static factors}
#' }
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @examples
#' common <- sim.factor.M2(500, 50)
#' @export
sim.factor.M2 <- function(n, p, q = 2){
  lags <- 1
  r <- q * (lags + 1)
  burnin <- 100
  u2 <- matrix(rnorm(q * (n + burnin)), nrow = q)
  D0 <- matrix(runif(q^2, 0, .3), nrow = q) 
  diag(D0) <- runif(q, .5, .8)
  D <- 0.7 * D0/norm(D0, type = '2')

  f <- matrix(0, nrow = q, ncol = n + burnin)
  f[, 1] <-  u2[, 1]
  for(tt in 2:(n + burnin)) f[, tt] <- D %*% f[, tt - 1] +  u2[, tt]
  f <- f[, -(1:(burnin - lags))]
  
  loadings <- matrix(rnorm(p * r, 0, 1), nrow = p) 
  chi <- matrix(0, p, n)
  for(ii in 0:lags) chi <- chi + loadings[, ii*q + 1:q] %*% f[, 1:n + lags - ii]
  return(list(data = chi, q = q, r = r))
  
}

#' @title Simulate data from a VAR(1) model
#' @description Simulate the idiosyncratic component following a VAR(1) model; see the reference for further details.
#' @param n sample size
#' @param p dimension
#' @param Gamma innovation covariance matrix
#' @return a list containing
#' \itemize{
#' \item{data}{ generated series}
#' \item{A}{ transition matrix}
#' \item{Gamma}{ innovation covariance matrix}
#' }
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @examples
#' idio <- sim.idio(500, 50)
#' @export
sim.idio <- function(n, p, Gamma = diag(1, p)){
  burnin <-100
  prob <- 1/p
  xi <- t(MASS::mvrnorm(n + burnin, mu = rep(0, p), Sigma = Gamma))
  A <- matrix(0, p, p)
  index <- sample(c(0, 1), p^2, TRUE, prob = c(1 - prob, prob))
  A[which(index == 1)] <- .275
  A <- A / norm(A, "2")
  
  for(tt in 2:(n + burnin)) xi[, tt] <- xi[, tt] + A %*% xi[, tt - 1]
  xi <- xi[, -(1:burnin)]

  return(list(data = xi, A = A, Gamma = Gamma))
  
}

