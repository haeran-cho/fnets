# simulate data

#' Simulate data from a dynamic factor model (Model 1)
#'
#' @param n sample size
#' @param p number of series
#' @param q dynamic dimension (default 2)
#' @param r factor number (default 4)
#' @param do.scale scale the output (default true)
#' @param K  transition matrix, dimension r by q (default null)
#' @param loadings  loading matrix, dimension p by r (default null)
#'
#' @return  A list containing
#' \itemize{
#' \item{\code{'data'}}{ generated series}
#' \item{\code{'shocks'}}{`q`-dimensional shock series}
#' \item{\code{'factors'}}{ `r`-dimensional factor series}
#' \item{\code{'K'}}{ transition matrix}
#' \item{\code{'loadings'}}{ factor loadings}
#' }
#' @export
#'
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @examples
#'     sim.factor.M1(100,10)
sim.factor.M1 <- function(n, p, q = 2, r = 4, do.scale = T, loadings=NULL, K = NULL){
  stopifnot(r>=q)
  burnin <- 100
  u2 <- matrix(rnorm(q * (n + burnin)), nrow = q)
  if(is.null(K) )K <- matrix(runif(r * q, -1, 1), nrow = r)
  D <- matrix(runif(r^2, -1, 1), nrow = r)
  D <- D/norm(D, type = '2')
  D <- runif(1, .4, .9) * D

  f <- matrix(0, nrow = r, ncol = n + burnin)
  f[, 1] <-  K %*% u2[, 1]
  for(t in 2:(n + burnin)) f[, t] <- D %*% f[, t - 1] + K %*% u2[, t]
  f <- f[, -(1:burnin)]
  u2 <- u2[, -(1:burnin)]

  if(is.null(loadings) )loadings <- matrix(runif(p * r, -1, 1), nrow = p)
  chi <- loadings %*% f

  if(do.scale) chi <- chi/apply(chi, 1, sd)
  return(list(data = chi, shocks = u2, factors = f, loadings = loadings, K=K))
}



#' Simulate data from a static factor model (Model 2)
#'
#' @param n sample size
#' @param p number of series
#' @param trunc.lags lag for MA representation
#' @param do.scale scale the output (default true)
#' @param a1,a2,alpha1,alpha2  generative parameters (default null, see reference)
#'
#' @return  A list containing
#' \itemize{
#' \item{\code{'data'}}{ generated series}
#' \item{\code{'shocks'}}{ 2-dimensional shock series}
#' \item{\code{'a1,a2,alpha1,alpha2'}}{ generative parameters}
#' }
#' @export
#'
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @examples
#'     sim.factor.M2(100,10)
sim.factor.M2 <- function(n, p, trunc.lags = 20, do.scale = T, a1 = NULL, a2 = NULL, alpha1 = NULL, alpha2 = NULL){
  stopifnot(n>=trunc.lags)
  chi <- matrix(0, p, n)
  u1 <- rnorm(n + trunc.lags)
  u2 <- rnorm(n + trunc.lags)
  uu <- rbind(u1,u2)
  a1 <- runif(p , -1,1)
  a2 <- runif(p , -1,1)
  alpha1 <- runif(p, -.8,.8)
  alpha2 <- runif(p, -.8,.8)
  for (jj in 1:p) {
    coeffs_1 <- alpha1[jj]*c(1,as.numeric( var.to.vma(as.matrix(a1[jj]), trunc.lags))) #gdfm::
    coeffs_2 <-  alpha2[jj]*c(1,as.numeric(var.to.vma(as.matrix(a2[jj]), trunc.lags)))
    for (t in 1:n) {
      chi[jj,t] <-  coeffs_1 %*% u1[(t+trunc.lags+1):t] + coeffs_2 %*% u2[(t+trunc.lags+1):t]
    }
  }
  if(do.scale) chi <- chi/apply(chi, 1, sd)
  return(list(data = chi, shocks = rbind(u1,u2),  a1 = a1, a2 = a2, alpha1 = alpha1, alpha2 = alpha2 ))
}



#' Simulate data from a (sparse) VAR(1) model
#'
#' @param n sample size
#' @param p number of series
#' @param A transition matrix, dimension p by p (default null)
#' @param cov generative covariance matrix (default identity)
#' @param prob edge probability
#' @param two_norm target 2-norm to scale A by (default NULL)
#' @param do.scale scale the output (default true)
#'
#' @return  A list containing
#' \itemize{
#' \item{\code{'data'}}{ generated series}
#' \item{\code{'A'}}{ transition matrix}
#' }
#' @export
#'
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @examples
#'     sim.idio(100,10, A=diag(0.3, 10))
sim.idio <- function(n, p, A = NULL, cov = diag(1,p), prob = 1/p, two_norm = NULL, do.scale = T){
  burnin <-100

  vep <- t(mvtnorm::rmvnorm(n + burnin, sigma = cov))
  if(is.null(A)) {
    A <- matrix(0, p, p)
    index <- sample(c(0,1), p^2, T, prob = c(1-prob,prob))
    A[which(index==1)] <- .275
  }
  if(!is.null(two_norm) ) A <- two_norm * A/ norm(A,"2")
  for(tt in 2:(n + burnin)) vep[, tt] <- vep[, tt] + A %*% vep[, tt - 1]
  vep <- vep[, -(1:burnin)]
  if(do.scale) vep <- vep/apply(vep, 1, sd)
  return(list(data = vep, A = A ))
}

