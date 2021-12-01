# simulate data

#' Simulate data from a factor model with MA loadings
#'
#' @param n sample size
#' @param p number of series
#' @param q dynamic dimension (default 2)
#' @param r factor number (default 4)
#' @param do_scale scale the output (default true)
#' @param K  generative parameter matrix, dimension r by q (default null)
#' @param lam  loading matrix, dimension p by r (default null)
#'
#' @return list of data and generative parameters
#' @export
#'
#' @examples simMAfactor(100,10)
simMAfactor <- function(n, p, q = 2, r = 4, do_scale = T, K = NULL, lam=NULL){
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

  if(is.null(lam) )lam <- matrix(runif(p * r, -1, 1), nrow = p)
  chi <- lam %*% f

  if(do_scale) chi <- chi/apply(chi, 1, sd)
  return(list(data = chi, shocks = u2, factors = f, loadings = lam, K=K))
}



#' Simulate data from a factor model with AR loadings
#'
#' @param n sample size
#' @param p number of series
#' @param MA_lag lag for MA representation
#' @param do_scale scale the ouptut (default true)
#' @param a1,a2,alpha1,alpha2  generative parameters (default null)
#'
#' @return list of data and generative parameters
#' @export
#'
#' @examples simARfactor(100,10)
simARfactor <- function(n, p, MA_lag = 20, do_scale = T, a1 = NULL, a2 = NULL, alpha1 = NULL, alpha2 = NULL){
  stopifnot(n>=MA_lag)
  chi <- matrix(0, p, n)
  u1 <- rnorm(n + MA_lag)
  u2 <- rnorm(n + MA_lag)
  uu <- rbind(u1,u2)
  a1 <- runif(p , -1,1)
  a2 <- runif(p , -1,1)
  alpha1 <- runif(p, -.8,.8)
  alpha2 <- runif(p, -.8,.8)
  for (jj in 1:p) {
    coeffs_1 <- alpha1[jj]*c(1,as.numeric( ARtoMA(array(a1[jj], c(1,1,1)), nlag = MA_lag))) #gdfm::
    coeffs_2 <-  alpha2[jj]*c(1,as.numeric(ARtoMA(array(a2[jj], c(1,1,1)), nlag = MA_lag)))
    for (t in 1:n) {
      chi[jj,t] <-  coeffs_1 %*% u1[(t+MA_lag):t] + coeffs_2 %*% u2[(t+MA_lag):t]
    }
  }
  if(do_scale) chi <- chi/apply(chi, 1, sd)
  return(list(data = chi, shocks = rbind(u1,u2),  a1 = a1, a2 = a2, alpha1 = alpha1, alpha2 = alpha2 ))
}



#' Simulate data from a VAR(1) model
#'
#' @param n sample size
#' @param p number of series
#' @param A generative parameter matrix, dimension p by p (default null)
#' @param cov generative covariance matrix (default identity)
#' @param prob edge probability 
#' @param two_norm target 2-norm to scale A by (default 0.7, specify NULL to avoid scaling)
#' @param do_scale scale the output (default true)
#'
#' @return list of data and A
#' @export
#'
#' @examples simARidio(100,10, A=diag(0.3, 10))
simARidio <- function(n, p, A = NULL, cov = diag(1,p), prob = 1/p, two_norm = NULL, do_scale = T){
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
  if(do_scale) vep <- vep/apply(vep, 1, sd)
  return(list(data = vep, A = A ))
}

