# simulate data

#' Simulate data from a dynamic factor model (Model 1) with factor number \code{r = q * lags}
#'
#' @param n sample size
#' @param p number of series
#' @param q dynamic dimension (default 2)
#' @param lags number of lags for which the observed series depends on the factor series (default 2)
#' @param do.scale scale the output (default \code{TRUE} )
#' @param D transition matrix, dimension q by q (default null)
#' @param loadings  loading matrix, dimension p by r (default null)
#'
#' @return  A list containing
#' \itemize{
#' \item{\code{'data'}}{ generated series}
#' \item{\code{'shocks'}}{ \code{q}-dimensional shock series}
#' \item{\code{'factors'}}{ \code{q}-dimensional factor series}
#' \item{\code{'D'}}{ transition matrix}
#' \item{\code{'loadings'}}{ factor loadings}
#' }
#' @export
#'
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @seealso \link{sim.factor.M2}
#' @examples
#'     sim.factor.M1(100,10)
sim.factor.M1 <- function(n, p, q = 2, lags = 2, do.scale = T, loadings=NULL, D = NULL){
  r <- q * lags
  burnin <- 100
  u2 <- matrix(rnorm(q * (n + burnin)), nrow = q)
  if(is.null(D)){
    D <- matrix(runif(q^2, 0, .3), nrow = q) #matrix(runif(r^2, -1, 1), nrow = r)
    diag(D) <- runif(q, .5,.8)
    D <- 0.7 * D/norm(D, type = '2')
  } else stopifnot(all(dim(D) == c(q,q)))

  #D <- runif(1, .4, .9) * D

  f <- matrix(0, nrow = q, ncol = n + burnin)
  f[, 1] <-  u2[, 1]
  for(t in 2:(n + burnin)) f[, t] <- D %*% f[, t - 1] +  u2[, t] #K %*%
  f <- f[, -(1:(burnin-lags+1) )]
  u2 <- u2[, -(1:burnin)]


  if(is.null(loadings))  loadings <- matrix(rnorm(p*r, 1,1), nrow = p) else stopifnot(all(dim(loadings) == c(p,r)))
  chi <- matrix(0, p, n)#
  for (ii in 1:lags) {
    chi <- chi + loadings[,(ii-1)*q + 1:q] %*% f[,1:n + lags - ii]
  }

  if(do.scale) chi <- chi/apply(chi, 1, sd)
  return(list(data = chi, shocks = u2, factors = f, loadings = loadings, D=D))
}



#' Simulate data from a static factor model (Model 2)
#'
#' @param n sample size
#' @param p number of series
#' @param trunc.lags lag for moving average representation
#' @param do.scale scale the output (default \code{TRUE} )
#' @param a1,a2,alpha1,alpha2  generative parameters (default null, see reference)
#'
#' @return  A list containing
#' \itemize{
#' \item{\code{'data'}}{ generated series}
#' \item{\code{'shocks'}}{ 2-dimensional shock series}
#' \item{\code{'a1','a2','alpha1','alpha2'}}{ generative parameters}
#' }
#' @export
#'
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @seealso \link{sim.factor.M1}
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
    coeffs_1 <- alpha1[jj]*c(1,as.numeric( var.to.vma(as.matrix(a1[jj]), trunc.lags)))
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
#' @param prob probability of an edge existing in the transition matrix, if \code{A} is \code{NULL}
#' @param two.norm target 2-norm to scale A by (default NULL)
#' @param do.scale scale the output (default \code{TRUE})
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
sim.idio <- function(n, p, A = NULL, cov = diag(1,p), prob = 1/p, two.norm = NULL, do.scale = T){
  burnin <-100
  prob <- max(0, min(1,prob))
  vep <- t(MASS::mvrnorm(n + burnin, mu = rep(0,p) , Sigma = cov))
  if(is.null(A)) {
    A <- matrix(0, p, p)
    index <- sample(c(0,1), p^2, T, prob = c(1-prob,prob))
    A[which(index==1)] <- .275
  }
  if(!is.null(two.norm) ) A <- two.norm * A/ norm(A,"2")
  for(tt in 2:(n + burnin)) vep[, tt] <- vep[, tt] + A %*% vep[, tt - 1]
  vep <- vep[, -(1:burnin)]
  if(do.scale) vep <- vep/apply(vep, 1, sd)
  return(list(data = vep, A = A ))
}

