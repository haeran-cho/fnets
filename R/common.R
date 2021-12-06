#' @title Blockwise VAR estimation under GDFM
#' @description internal function
#' @references Forni, M., Hallin, M., Lippi, M., & Zaffaroni, P. (2017). Dynamic factor models with infinite-dimensional factor space: Asymptotic analysis. Journal of Econometrics, 199(1), 74--92.
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @keywords internal
# #' @export
common.var.estimation <- function(xx, Gamma_c, q, var.order = NULL, max.var.order = NULL, trunc.lags, n.perm){

  n <- dim(xx)[2]; p <- dim(xx)[1]
  mm <- (dim(Gamma_c)[3] - 1)/2
  N <- p %/% (q + 1)
  if(!is.null(var.order)) max.var.order <- var.order
  if(is.null(max.var.order)) max.var.order <- min(max(1, ceiling(10 * log(n, 10)/(q + 1)^2)), 10)

  if(q < 1) warning(paste0('There should be at least one factor for common component estimation!'))

  if(q >= 1){
    # approach 1: for each permutation, obtain IRF and u (with cholesky identification), average the output, what FHLZ 2017 do
    irf.array <- array(0, dim = c(p, q, trunc.lags + 2, n.perm))
    u.array <- array(0, dim = c(q, n, n.perm))

    # # approach 2: average A^(-1), average Z, obtain IRF and u (with cholesky identification)
    # inv.var.est <- array(0, dim = c(p, p, trunc.lags + 2))
    # res <- array(0, dim = c(p, n))
    # res[, 1:var.order] <- NA

    for(ii in 1:n.perm){
      if(ii == 1) perm.index <- 1:p else perm.index <- sample(p, p)
      Gamma_c_perm <- Gamma_c[perm.index, perm.index, ]

      A <- list()
      z <- xx[perm.index, ]; z[, 1:max.var.order] <- NA
      for(jj in 1:N){
        if(jj == N) block <- ((jj - 1) * (q + 1) + 1):p else block <- ((jj - 1) * (q + 1) + 1):(jj * (q + 1))
        pblock <- perm.index[block]
        nblock <- length(block)

        if(is.null(var.order)){
          bic <- common.bic(Gamma_c_perm, block, n, max.var.order)
          s <- which.min(bic[-1])
          # print(s)
        } else s <- var.order

        tmp <- common.yw.est(Gamma_c_perm, block, s)$A
        A <- c(A, list(tmp))
        for(ll in 1:s) z[block, (max.var.order + 1):n] <- z[block, (max.var.order + 1):n] - tmp[, nblock * (ll - 1) + 1:nblock] %*% xx[pblock, (max.var.order + 1):n - ll]
      }
      # res[perm.index, (var.order + 1):n] <- res[perm.index, (var.order + 1):n] + z[, (var.order + 1):n]/n.perm

      zz <- z[, (max.var.order + 1):n] %*% t(z[, (max.var.order + 1):n]) / (n - max.var.order)
      svz <- svd(zz, nu = q, nv = 0) # plot(svz$d)
      R <- as.matrix(t(t(svz$u) * sqrt(svz$d[1:q])))
      u <- t(svz$u) %*% z / sqrt(svz$d[1:q])

      tmp.irf <- irf.array[,,, 1, drop=FALSE] * 0
      for(jj in 1:N){
        if(jj == N) block <- ((jj - 1) * (q + 1) + 1):p else block <- ((jj - 1) * (q + 1) + 1):(jj * (q + 1))
        pblock <- perm.index[block]
        invA <- var.to.vma(A[[jj]], trunc.lags + 1)
        for(ll in 1:(trunc.lags + 2)){
          # inv.var.est[pblock, pblock, ll] <- inv.var.est[pblock, pblock, ll] + invA[,, ll]/n.perm
          tmp.irf[pblock,, ll,  ] <- invA[,, ll] %*% R[block, ]
        }
      }
      B0 <- tmp.irf[1:q,, 1,]
      if(all(B0==0)) H<- as.matrix(B0) else {
        C0 <- t(chol(B0%*%t(B0))) # cholesky identification
        H <- solve(B0) %*% C0
      }
      for(ll in 1:(trunc.lags + 2)) irf.array[,, ll, ii ] <- tmp.irf[,, ll,] %*% H
      u.array[,, ii] <- t(H) %*% u
    }

    # irf.est <- array(0, dim = c(dim(irf.array)[1:3], 2))
    # u.est <- array(0, dim = c(q, n, 2))

    # approach 1
    irf.est <- apply(irf.array, c(1, 2, 3), mean)
    u.est <- apply(u.array, c(1, 2), mean)

    # # approach 2
    # zz <- res[, (var.order + 1):n] %*% t(res[, (var.order + 1):n]) / (n - var.order)
    # svz <- svd(zz, nu = q, nv = 0)
    # R <- t(t(svz$u) * sqrt(svz$d[1:q]))
    # u <- t(svz$u) %*% z / sqrt(svz$d[1:q])
    #
    # tmp.irf <- irf.est[,,, 2] * 0
    # for(ll in 1:(trunc.lags + 2)) tmp.irf[,, ll] <- inv.var.est[,, ll] %*% R
    # B0 <- tmp.irf[1:q,, 1]
    # C0 <- t(chol(B0%*%t(B0))) # cholesky identification
    # H <- solve(B0) %*% C0
    # for(ll in 1:(trunc.lags + 2)) irf.est[,, ll, 2] <- tmp.irf[,, ll] %*% H
    # u.est[,, 2] <- t(H) %*% u

    out <- list(irf.array = irf.array, u.array = u.array, irf.est = irf.est, u.est = u.est)
    return(out)
  }
}

#' @title Prediction for the factor-driven common component
#' @description Predicts common component from a \code{fnets} object for new data
#' @param object \code{fnets} object
#' @param x input time series matrix, with each row representing a time series
#' @param h forecast horizon
#' @param common.method which of "static" or "var" to forecast the common component with
#' @param r factor number, if r=NULL this is selected using the maximal eigenratio
#' @return A list containing
#' \itemize{
#' \item{\code{'is'}}{ x in-sample estimation}
#' \item{\code{'fc'}}{ x forecast}
#' \item{\code{'r'}}{ factor number}
#' \item{\code{'h'}}{ forecast horizon}
#' }
#' @example examples/predict.R
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @references Forni, M., Hallin, M., Lippi, M., & Reichlin, L. (2005). The generalized dynamic factor model: one-sided estimation and forecasting. Journal of the American Statistical Association, 100(471), 830--840.
#' @references Forni, M., Hallin, M., Lippi, M., & Zaffaroni, P. (2017). Dynamic factor models with infinite-dimensional factor space: Asymptotic analysis. Journal of Econometrics, 199(1), 74--92.
#' @export
common.predict <- function(object, x, h = 1, common.method = c('static', 'var'), r = NULL){

  xx <- x - object$mean.x
  p <- dim(x)[1]
  common.method <- match.arg(common.method, c('static', 'var'))
  if(object$q < 1){
    warning(paste0('There should be at least one factor for common component estimation!'))
    pre <- list(is = 0 * x, fc = matrix(0, nrow = p, ncol = h))
  }
  if(object$q >= 1){
    if(common.method == 'static') pre <- common.static.predict(xx = xx, Gamma_c = object$acv$Gamma_c, q = object$q, r = r, h = h)
    if(common.method == 'var') pre <- common.var.predict(xx = xx, cve = object$common.var, h = h)
  }
  return(pre)

}

# #' @export
#' @description internal function
#' @keywords internal
common.static.predict <- function(xx, Gamma_c, q, r = NULL, max.r = NULL, h = 1){

  p <- dim(xx)[1]; n <- dim(xx)[2]
  if(is.null(max.r)) max.r <- max(q, min(50, round(sqrt(min(n, p)))))
  if(h >= dim(Gamma_c)[3]){
    warning(paste0('At most ', (dim(Gamma_c)[3] - 1)/2, '-step ahead forecast is available!'))
    h <- (dim(Gamma_c)[3] - 1)/2
  }

  sv <- svd(Gamma_c[,, 1], nu = max.r, nv = 0)
  if(is.null(r)) r <- which.max(sv$d[q:max.r]/sv$d[1 + q:max.r]) + q - 1
                                #if(is.null(r)) r <- which.max(sv$d[1:max.r]/sv$d[1 + 1:max.r])

  is <- sv$u[, 1:r, drop = FALSE] %*% t(sv$u[, 1:r, drop = FALSE]) %*% xx
  if(h >= 1){
    fc <- matrix(0, nrow = p, ncol = h)
    proj.x <- t(t(sv$u[, 1:r, drop = FALSE])/sv$d[1:r]) %*% t(sv$u[, 1:r, drop = FALSE]) %*% xx[, n]
    for(hh in 1:h) fc[, hh] <- t(Gamma_c[,, hh + 1]) %*% proj.x
  } else fc <- NA

  out <- list(is = is, fc = fc, r = r)
  return(out)

}

# #' @export
#' @description internal function
#' @keywords internal
common.var.predict <- function(xx, cve, h = 1){
  p <- dim(xx)[1]; n <- dim(xx)[2]
  trunc.lags <- dim(cve$irf.est)[3]
  if(h >= trunc.lags + 1){
    warning(paste0('At most ', trunc.lags, '-step ahead forecast is available!'))
    h <- trunc.lags
  }

  irf <- cve$irf.est
  u <- cve$u.est
  trunc.lags <- dim(irf)[3] - 1

  is <- xx * 0
  is[, 1:trunc.lags] <- NA
  for(ll in 1:(trunc.lags + 1)) is[, (trunc.lags + 1):n] <- is[, (trunc.lags + 1):n] + as.matrix(irf[,, ll]) %*% u[, (trunc.lags + 1):n - ll + 1,drop=FALSE]

  if(h >= 1){
    fc <- matrix(0, nrow = p, ncol = h)
    for(hh in 1:h) for(ll in 1:(trunc.lags + 1 - hh)) fc[, hh] <- fc[, hh] + as.matrix(irf[,, ll+hh]) %*% u[, n - ll + 1,drop=FALSE]
  } else fc <- NA

  out <- list(is = is, fc = fc, h = h)
  return(out)

}

# #' @export
#' @description internal function
#' @keywords internal
var.to.vma <- function(A, trunc.lags){

  d <- dim(A)[1]; s <- dim(A)[2]
  l <- s/d
  B <- array(0, dim = c(d, d, trunc.lags + 1))
  B[,, 1] <- diag(1, d)
  for(ii in 1:trunc.lags){
    for(jj in 1:min(ii, l)) B[,, ii + 1] <- B[,, ii + 1] + B[,, ii - jj + 1] %*% A[, (jj - 1) * d + 1:d]
  }
  B

}

 # #' @export
#' @description internal function
#' @keywords internal
common.yw.est <- function(Gcp, block, var.order){
  nblock <- length(block)
  B <- matrix(0, nrow = nblock, ncol = nblock * var.order)
  C <- matrix(0, nrow = nblock * var.order, ncol = nblock * var.order)
  for(ll in 1:var.order){
    B[, nblock * (ll - 1) + 1:nblock] <- t(Gcp[block, block, 1 + ll])
    for(lll in 1:var.order){
      if(ll >= lll){
        C[nblock * (ll - 1) + 1:nblock, nblock * (lll - 1) + 1:nblock] <- Gcp[block, block, 1 + ll - lll]
      } else C[nblock * (ll - 1) + 1:nblock, nblock * (lll - 1) + 1:nblock] <- t(Gcp[block, block, 1 + lll - ll])
    }
  }
  A <- B %*% solve(C, symmetric = TRUE)
  out <- list(A = A, B = B, C = C)
  return(out)
}

# #' @export
#' @description internal function
#' @keywords internal
common.bic <- function(Gcp, block, len, max.var.order = 5){

  nblock <- length(block)
  bic <- rep(0, max.var.order + 1)
  bic[1] <- log(det(Gcp[block, block, 1]))
  for(ii in 1:max.var.order){
    cye <- common.yw.est(Gcp, block, ii)
    G0 <- Gcp[block, block, 1] - cye$B %*% t(cye$A) - cye$A %*% t(cye$B) + cye$A %*% cye$C %*% t(cye$A)
    bic[ii + 1] <- log(det(G0)) + 2 * log(len) * ii * nblock^2/len
  }
  bic

}
