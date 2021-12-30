#' @title Nonparametric estimation of long-run partial correlations of factor-adjusted VAR processes
#' @description Returns a nonparametric estimate of long-run partial correlations of the VAR process
#' from the inverse of long-run covariance matrix obtained via constrained \code{l1}-minimisation.
#' @param object \code{fnets} object
#' @param x input time series matrix; with each row representing a variable
#' @param eta regularisation parameter; if \code{eta = NULL}, it is selected by cross validation
#' @param cv.args a list specifying arguments for the cross validation procedure
#' for selecting the tuning parameter involved in long-run partial correlation matrix estimation. It contains:
#' \itemize{
#'    \item{n.folds}{ number of folds}
#'    \item{path.length}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{do.plot}{ whether to plot the output of the cross validation step}
#' }
#' @param correct.zero whether to correct for any zero-entries in the diagonals of the inverse of long-run covariance matrix
#' @param n.cores number of cores to use for parallel communing, see \code{\link[parallel]{makePSOCKcluster}}
#' @return a list containing
#' \itemize{
#' \item{Omega}{ estimated inverse of the long-run covariance matrix}
#' \item{lrpc}{ estimated long-run partial correlation matrix}
#' \item{eta}{ regularisation parameter}
#' }
#' @example R/examples/nonparlrpc.R
#' @export
npar.lrpc <- function(object, x, eta = NULL,
                      cv.args = list(n.folds = 1, path.length = 10, do.plot = FALSE),
                      correct.zero = TRUE, n.cores = min(parallel::detectCores() - 1, 3)){

  xx <- x - object$mean.x
  p <- dim(x)[1]
  GG <- Re(object$spec$Sigma_i[,, 1])

  if(is.null(eta)){
    dcv <- direct.cv(object, xx, target = 'spec', symmetric = 'min',
                     n.folds = cv.args$n.folds, path.length = cv.args$path.length, 
                     q = object$q, kern.const = object$kern.const, n.cores = n.cores,
                     do.plot = cv.args$do.plot)
    eta <- dcv$eta
  }
  DD <- direct.inv.est(GG, eta = eta, symmetric = 'min',
                       correct.zero = correct.zero, n.cores = n.cores)$DD
  lrpc <- - t(t(DD)/sqrt(diag(DD)))/sqrt(diag(DD))
  out <- list(Omega = DD, lrpc = lrpc, eta = eta)

  return(out)

}

#' @title Parametric estimation of long-run partial correlations of factor-adjusted VAR processes
#' @description Returns a parametric estimate of long-run partial correlations of the VAR process
#' from the VAR parameter estimates and the inverse of innovation covariance matrix obtained via constrained \code{l1}-minimisation.
#' @details See Barigozzi, Cho and Owens (2021) for further details.
#' @param object \code{fnets} object
#' @param x input time series matrix; with each row representing a variable
#' @param eta regularisation parameter; if \code{eta = NULL}, it is selected by cross validation
#' @param cv.args a list specifying arguments for the cross validation procedure
#' for selecting the tuning parameter involved in long-run partial correlation matrix estimation. It contains:
#' \itemize{
#'    \item{n.folds}{ number of folds}
#'    \item{path.length}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{do.plot}{ whether to plot the output of the cross validation step}
#' }
#' @param correct.zero whether to correct for any zero-entries in the diagonals of the inverse of long-run covariance matrix
#' @param n.cores number of cores to use for parallel communing, see \code{\link[parallel]{makePSOCKcluster}}
#' @return a list containing
#' \itemize{
#' \item{Delta}{ estimated inverse of the innovation covariance matrix}
#' \item{Omega}{ estimated inverse of the long-run covariance matrix}
#' \item{pc}{ estimated innovation partial correlation matrix}
#' \item{lrpc}{ estimated long-run partial correlation matrix}
#' \item{eta}{ regularisation parameter}
#' }
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @example R/examples/paramlrpc.R
#' @export
par.lrpc <- function(object, x, eta = NULL,
                       cv.args = list(n.folds = 1, path.length = 10, do.plot = FALSE),
                       correct.zero = TRUE,
                       n.cores = min(parallel::detectCores() - 1, 3)){

  xx <- x - object$mean.x
  p <- dim(x)[1]

  GG <- object$idio.var$Gamma
  A <- t(object$idio.var$beta)
  d <- dim(A)[2]/p

  A1 <- diag(1, p)
  for(ll in 1:d) A1 <- A1 - A[, (ll - 1) * p + 1:p]

  if(is.null(eta)){
    dcv <- direct.cv(object, xx, target = 'acv', symmetric = 'min',
                     n.folds = cv.args$n.folds, path.length = cv.args$path.length, 
                     q = object$q, kern.const = object$kern.const, n.cores = n.cores,
                     do.plot = cv.args$do.plot)
    eta <- dcv$eta
  }
  Delta <- direct.inv.est(GG, eta = eta, symmetric = 'min',
                          correct.zero = correct.zero, n.cores = n.cores)$DD
  Omega <- 2 * pi * t(A1) %*% Delta %*% A1
  pc <- - t(t(Delta)/sqrt(diag(Delta)))/sqrt(diag(Delta))
  lrpc <- - t(t(Omega)/sqrt(diag(Omega)))/sqrt(diag(Omega))
  out <- list(Delta = Delta, Omega = Omega, pc = pc, lrpc = lrpc, eta = eta)

  return(out)

}

#' @keywords internal
direct.cv <- function(object, xx, target = c('spec', 'acv'), symmetric = c('min', 'max', 'avg', 'none'),
                      n.folds = 1, path.length = 10, q = 0, kern.const = 4, n.cores = min(parallel::detectCores() - 1, 3),
                      do.plot = FALSE){

  n <- ncol(xx)
  p <- nrow(xx)

  target <- match.arg(target, c('spec', 'acv'))
  if(target == 'spec'){
    GG <- Re(object$spec$Sigma_i[,, 1])
    eta.max <- max(abs(GG))
    eta.path <- round(exp(seq(log(eta.max), log(eta.max * .01), length.out = path.length)), digits = 10)
  }
  if(target == 'acv'){
    A <- t(object$idio.var$beta)
    d <- dim(A)[2]/p
    GG <- object$idio.var$Gamma
    eta.max <- max(abs(GG))
    eta.path <- round(exp(seq(log(eta.max), log(eta.max * .01), length.out = path.length)), digits = 10)
  }

  cv.err <- rep(0, length = path.length)
  ind.list <- split(1:n, ceiling(n.folds*(1:n)/n))
  for(fold in 1:n.folds){
    train.ind <- 1:ceiling(length(ind.list[[fold]]) * .5)
    train.x <- xx[, ind.list[[fold]][train.ind]]
    test.x  <- xx[, ind.list[[fold]][- train.ind]]
    if(target == 'spec'){
      train.GG <- Re(dyn.pca(train.x, q = q, kern.const = kern.const)$spec$Sigma_i[,, 1])
      test.GG <- Re(dyn.pca(test.x, q = q, kern.const = kern.const)$spec$Sigma_i[,, 1])
    }
    if(target == 'acv'){
      train.G0 <- dyn.pca(train.x, q = q, kern.const = kern.const)$acv$Gamma_i
      test.G0 <- dyn.pca(test.x, q = q, kern.const = kern.const)$acv$Gamma_i
      train.GG <- train.G0[,, 1]
      test.GG <- test.G0[,, 1]
      for(ll in 1:d){
        train.GG <- train.GG - A[, (ll - 1) * p + 1:p] %*% train.G0[,, ll + 1]
        test.GG <- test.GG - A[, (ll - 1) * p + 1:p] %*% test.G0[,, ll + 1]
      }
    }

    for(ii in 1:path.length){
      DD <- direct.inv.est(train.GG, eta = eta.path[ii], symmetric = symmetric, n.cores = n.cores)$DD
      DG <- DD %*% test.GG
      sv <- svd(DG, nu = 0, nv = 0)
      cv.err[ii] <- cv.err[ii] + sum(sv$d) - sum(log(sv$d)) - p
    }
  }

  eta.min <- eta.path[which.min(cv.err)]

  if(do.plot){
    plot(eta.path, cv.err, type = 'b', col = 2, pch = 2, log = 'x', 
         xlab = 'eta (log scale)', ylab = 'CV error', main = 'CV for (LR)PC matrix estimation')
    abline(v = eta.min)
  }
  
  out <- list(eta = eta.min, cv.error = cv.err, eta.path = eta.path)
  return(out)

}

#' @keywords internal
direct.inv.est <- function(GG, eta = NULL, symmetric = c('min', 'max',  'avg', 'none'),
                           correct.zero = FALSE,
                           n.cores = min(parallel::detectCores() - 1, 3)){

  p <- dim(GG)[1]
  f.obj <- rep(1, 2 * p)
  f.con <- rbind(-GG, GG)
  f.con <- cbind(f.con,-f.con)
  f.dir <- rep('<=', 2 * p)

  cl <- parallel::makePSOCKcluster(n.cores)
  doParallel::registerDoParallel(cl)

  DD <- foreach::foreach(ii = 1:p, .combine = 'cbind', .multicombine = TRUE, .export = c('lp')) %dopar% {
    ee <- rep(0, p)
    ee[ii] <- 1
    b1 <- rep(eta, p) - ee
    b2 <- rep(eta, p) + ee
    f.rhs <- c(b1, b2)
    lpout <- lpSolve::lp('min', f.obj, f.con, f.dir, f.rhs)
    lpout$solution[1:p] - lpout$solution[-(1:p)]
  }
  parallel::stopCluster(cl)

  DD <- make.symmetric(DD, symmetric)

  if(correct.zero){
    tmp <- gen.inverse(GG)
    ind <- which(diag(DD) == 0)
    diag(DD)[ind] <- tmp[ind]
  }

  out <- list(DD = DD, eta = eta, symmetric = symmetric)
  return(out)

}

#' @keywords internal
gen.inverse <- function(GG){

  p <- dim(GG)[1]
  sv <- svd(GG)
  L <- GG * 0
  diag(L)[sv$d > 0] <- sv$d[sv$d > 0]
  return(diag(sv$u %*% L %*% t(sv$u)))

}
