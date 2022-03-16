#' @title Nonparametric estimation of long-run partial correlations of factor-adjusted VAR processes
#' @description Returns a nonparametric estimate of long-run partial correlations of the VAR process
#' from the inverse of long-run covariance matrix obtained via constrained \code{l1}-minimisation.
#' @param object \code{fnets} object
#' @param x input time series matrix; with each row representing a variable
#' @param eta regularisation parameter; if \code{eta = NULL}, it is selected by cross validation
#' @param cv.args a list specifying arguments for the cross validation procedure
#' for selecting the tuning parameter involved in long-run partial correlation matrix estimation. It contains:
#' \itemize{
#'    \item{\code{n.folds}}{ number of folds}
#'    \item{\code{path.length}}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{\code{do.plot}}{ whether to plot the output of the cross validation step}
#' }
#' @param do.correct whether to correct for any negative entries in the diagonals of the inverse of long-run covariance matrix
#' @param n.cores number of cores to use for parallel computing, see \link[parallel]{makePSOCKcluster}
#' @return a list containing
#' \item{Omega}{ estimated inverse of the long-run covariance matrix}
#' \item{lrpc}{ estimated long-run partial correlation matrix}
#' \item{eta}{ regularisation parameter}
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.common1(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#' out <- fnets(x, q = NULL, idio.method = 'lasso', lrpc.method = 'none')
#' nlrpc <- npar.lrpc(out, x, cv.args = list(n.folds = 1, path.length = 10, do.plot = TRUE))
#' out$lrpc <- nlrpc
#' out$lrpc.method <- 'npar'
#' plot(out, type = 'lrpc', display = 'heatmap', threshold = .05)
#' }
#' @importFrom parallel detectCores
#' @export
npar.lrpc <- function(object, x, eta = NULL,
                      cv.args = list(n.folds = 1, path.length = 10, do.plot = FALSE),
                      do.correct = TRUE, n.cores = min(parallel::detectCores() - 1, 3)){

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
                       do.correct = do.correct, n.cores = n.cores)$DD
  lrpc <- - t(t(DD)/sqrt(diag(DD)))/sqrt(diag(DD))
  out <- list(Omega = DD, lrpc = lrpc, eta = eta)

  return(out)

}

#' @title Parametric estimation of long-run partial correlations of factor-adjusted VAR processes
#' @description Returns a parametric estimate of long-run partial correlations of the VAR process
#' from the VAR parameter estimates and the inverse of innovation covariance matrix obtained via constrained \code{l1}-minimisation.
#' @details See Barigozzi, Cho and Owens (2021) for further details, and Cai, Liu and Zhou (2016) for further details on the adaptive estimation procedure.
#' @param object \code{fnets} object
#' @param x input time series matrix; with each row representing a variable
#' @param eta regularisation parameter; if \code{eta = NULL}, it is selected by cross validation
#' @param cv.args a list specifying arguments for the cross validation procedure
#' for selecting the tuning parameter involved in long-run partial correlation matrix estimation. It contains:
#' \itemize{
#'    \item{\code{n.folds}}{ number of folds}
#'    \item{\code{path.length}}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{\code{do.plot}}{ whether to plot the output of the cross validation step}
#' }
#' @param adaptive whether to use the adaptive estimation procedure
#' @param eta.1 regularisation parameter for Step 1 of the adaptive estimation procedure; if \code{eta.1 = NULL}, defaults to \code{2 * sqrt(log(dim(x)[1])/dim(x)[2])}
#' @param do.correct whether to correct for any negative entries in the diagonals of the inverse of long-run covariance matrix
#' @param n.cores number of cores to use for parallel computing, see \link[parallel]{makePSOCKcluster}
#' @return a list containing
#' \item{Delta}{ estimated inverse of the innovation covariance matrix}
#' \item{Omega}{ estimated inverse of the long-run covariance matrix}
#' \item{pc}{ estimated innovation partial correlation matrix}
#' \item{lrpc}{ estimated long-run partial correlation matrix}
#' \item{eta}{ regularisation parameter}
#' \item{adaptive}{ was the adaptive procedure used}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#'
#' Cai, T. T., Liu, W., & Zhou, H. H. (2016). Estimating sparse precision matrix: Optimal rates of convergence and adaptive estimation. The Annals of Statistics, 44(2), 455-488.
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.common1(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#' out <- fnets(x, q = NULL, idio.method = 'lasso', lrpc.method = 'none')
#' plrpc <- par.lrpc(out, x, cv.args = list(n.folds = 1, path.length = 10, do.plot = TRUE))
#' out$lrpc <- plrpc
#' out$lrpc.method <- 'par'
#' plot(out, type = 'pc', display = 'network', threshold = .05)
#' plot(out, type = 'lrpc', display = 'heatmap', threshold = .05)
#' }
#' @importFrom parallel detectCores
#' @export
par.lrpc <- function(object, x, eta = NULL,
                     cv.args = list(n.folds = 1, path.length = 10, do.plot = FALSE),
                     adaptive = FALSE, eta.1 = NULL,
                     do.correct = TRUE,
                     n.cores = min(parallel::detectCores() - 1, 3)){

  xx <- x - object$mean.x
  p <- dim(x)[1]
  n <- dim(x)[2]

  GG <- object$idio.var$Gamma
  A <- t(object$idio.var$beta)
  d <- dim(A)[2]/p

  A1 <- diag(1, p)
  for(ll in 1:d) A1 <- A1 - A[, (ll - 1) * p + 1:p]

  if(is.null(eta)){
    dcv <- direct.cv(object, xx, target = 'acv', symmetric = 'min',
                     n.folds = cv.args$n.folds, path.length = cv.args$path.length,
                     q = object$q, kern.const = object$kern.const, n.cores = n.cores,
                     adaptive = adaptive, eta.1 = eta.1,
                     do.plot = cv.args$do.plot)
    eta <- dcv$eta
  }
  if(adaptive){
    if(is.null(eta.1)){
      eta.1 <- 2 * sqrt(log(p)/n)
    }
    Delta <- adaptive.direct.inv.est(GG, n, eta = eta, eta.1 = eta.1, symmetric = 'min',
                                     do.correct = do.correct, n.cores = n.cores)$DD
  } else Delta <- direct.inv.est(GG, eta = eta, symmetric = 'min',
                          do.correct = do.correct, n.cores = n.cores)$DD
  Omega <- 2 * pi * t(A1) %*% Delta %*% A1
  if(do.correct) Omega <- correct.diag(Re(object$spec$Sigma_i[,, 1]), Omega)

  pc <- - t(t(Delta)/sqrt(diag(Delta)))/sqrt(diag(Delta))
  lrpc <- - t(t(Omega)/sqrt(diag(Omega)))/sqrt(diag(Omega))
  out <- list(Delta = Delta, Omega = Omega, pc = pc, lrpc = lrpc, eta = eta, adaptive = adaptive)

  return(out)

}

#' @keywords internal
#' @importFrom parallel detectCores
#' @importFrom graphics abline
direct.cv <- function(object, xx, target = c('spec', 'acv'), symmetric = c('min', 'max', 'avg', 'none'),
                      n.folds = 1, path.length = 10, q = 0, kern.const = 4, n.cores = min(parallel::detectCores() - 1, 3),
                      adaptive = FALSE, eta.1 = NULL,
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
    if(adaptive) eta.max.2 <- 2*sqrt(log(p)/n) else eta.max.2 <- eta.max
    eta.path <- round(exp(seq(log(eta.max.2), log(eta.max * .01), length.out = path.length)), digits = 10)
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
      train.G0 <- dyn.pca(train.x, q = q, kern.const = kern.const, mm = d)$acv$Gamma_i
      test.G0 <- dyn.pca(test.x, q = q, kern.const = kern.const, mm = d)$acv$Gamma_i
      train.GG <- train.G0[,, 1]
      test.GG <- test.G0[,, 1]
      for(ll in 1:d){
        train.GG <- train.GG - A[, (ll - 1) * p + 1:p] %*% train.G0[,, ll + 1]
        test.GG <- test.GG - A[, (ll - 1) * p + 1:p] %*% test.G0[,, ll + 1]
      }
    }

    for(ii in 1:path.length){
      if(adaptive) DD <- adaptive.direct.inv.est(train.GG, n=n, eta = eta.path[ii], eta.1 = eta.1, symmetric = symmetric, n.cores = n.cores)$DD
        else DD <- direct.inv.est(train.GG, eta = eta.path[ii], symmetric = symmetric, n.cores = n.cores)$DD
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
#' @importFrom parallel makePSOCKcluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lpSolve lp
direct.inv.est <- function(GG, eta = NULL, symmetric = c('min', 'max',  'avg', 'none'),
                           do.correct = FALSE,
                           n.cores = min(parallel::detectCores() - 1, 3)){

  p <- dim(GG)[1]
  f.obj <- rep(1, 2 * p)
  f.con <- rbind(-GG, GG)
  f.con <- cbind(f.con,-f.con)
  f.dir <- rep('<=', 2 * p)

  cl <- parallel::makePSOCKcluster(n.cores)
  doParallel::registerDoParallel(cl)

  ii <- 1
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
  if(do.correct) DD <- correct.diag(GG, DD)

  out <- list(DD = DD, eta = eta, symmetric = symmetric)
  return(out)

}

#' @keywords internal
#' @importFrom parallel makePSOCKcluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lpSolve lp
adaptive.direct.inv.est <- function(GG, n, eta = NULL, eta.1 = NULL, symmetric = c('min', 'max',  'avg', 'none'),
                                    do.correct = FALSE,
                                    n.cores = min(parallel::detectCores() - 1, 3)){
  p <- dim(GG)[1]
  f.obj <- rep(1, 2 * p)
  GG.n <- GG + diag(1/n,p) #add ridge
  dGG <- pmax(diag(GG),0)
  f.dir <- rep('<=', 2 * p)

  f.con.0 <- rbind(-GG.n, GG.n) #initialise
  f.con.0 <- cbind(f.con.0,-f.con.0)
  ## Step 1 //
  cl <- parallel::makePSOCKcluster(n.cores)
  doParallel::registerDoParallel(cl)
  if(is.null(eta.1))  eta.1 <- 2 * sqrt(log(p)/n)
  ii <- 1
  step1.index <- which(dGG <= sqrt(n/log(p)))
  f.con.1 <- rbind(f.con.0, 0)
  f.dir.1 <- c(f.dir, "==") #diagonals are positive
  DD.1 <- foreach::foreach(ii = step1.index, .combine = 'cbind', .multicombine = TRUE, .export = c('lp')) %dopar% {
    f.con.ii <- f.con.1
    ii.replace <- eta.1 * pmax(dGG, dGG[ii])
    f.con.ii[1:(2*p),ii] <- f.con.ii[1:(2*p),ii] - ii.replace  #mutate cols
    f.con.ii[1:(2*p),ii+p] <- f.con.ii[1:(2*p),ii+p] - ii.replace  #mutate cols
    f.con.ii[2*p+1,ii+p] <- 1 #diagonals are positive
    ee <- rep(0,p)
    ee[ii] <- 1
    b1 <- -ee
    b2 <- ee
    f.rhs <-  c(b1, b2,0)
    lpout <- lpSolve::lp('min', f.obj, f.con.ii, f.dir.1, f.rhs)
    lpout$solution[1:p] - lpout$solution[p+(1:p)]
  }
  dDD.1 <- diag(DD.1)
  dDD.1[!step1.index] <- sqrt(log(p)/n)
  ## Step 2 //
  if(is.null(eta)){
    eta <- 2 * sqrt(log(p)/n)
  }
  ii <- 1
  DD.2 <- foreach::foreach(ii = 1:p, .combine = 'cbind', .multicombine = TRUE, .export = c('lp')) %dopar% {
    ee <- rep(0, p)
    ee[ii] <- 1
    bb <- eta*sqrt(dGG)*sqrt(dDD.1[ii])
    b1 <- bb - ee
    b2 <- bb + ee
    f.rhs <- c(b1, b2)
    lpout <- lpSolve::lp('min', f.obj, f.con.0, f.dir, f.rhs)
    lpout$solution[1:p] - lpout$solution[-(1:p)]
  }
  parallel::stopCluster(cl)
  DD.2 <- make.symmetric(DD.2, symmetric)
  if(do.correct){
    tmp <- gen.inverse(GG)
    ind <- which(diag(DD.2) == 0)
    diag(DD.2)[ind] <- tmp[ind]
  }
  out <- list(DD = DD.2, eta = eta, symmetric = symmetric)
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

#' @keywords internal
correct.diag <- function(GG, DD){

  p <- dim(GG)[1]
  tmp <- gen.inverse(GG)
  ind <- which(diag(DD) <= 0)
  diag(DD)[ind] <- tmp[ind]

  # ind0 <- setdiff(1:p, ind)
  # mat <- t(t(DD[ind0, ind0])/sqrt(diag(DD)[ind0]))/sqrt(diag(DD)[ind0])
  # ind <- c(ind, ind0[apply(abs(mat), 1, max) > 1])
  # ind0 <- setdiff(1:p, ind)
  # if(length(ind) > 0){
  #   mat <- t(t(DD[ind, ind])/sqrt(diag(DD)[ind]))/sqrt(diag(DD)[ind])
  #   while(max(abs(mat)) - 1 > 1e-10){
  #     for(ii in ind[which(apply(abs(mat), 1, max) > 1)]){
  #       ind1 <- setdiff(ind, ii)
  #       DD[ii, ii] <- max(tmp[ii], (DD[ii, ind1]/sqrt(diag(DD)[ind1]))^2)
  #     }
  #     mat <- t(t(DD[ind, ind])/sqrt(diag(DD)[ind]))/sqrt(diag(DD)[ind])
  #   }
  # }

  return(DD)

}
