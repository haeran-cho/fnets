#' @title \code{l1}-regularised Yule-Walker estimation for VAR processes
#' @description Estimates the VAR parameter matrices via \code{l1}-regularised Yule-Walker estimation
#' and innovation covariance matrix via constrained \code{l1}-minimisation.
#' @details Further information can be found in Barigozzi, Cho and Owens (2021).
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param method a string specifying the method to be adopted for VAR process estimation; possible values are:
#' \itemize{
#'    \item{\code{"lasso"}}{ Lasso-type \code{l1}-regularised \code{M}-estimation}
#'    \item{\code{"ds"}}{ Dantzig Selector-type constrained \code{l1}-minimisation}
#' }
#' @param var.order order of the VAR process; if a vector of integers is supplied, the order is chosen via cross validation
#' @param lambda regularisation parameter; if \code{lambda = NULL}, cross validation is employed to select the parameter
#' @param cv.args a list specifying arguments for the cross validation procedure
#' for selecting the regularisation parameter (and VAR order). It contains:
#' \itemize{
#'    \item{\code{n.folds}}{ number of folds}
#'    \item{\code{path.length}}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{\code{do.plot}}{ whether to plot the output of the cross validation step}
#' }
#' @param n.iter maximum number of descent steps; applicable when \code{method = "lasso"}
#' @param tol numerical tolerance for increases in the loss function; applicable when \code{method = "lasso"}
#' @param n.cores number of cores to use for parallel computing, see \link[parallel]{makePSOCKcluster}; applicable when \code{method = "ds"}
#' @return a list which contains the following fields:
#' \item{beta}{ estimate of VAR parameter matrix; each column contains parameter estimates for the regression model for a given variable}
#' \item{Gamma}{ estimate of the innovation covariance matrix}
#' \item{lambda}{ regularisation parameter}
#' \item{convergence}{ returned when \code{method = "lasso"}; indicates whether a convergence criterion is met}
#' \item{var.order}{ VAR order}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#' @example R/examples/var_ex.R
#' @importFrom parallel detectCores
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @export
fnets.var  <- function(x, center = TRUE, method = c('lasso', 'ds'),
                       lambda = NULL, var.order = 1,
                       cv.args = list(n.folds = 1, path.length = 10, do.plot = FALSE),
                       n.iter = 100, tol = 1e-5, n.cores = min(parallel::detectCores() - 1, 3)){
  p <- dim(x)[1]
  n <- dim(x)[2]

  method <- match.arg(method, c('lasso', 'ds'))
  if(center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x
  dpca <- dyn.pca(xx, q = 0)
  acv <- dpca$acv

  icv <- yw.cv(xx, method = method,
               lambda.max = NULL, var.order = var.order,
               n.folds = cv.args$n.folds, path.length = cv.args$path.length,
               q = 0, kern.const = 4, do.plot = cv.args$do.plot)
  mg <- make.gg(acv$Gamma_i, icv$var.order)
  gg <- mg$gg; GG <- mg$GG

  if(method == 'lasso') ive <- var.lasso(GG, gg, lambda = icv$lambda, symmetric = 'min', n.iter = n.iter, tol = tol)
  if(method == 'ds') ive <- var.dantzig(GG, gg, lambda = icv$lambda, symmetric = 'min', n.cores = n.cores)
  ive$var.order <- icv$var.order
  ive$mean.x <- mean.x

  return(ive)

}

#' @title Lasso-type estimator of VAR processes via \code{l1}-regularised \code{M}-estimation
#' @keywords internal
var.lasso <- function(GG, gg, lambda, symmetric = 'min', n.iter = 100, tol = 1e-5){

  backtracking <- TRUE
  p <- ncol(gg)
  d <- nrow(gg)/ncol(gg)

  ii <- 0
  t.new <- t <- 1
  x <- gg * 0
  x.new <- y <- x
  diff.val <- tol - 1

  if(backtracking){
    L <- norm(GG, "F") / 5
    gamma <- 2
  } else L <- norm(GG, "F")

  obj.val <- rel.err <- c()
  while(ii < n.iter & diff.val < tol){
    ii <- ii + 1
    if(backtracking){
      L.bar <- L
      found <- FALSE
      while(!found){
        prox <- prox.func(y, lambda, L = 2 * L.bar, GG, gg)
        if(f.func(GG, gg, prox) <= Q.func(prox, y, L.bar, GG, gg)){
          found <- TRUE
        } else{
          L.bar <- L.bar * gamma
        }
      }
      L <- L.bar
    } else prox <- prox.func(y, lambda, L = 2 * L, GG, gg)

    x <- x.new
    x.new <- prox
    t <- t.new
    t.new <- (1 + sqrt(1 + 4*t^2))/2
    y <- x.new + (t - 1) / t.new * (x.new - y)

    obj.val <- c(obj.val, f.func(GG, gg, x.new) + lambda * sum(abs(x.new)))
    if(ii > 1) diff.val <- obj.val[ii] - obj.val[ii - 1]
  }
  if(ii == n.iter) warning("lasso estimation did not converge")

  A <- t(x.new)
  Gamma <- GG[1:p, 1:p]
  for(ll in 1:d) Gamma <- Gamma - A[, (ll - 1) * p + 1:p] %*% gg[(ll - 1) * p + 1:p, ]
  Gamma <- make.symmetric(Gamma, symmetric)
  out <- list(beta = x.new, Gamma = Gamma, lambda = lambda, convergence = (abs(diff.val) <= abs(obj.val[1]) * tol))

  return(out)

}

#' @title Dantzig selector-type estimator of VAR processes via constrained \code{l1}-minimisation
#' @importFrom parallel makePSOCKcluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lpSolve lp
#' @keywords internal
var.dantzig <- function(GG, gg, lambda, symmetric = 'min', n.cores = min(parallel::detectCores() - 1, 3)){

  p <- dim(gg)[2]
  d <- dim(gg)[1]/dim(gg)[2]
  beta <- gg * 0

  f.obj <- rep(1, 2 * p * d)
  f.con <- rbind(-GG, GG)
  f.con <- cbind(f.con,-f.con)
  f.dir <- rep('<=', 2 * p * d)

  cl <- parallel::makePSOCKcluster(n.cores)
  doParallel::registerDoParallel(cl)

  ii <- 1
  beta <- foreach::foreach(ii = 1:p, .combine = 'cbind', .multicombine = TRUE, .export = c('lp')) %dopar% {
    b1 <- rep(lambda, p * d) - gg[, ii]
    b2 <- rep(lambda, p * d) + gg[, ii]
    f.rhs <- c(b1, b2)
    lpout <- lpSolve::lp('min', f.obj, f.con, f.dir, f.rhs)
    lpout$solution[1:(p * d)] - lpout$solution[-(1:(p * d))]
  }
  parallel::stopCluster(cl)

  A <- t(beta)
  Gamma <- GG[1:p, 1:p]
  for(ll in 1:d) Gamma <- Gamma - A[, (ll - 1) * p + 1:p] %*% gg[(ll - 1) * p + 1:p, ]
  Gamma <- make.symmetric(Gamma, symmetric)

  out <- list(beta = beta, Gamma = Gamma, lambda = lambda)
  return(out)

}

#' @title Cross validation for factor-adjusted VAR estimation
#' @importFrom graphics abline legend matplot
#' @keywords internal
yw.cv <- function(xx, method = c('lasso', 'ds'),
                  lambda.max = NULL, var.order = 1,
                  n.folds = 1, path.length = 10,
                  q = 0, kern.const = 4, do.plot = FALSE){

  n <- ncol(xx)
  p <- nrow(xx)

  if(is.null(lambda.max)) lambda.max <- max(abs(xx %*% t(xx)/n)) * 1
  lambda.path <- round(exp(seq(log(lambda.max), log(lambda.max * .0001), length.out = path.length)), digits = 10)

  cv.err.mat <- matrix(0, nrow = path.length, ncol = length(var.order))
  dimnames(cv.err.mat)[[1]] <- lambda.path
  dimnames(cv.err.mat)[[2]] <- var.order
  ind.list <- split(1:n, ceiling(n.folds*(1:n)/n))

  for(fold in 1:n.folds){
    train.ind <- 1:ceiling(length(ind.list[[fold]]) * .5)
    train.x <- xx[, ind.list[[fold]][train.ind]]
    test.x  <- xx[, ind.list[[fold]][- train.ind]]
    train.acv <- dyn.pca(train.x, q = q, kern.const = kern.const, mm = max(var.order))$acv$Gamma_i
    test.acv <- dyn.pca(test.x, q = q, kern.const = kern.const, mm = max(var.order))$acv$Gamma_i

    for(jj in 1:length(var.order)){
      mg <- make.gg(train.acv, var.order[jj])
      gg <- mg$gg; GG <- mg$GG
      mg <- make.gg(test.acv, var.order[jj])
      test.gg <- mg$gg; test.GG <- mg$GG
      for(ii in 1:path.length){
        if(method == 'ds') train.beta <- var.dantzig(GG, gg, lambda = lambda.path[ii])$beta
        if(method == 'lasso') train.beta <- var.lasso(GG, gg, lambda = lambda.path[ii])$beta
        beta.gg <- t(train.beta) %*% test.gg
        cv.err.mat[ii, jj] <- cv.err.mat[ii, jj] +
          sum(diag(test.acv[,, 1] - beta.gg - t(beta.gg) + t(train.beta) %*% test.GG %*% (train.beta) ))
      }
    }
  }
  cv.err.mat[cv.err.mat < 0] <- Inf
  lambda.min <- min(lambda.path[apply(cv.err.mat, 1, min) == min(apply(cv.err.mat, 1, min))])
  order.min <- min(var.order[apply(cv.err.mat, 2, min) == min(apply(cv.err.mat, 2, min))])

  if(do.plot){
    matplot(lambda.path, cv.err.mat, type = 'b', col = 2:(max(var.order) + 1), pch = 2:(max(var.order) + 1),
            log = 'x', xlab = 'lambda (log scale)', ylab = 'CV error', main = 'CV for VAR parameter estimation')
    abline(v = lambda.min)
    legend('topleft', legend = var.order, col = 2:(max(var.order) + 1), pch = 2:(max(var.order) + 1), lty = 1)
  }

  out <- list(lambda = lambda.min, var.order = order.min, cv.error = cv.err.mat, lambda.path = lambda.path)
  return(out)

}

#' @title Forecasting idiosyncratic VAR process
#' @description Produces forecasts of the idiosyncratic VAR process
#' for a given forecasting horizon by estimating the best linear predictors
#' @param object \code{fnets} object
#' @param x input time series matrix, with each row representing a variable
#' @param cpre output of \link[fnets]{common.predict}
#' @param h forecast horizon
#' @return a list containing
#' \item{is}{ in-sample estimator of the idiosyncratic component}
#' \item{fc}{ forecasts of the idiosyncratic component for a given forecasting horizon \code{h}}
#' \item{h}{ forecast horizon}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @examples
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.common1(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#' out <- fnets(x, q = NULL, idio.var.order = 1, idio.method = "lasso", lrpc.method = "none")
#' cpre <- common.predict(out, x, h = 1, common.method = 'restricted', r = NULL)
#' ipre <- idio.predict(out, x, cpre, h = 1)
#' @export
idio.predict <- function(object, x, cpre, h = 1){

  p <- dim(x)[1]; n <- dim(x)[2]
  xx <- x - object$mean.x
  beta <- object$idio.var$beta
  d <- dim(beta)[1]/p
  A <- t(beta)

  is <- xx - cpre$is
  if(h >= 1){
    fc <- matrix(0, nrow = p, ncol = h)
    for(ii in 1:h){
      for(ll in 1:d) fc[, ii] <- fc[, ii] + A[, p * (ll - 1) + 1:p] %*% is[, n + ii - ll]
      is <- cbind(is, fc[, ii])
    }
  } else fc <- NA

  out <- list(is = is[, 1:n], fc = fc, h = h)
  return(out)

}

#' @keywords internal
make.gg <- function(acv, d){

  p <- dim(acv)[1]
  gg <- matrix(0, nrow = p * d, ncol = p)
  GG <- matrix(0, p * d, p * d)
  for(ll in 1:d){
    gg[(ll - 1) * p + 1:p, ] <- acv[,, ll + 1]
    for(lll in ll:d){
      GG[(ll - 1) * p + 1:p, (lll - 1) * p + 1:p] <- t(acv[,, 1 + lll - ll])
      GG[(lll - 1) * p + 1:p, (ll - 1) * p + 1:p] <- acv[,, 1 + lll - ll]
    }
  }
  out <- list(gg = gg, GG = GG)
  return(out)

}

#' @keywords internal
f.func <- function(GG, gg, A){
  return(.5 * sum(diag(t(A) %*% GG %*% A - 2 * t(A) %*% gg)))
}

#' @keywords internal
gradf.func <- function(GG, gg, A){
  return(GG %*% (A) - gg)
}

#' @keywords internal
Q.func <- function(A, A.up, L, GG, gg){
  Adiff <- (A - A.up)
  return(f.func(GG, gg, A.up) + sum(Adiff * gradf.func(GG, gg, A.up)) + 0.5 * L * norm(Adiff, "F")^2)
}

#' @keywords internal
prox.func <- function(B, lambda, L, GG, gg){
  b <- B - (1/L) * gradf.func(GG, gg, B)
  sgn <- sign(b)
  ab <- abs(b)
  sub <- ab - 2*lambda/L
  sub[sub < 0] <- 0
  out <- sub * sgn
  return(as.matrix(out))
}
