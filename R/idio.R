#library(lpSolve)
#library(foreach)


# var.fit <- function(){
#   ## idio estimation
#   mg <- make.gg(acv$Gamma_i, idio.var.order)
#   gg <- mg$gg; GG <- mg$GG
#
#   icv <- idio.cv(xx, lambda.max = max(abs(GG)), var.order = idio.var.order, idio.method = idio.method,
#                  path.length = idio.cv.args$path.length, n.folds = idio.cv.args$n.folds,
#                  q = q, kern.bandwidth.const = kern.bandwidth.const, cv.plot = idio.cv.args$n.folds)
#   if(idio.method == 'lasso') ive <- var.lasso(GG, gg, lambda = icv$lambda, symmetric = idio.cv.args$symmetric)
#   if(idio.method == 'ds') ive <- var.dantzig(GG, gg, lambda = icv$lambda, symmetric = idio.cv.args$symmetric)
# }

#' @title Lasso-type Yule-Walker estimation for VAR processes
#' @description Returns parameter estimates for the idiosyncratic VAR and the corresponding Gamma matrix
#' @details
#'  Further information can be found in Barigozzi, Cho and Owens (2021).
#'
#' @param GG,gg output from `make.gg`
#' @param lambda regularisation parameter
#' @param symmetric type of symmetry to enforce on `Gamma`, one of 'min', 'max', 'avg', 'none'
#' @param niter maximum number of descent steps
#' @param tol numerical tolerance for increases in the loss function
#' @param do.plot return a plot of the loss function against descent steps
#' @return A list which contains the following fields:
#' \itemize{
#' \item{beta: VAR parameters}
#' \item{lambda: regularisation parameter}
#' \item{Gamma: Estimated noise covariance}
#' \item{loss: Objective function value}
#' }
#' @example examples/idio.R
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @export
var.lasso <- function(GG, gg, lambda, symmetric = 'min', niter = 100, tol = 0, do.plot = FALSE){
  backtracking = TRUE
  p <- ncol(gg)
  d <- nrow(gg)/ncol(gg)

  GtG <- t(GG) %*% GG
  Gtg <- t(GG) %*% gg

  ii<-0
  tnew <- t <- 1
  beta1 <- gg*0
  beta.mid <-  beta.up <- prox <- gg*0#matrix(0, nrow = p, ncol = p)
  diff <- tol-1

  if(backtracking){
    L <- norm(GG, "F")^2 / 5
    gamma <- 2
  }else{
    L <- norm(GG, "F")^2
  }

  obj.val = rel.err <- c()
  while(ii < niter & diff < tol){
    ii <- ii+1
    if(backtracking){
      L.bar <- L
      found <- FALSE
      while(!found){
        # prox1 <- prox.sparse.func(y1, y, A, b, 2*L.bar, lambda, AtA, Atb)
        # GG %*% x - Atb
        prox <- fnsl.update(beta.up, beta1, lambda, eta =2*L.bar, GtG, Gtg)
        if(f.func(GG,gg,prox) <= Q.func(prox,beta.up,L.bar,GG,gg, GtG, Gtg)){
          found <- TRUE
        }else{
          L.bar <- L.bar * gamma
        }
      }
      L <- L.bar
    } else prox <- fnsl.update(beta.up, beta1, lambda, eta =2*L, GtG,gtg)
    beta1 <- beta.mid
    beta.mid <- prox
    t = tnew
    tnew <- (1 + sqrt(1 + 4*t^2))/2
    beta.up <- beta.mid + (t - 1) / tnew * (beta.mid - beta1)

    obj.val <- c(obj.val, f.func(GG, gg, beta.mid) + lambda*sum(abs(beta.mid)))
    #if(obj.val[ii]==min(obj.val)) beta.min <- beta.mid #pick minimiser
    if(ii>1) diff <- obj.val[ii] - obj.val[ii-1] #max(abs(beta.mid - beta1))
  }

  A <- t(beta.mid)
  Gamma <- GG[1:p, 1:p]
  for(ll in 1:d) Gamma <- Gamma - A[,(ll - 1) * p + 1:p ] %*% gg[(ll - 1) * p + 1:p, ]
  Gamma <- make.symmetric(Gamma, symmetric)

  if(do.plot) plot.ts(obj.val);# abline(v = which.min(obj.val), col = "red")
  out <- list(beta = beta.mid, lambda = lambda, Gamma = Gamma, loss = obj.val)
  return(out)
}

# #' @export
#' @description internal function
#' @keywords internal
f.func <- function(GG, gg, A){
return(0.5 * norm(   (GG %*% (A) - gg ) , "F")^2) ##
}
# #' @export
#' @description internal function
#' @keywords internal
gradf.func <- function(GtG, Gtg, A){
  return( (GtG %*% (A) - Gtg ) ) # GG %*% here??
}
# #' @export
#' @description internal function
#' @keywords internal
Q.func <- function(A, A.up, L, GG, gg, GtG, Gtg){
  Adiff <- (A - A.up)
  return(f.func(GG, gg, A.up) + sum(Adiff * gradf.func(GtG, Gtg,A.up)) + 0.5 * L * norm(Adiff, "F")^2 )
}
# #' @export
#' @description internal function
#' @keywords internal
fnsl.update <- function(B, B_md, lambda, eta, GtG, Gtg){
  b <- B - (1/eta) * (GtG %*% (B_md) - Gtg ) ##
  sgn <- sign(b)
  ab <- abs(b)
  sub <- ab - 2*lambda/eta
  sub[sub < 0] <- 0
  out <- sub * sgn
  return(as.matrix(out))
}


#' @title Dantzig selector-type Yule-Walker estimation for VAR processes
#' @description Returns parameter estimates for the idiosyncratic VAR and the corresponding Gamma matrix
#' @details
#'  Further information can be found in Barigozzi, Cho and Owens (2021).
#'
#' @param GG,gg output from `make.gg`
#' @param lambda regularisation parameter
#' @param symmetric type of symmetry to enforce on `Gamma`, one of 'min', 'max', 'avg', 'none'
#' @param n.cores number of cores to use for parallel computing
#' @return A list which contains the following fields:
#' \itemize{
#' \item{beta: VAR parameters}
#' \item{lambda: regularisation parameter}
#' \item{Gamma: Estimated noise covariance}
#' }
#' @example examples/idio.R
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @export
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

  beta <- foreach::foreach(ii = 1:p, .combine = 'cbind', .multicombine = TRUE, .export = c('lp')) %dopar% {
    b1 <- rep(lambda, p * d) - gg[, ii]
    b2 <- rep(lambda, p * d) + gg[, ii]
    f.rhs <- c(b1, b2)
    lpout <- lp('min', f.obj, f.con, f.dir, f.rhs)
    lpout$solution[1:(p * d)] - lpout$solution[-(1:(p * d))]
  }
  parallel::stopCluster(cl)

  A <- t(beta)
  Gamma <- GG[1:p, 1:p]
  for(ll in 1:d) Gamma <- Gamma - A[, (ll - 1) * p + 1:p] %*% gg[(ll - 1) * p + 1:p, ]
  Gamma <- make.symmetric(Gamma, symmetric)

  out <- list(beta = (beta), lambda = lambda, Gamma = Gamma)
  return(out)
}

####

#' @title Cross-validation for l1-regularised VAR estimation
#' @description selects the prediction-optimal regularisation parameter for the estimation of the idiosyncratic VAR
#' @details
#'  Further information can be found in Barigozzi, Cho and Owens (2021).
#'
#' @param xx centred input time series matrix, with each row representing a time series
#' @param lambda.max maximum regularisation parameter, if NULL this is set to the smallest which sets all entries to 0
#' @param var.order vector of VAR orders to consider
#' @param idio.method estimation method, one of "lasso" or "ds"
#' @param path.length number of regularisation parameters to consider
#' @param n.folds number of CV folds
#' @param q factor number
#' @param kern.bandwidth.const constant to determine bandwidth size
#' @param cv.plot return a plot of the CV error against regularisation parameters, stratified by VAR order
#' @return A list which contains the following fields:
#' \itemize{
#' \item{\code{'lambda'}}{ minimising argument}
#' \item{\code{'var.order'}}{ minimising order}
#' \item{\code{'cv.error'}}{ matrix of errors}
#' \item{\code{'lambda.path'}}{ candidate lambda values}
#' }
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @examples
#' idio.cv(sample.data, idio.method = "lasso", q=2)
#' @export
idio.cv <- function(xx, lambda.max = NULL, var.order = 1, idio.method = c('lasso', 'ds'),
                    path.length = 10, n.folds = 1,
                    q = 0, kern.bandwidth.const = 4, cv.plot = TRUE){

  n <- ncol(xx)
  p <- nrow(xx)

  if(is.null(lambda.max)) lambda.max <- max(abs(xx %*% t(xx)/n)) * 1
  lambda.path <- round(exp(seq(log(lambda.max), log(lambda.max * .0001), length.out = path.length)), digits = 10)

  cv.err.mat <- matrix(0, nrow = path.length, ncol = length(var.order))
  ind.list <- split(1:n, ceiling(n.folds*(1:n)/n))

  for(fold in 1:n.folds){
    train.ind <- 1:ceiling(length(ind.list[[fold]]) * .5)
    train.x <- xx[, ind.list[[fold]][train.ind]]
    test.x  <- xx[, ind.list[[fold]][- train.ind]]
    train.acv <- dyn.pca(train.x, q = q, kern.bandwidth.const = kern.bandwidth.const)$acv$Gamma_i
    test.acv <- dyn.pca(test.x, q = q, kern.bandwidth.const = kern.bandwidth.const)$acv$Gamma_i

    for(jj in 1:length(var.order)){
      mg <- make.gg(train.acv, var.order[jj])
      gg <- mg$gg; GG <- mg$GG
      mg <- make.gg(test.acv, var.order[jj])
      test.gg <- mg$gg; test.GG <- mg$GG
      for(ii in 1:path.length){
        if(idio.method == 'ds') train.beta <- var.dantzig(GG, gg, lambda = lambda.path[ii])$beta
        if(idio.method == 'lasso') train.beta <- var.lasso(GG, gg, lambda = lambda.path[ii])$beta
        # train.beta[1:5,1:5]
        beta.gg <- t(train.beta) %*% test.gg
        cv.err.mat[ii, jj] <- cv.err.mat[ii, jj] + sum(diag(test.acv[,, 1] - beta.gg - t(beta.gg) + t(train.beta) %*% test.GG %*% (train.beta) )) #test.GG
      }
    }
  }

  cv.err.mat
  lambda.min <- lambda.path[which.min(apply(cv.err.mat, 1, min))]
  order.min <- var.order[which.min(apply(cv.err.mat, 2, min))]

  if(cv.plot){
  matplot(lambda.path, cv.err.mat, type = 'b', col = 2:( max(var.order) + 1), pch = 2:(max(var.order) + 1), log = 'x', xlab = 'lambda (log scale)', ylab = 'CV error')
  abline(v = lambda.min)
  legend('topleft', legend = var.order, col = 2:(max(var.order) + 1), pch = 2:(max(var.order) + 1), lty = 1)
  }


  out <- list(lambda = lambda.min, var.order = order.min,  cv.error = cv.err.mat, lambda.path = lambda.path)
  return(out)

}

####

#' @title Prediction for the idiosyncratic VAR process
#' @description Predictsidiosyncratic components from a fnets object for new data
#' @param object fnets object
#' @param x input time series matrix, with each row representing a time series
#' @param cpre estimated common component
#' @param h forecast horizon
#' @return A list containing
#' \itemize{
#' \item{\code{'is'}}{ in-sample estimation}
#' \item{\code{'fc'}}{ forecast}
#' \item{\code{'h'}}{ forecast horizon}
#' }
#' @example examples/predict.R
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @export
idio.predict <- function(object, x, cpre, h = 1){

  xx <- x - object$mean.x
  p <- dim(x)[1]; n <- dim(x)[2]
  beta <- object$idio.var$beta
  d <- dim(beta)[1]/p
  A <- t(beta)

  is <- xx - cpre$is
  if(h >= 1){
    fc <- matrix(0, nrow = p, ncol = h)
    for(ii in 1:h){
      for(ll in 1:d) fc[, ii] <- fc[, ii] + A[, p * (ll - 1) + 1:p] %*% is[, n + ii - ll]
      is <- cbind(is, fc)
    }
  } else fc <- NA

  out <- list(is = is[, 1:n], fc = fc, h = h)
  return(out)
}

## misc


#' @title Construct `gg` and `GG` matrices for Yule-Walker estimation
#' @param acv autocovariance array
#' @param d order of VAR
#' @export
# #' @keywords internal
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



