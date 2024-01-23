#' @title \code{l1}-regularised Yule-Walker estimation for VAR processes
#' @description Estimates the VAR parameter matrices via \code{l1}-regularised Yule-Walker estimation
#' and innovation covariance matrix via constrained \code{l1}-minimisation.
#' @details Further information can be found in Barigozzi, Cho and Owens (2024+).
#' @param x input time series each column representing a time series variable; it is coerced into a \link[stats]{ts} object
#' @param center whether to de-mean the input \code{x}
#' @param method a string specifying the method to be adopted for VAR process estimation; possible values are:
#' \describe{
#'    \item{\code{"lasso"}}{ Lasso-type \code{l1}-regularised \code{M}-estimation}
#'    \item{\code{"ds"}}{ Dantzig Selector-type constrained \code{l1}-minimisation}
#' }
#' @param var.order order of the VAR process; if a vector of integers is supplied, the order is chosen via \code{tuning}
#' @param lambda \code{l1}-regularisation parameter; if \code{lambda = NULL}, \code{tuning} is employed to select the parameter
#' @param tuning.args a list specifying arguments for \code{tuning}
#' for selecting the regularisation parameter (and VAR order). It contains:
#' \describe{
#'    \item{\code{tuning}}{a string specifying the selection procedure for \code{var.order} and \code{lambda}; possible values are:
#'    \code{"cv"} for cross validation, and \code{"bic"} for information criterion}
#'    \item{\code{n.folds}}{ if \code{tuning = "cv"}, positive integer number of folds}
#'    \item{\code{penalty}}{ if \code{tuning = "bic"}, penalty multiplier between 0 and 1; if \code{penalty = NULL}, it is set to \code{1/(1+exp(dim(x)[1])/dim(x)[2]))}} by default
#'    \item{\code{path.length}}{ positive integer number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#' }
#' @param do.threshold whether to perform adaptive thresholding of VAR parameter estimator with \link[fnets]{threshold}
#' @param n.iter maximum number of descent steps, by default depends on \code{var.order}; applicable when \code{method = "lasso"}
#' @param tol numerical tolerance for increases in the loss function; applicable when \code{method = "lasso"}
#' @param n.cores number of cores to use for parallel computing, see \link[parallel]{makePSOCKcluster}; applicable when \code{method = "ds"}
#' @return a list which contains the following fields:
#' \item{beta}{ estimate of VAR parameter matrix; each column contains parameter estimates for the regression model for a given variable}
#' \item{Gamma}{ estimate of the innovation covariance matrix}
#' \item{lambda}{ \code{l1}-regularisation parameter}
#' \item{var.order}{ VAR order}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#' @example R/examples/var_ex.R
#' @importFrom parallel detectCores
#' @importFrom stats as.ts
#' @references Barigozzi, M., Cho, H. & Owens, D. (2024+) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. Journal of Business & Economic Statistics (to appear).
#' @references Owens, D., Cho, H. & Barigozzi, M. (2024+) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. The R Journal (to appear).
#' @export
fnets.var <- function(x,
                      center = TRUE,
                      method = c("lasso", "ds"),
                      lambda = NULL,
                      var.order = 1,
                      tuning.args = list(
                        tuning = c("cv", "bic"),
                        n.folds = 1,
                        penalty = NULL,
                        path.length = 10
                      ),
                      do.threshold = FALSE,
                      n.iter = NULL,
                      tol = 0,
                      n.cores = 1) {
  x <- t(as.ts(x))
  p <- dim(x)[1]
  n <- dim(x)[2]

  if(!is.null(lambda)) lambda <- max(0, tol)
  if(!is.null(n.iter)) n.iter <- posint(n.iter)
  tol <- max(0, tol)
  n.cores <- posint(n.cores)
  tuning.args <- check.list.arg(tuning.args)

  method <- match.arg(method, c("lasso", "ds"))
  tuning <- match.arg(tuning.args$tuning, c("cv", "bic"))

  args <- as.list(environment())
  args$x <- t(args$x)

  ifelse(center,mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
  xx <- x - mean.x
  dpca <- dyn.pca(xx, q = 0)
  acv <- dpca$acv

  ive <- fnets.var.internal(xx, acv, method = method,
                            lambda = lambda,
                            var.order = var.order,
                            tuning.args = tuning.args,
                            do.threshold = do.threshold,
                            n.iter = n.iter,
                            tol = tol,
                            n.cores = n.cores)
  ive$mean.x <- mean.x
  args$do.lrpc <- ive$do.lrpc <- FALSE
  ive$q <- 0

  attr(ive, "class") <- "fnets"
  attr(ive, "factor") <- "none"
  attr(ive, "args") <- args
  return(ive)
}

#' @title internal function for \code{fnets.var}
#' @keywords internal
fnets.var.internal <- function(xx,
                               acv,
                               method = c("lasso", "ds"),
                               lambda = NULL,
                               var.order = 1,
                               tuning.args = list(
                                 tuning = c("cv", "bic"),
                                 n.folds = 1,
                                 penalty = NULL,
                                 path.length = 10
                               ),
                               do.threshold = FALSE,
                               n.iter = NULL,
                               tol = 0,
                               n.cores = 1){
  for (ii in 1:length(var.order)) var.order[ii] <- posint(var.order[ii])

  method <- match.arg(method, c("lasso", "ds"))
  tuning <- match.arg(tuning.args$tuning, c("cv", "bic"))

  if(tuning == "cv") {
    icv <- yw.cv(
      xx,
      method = method,
      lambda.max = NULL,
      var.order = var.order,
      n.folds = tuning.args$n.folds,
      path.length = tuning.args$path.length,
      q = 0,
      kern.bw = NULL,
      n.cores = n.cores
    )
  }

  if(tuning == "bic") {
    icv <- yw.ic(
      xx,
      method = method,
      lambda.max = NULL,
      var.order = var.order,
      penalty = tuning.args$penalty,
      path.length = tuning.args$path.length,
      q = 0,
      kern.bw = NULL
    )
  }

  mg <- make.gg(acv$Gamma_i, icv$order.min)
  gg <- mg$gg
  GG <- mg$GG

  if(method == "lasso"){
    if(is.null(n.iter)) n.iter <- icv$order.min*100
    ive <-
      var.lasso(
        GG,
        gg,
        lambda = icv$lambda,
        symmetric = "min",
        n.iter = n.iter,
        tol = tol
      )
  }
  if(method == "ds")
    ive <-
    var.dantzig(
      GG,
      gg,
      lambda = icv$lambda,
      symmetric = "min",
      n.cores = n.cores
    )
  ive$var.order <- icv$order.min
  if(do.threshold)
    ive$beta <- threshold(ive$beta)$thr.mat
  attr(ive, "data") <- icv
  return(ive)
}

#' @title Lasso-type estimator of VAR processes via \code{l1}-regularised \code{M}-estimation
#' @keywords internal
var.lasso <-
  function(GG,
           gg,
           lambda,
           symmetric = "min",
           n.iter = 100,
           tol = 0) {
    backtracking <- TRUE
    p <- ncol(gg)
    d <- nrow(gg) / ncol(gg)

    ii <- 0
    t.new <- t <- 1
    x <- gg * 0
    x.new <- y <- x
    diff.val <- tol - 1

    if(backtracking) {
      L <- norm(GG, "F") / 5
      gamma <- 2
    } else {
      L <- norm(GG, "F")
    }

    obj.val <- rel.err <- c()
    while (ii < n.iter & diff.val < tol) {
      ii <- ii + 1
      if(backtracking) {
        L.bar <- L
        found <- FALSE
        while (!found) {
          prox <- prox.func(y, lambda, L = 2 * L.bar, GG, gg)
          ifelse(f.func(GG, gg, prox) <= Q.func(prox, y, L.bar, GG, gg),
                 found <- TRUE,
                 L.bar <- L.bar * gamma
                 )
        }
        L <- L.bar
      } else {
        prox <- prox.func(y, lambda, L = 2 * L, GG, gg)
      }

      x <- x.new
      x.new <- prox
      t <- t.new
      t.new <- (1 + sqrt(1 + 4 * t^2)) / 2
      y <- x.new + (t - 1) / t.new * (x.new - y)

      obj.val <-
        c(obj.val, f.func(GG, gg, x.new) + lambda * sum(abs(x.new)))
      if(ii > 1)
        diff.val <- obj.val[ii] - obj.val[ii - 1]
    }
    A <- t(x.new)
    Gamma <- GG[1:p, 1:p]
    for (ll in 1:d)
      Gamma <-
      Gamma - A[, (ll - 1) * p + 1:p] %*% gg[(ll - 1) * p + 1:p, ]
    Gamma <- make.symmetric(Gamma, symmetric)
    out <-
      list(
        beta = x.new,
        Gamma = Gamma,
        lambda = lambda
      )

    return(out)
  }

#' @title Dantzig selector-type estimator of VAR processes via constrained \code{l1}-minimisation
#' @importFrom parallel makePSOCKcluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lpSolve lp
#' @keywords internal
var.dantzig <-
  function(GG,
           gg,
           lambda,
           symmetric = "min",
           n.cores = 1) {
    p <- dim(gg)[2]
    d <- dim(gg)[1] / dim(gg)[2]
    beta <- gg * 0

    f.obj <- rep(1, 2 * p * d)
    f.con <- rbind(-GG, GG)
    f.con <- cbind(f.con, -f.con)
    f.dir <- rep("<=", 2 * p * d)

    cl <- parallel::makePSOCKcluster(n.cores)
    doParallel::registerDoParallel(cl)

    ii <- 1
    beta <-
      foreach::foreach(
        ii = 1:p,
        .combine = "cbind",
        .multicombine = TRUE,
        .export = c("lp")
      ) %dopar% {
        b1 <- rep(lambda, p * d) - gg[, ii]
        b2 <- rep(lambda, p * d) + gg[, ii]
        f.rhs <- c(b1, b2)
        lpout <- lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
        lpout$solution[1:(p * d)] - lpout$solution[-(1:(p * d))]
      }
    parallel::stopCluster(cl)

    A <- t(beta)
    Gamma <- GG[1:p, 1:p]
    for (ll in 1:d)
      Gamma <-
      Gamma - A[, (ll - 1) * p + 1:p] %*% gg[(ll - 1) * p + 1:p, ]
    Gamma <- make.symmetric(Gamma, symmetric)

    out <- list(beta = beta,
                Gamma = Gamma,
                lambda = lambda)
    return(out)
  }

#' @title Cross validation for factor-adjusted VAR estimation
#' @keywords internal
yw.cv <- function(xx,
                  method = c("lasso", "ds"),
                  lambda.max = NULL,
                  var.order = 1,
                  n.folds = 1,
                  path.length = 10,
                  q = 0,
                  kern.bw = NULL,
                  n.cores = 1) {
  n <- ncol(xx)
  p <- nrow(xx)

  if(is.null(kern.bw))
    kern.bw <- 4 * floor((n / log(n))^(1 / 3))
  if(is.null(lambda.max))
    lambda.max <- max(abs(xx %*% t(xx) / n)) * 1
  lambda.path <-
    round(exp(seq(
      log(lambda.max), log(lambda.max * .0001), length.out = path.length
    )), digits = 10)

  cv.err.mat <-
    matrix(0, nrow = path.length, ncol = length(var.order))
  dimnames(cv.err.mat)[[1]] <- lambda.path
  dimnames(cv.err.mat)[[2]] <- var.order
  ind.list <- split(1:n, ceiling(n.folds * (1:n) / n))

  for (fold in 1:n.folds) {
    train.ind <- 1:ceiling(length(ind.list[[fold]]) * .5)
    train.x <- xx[, ind.list[[fold]][train.ind]]
    test.x <- xx[, ind.list[[fold]][-train.ind]]
    train.acv <-
      dyn.pca(train.x,
              q = q,
              kern.bw = kern.bw,
              mm = max(var.order))$acv$Gamma_i
    test.acv <-
      dyn.pca(test.x,
              q = q,
              kern.bw = kern.bw,
              mm = max(var.order))$acv$Gamma_i

    for (jj in 1:length(var.order)) {
      mg <- make.gg(train.acv, var.order[jj])
      gg <- mg$gg
      GG <- mg$GG
      mg <- make.gg(test.acv, var.order[jj])
      test.gg <- mg$gg
      test.GG <- mg$GG
      for (ii in 1:path.length) {
        if(method == "ds")
          train.beta <-
            var.dantzig(GG, gg, lambda = lambda.path[ii], n.cores = n.cores)$beta
        if(method == "lasso")
          train.beta <-
            var.lasso(GG, gg, lambda = lambda.path[ii])$beta
        beta.gg <- t(train.beta) %*% test.gg
        cv.err.mat[ii, jj] <- cv.err.mat[ii, jj] +
          sum(diag(
            test.acv[, , 1] - beta.gg - t(beta.gg) + t(train.beta) %*% test.GG %*% (train.beta)
          ))
      }
    }
  }
  cv.err.mat[cv.err.mat < 0] <- Inf
  lambda.min <-
    min(lambda.path[apply(cv.err.mat, 1, min) == min(apply(cv.err.mat, 1, min))])
  order.min <-
    min(var.order[apply(cv.err.mat, 2, min) == min(apply(cv.err.mat, 2, min))])



  out <-
    list(
      lambda = lambda.min,
      order.min = order.min,
      error = cv.err.mat,
      lambda.path = lambda.path,
      var.order = var.order
    )
  return(out)
}


# ebic

#' @title logarithmic factorial of `n`
#' @keywords internal
logfactorial <- function(n)
  sum(log(1:max(n, 1)))

#' @title full likelihood
#' @keywords internal
f.func.full <- function(GG, gg, A) {
  return(0.5 * sum(diag(
    GG[1:ncol(A), 1:ncol(A)] + t(A) %*% GG %*% A - t(A) %*% gg - t(gg) %*% A
  )))
}

f.func.mat <- function(GG, gg, A) {
  return(0.5 * (diag(
    GG[1:ncol(A), 1:ncol(A)] + t(A) %*% GG %*% A - t(A) %*% gg - t(gg) %*% A
  )))
}

#' @title extended Bayesian Information Criterion
#' @keywords internal
ebic <- function(object, n, penalty = 0) {
  penalty <- max( min(1,penalty) , 0)
  beta <- object$idio.var$beta
  p <- ncol(beta)
  d <- nrow(beta) / ncol(beta)
  sparsity <- sum(beta != 0)
  mg <- make.gg(object$acv$Gamma_i, d)
  gg <- mg$gg
  GG <- mg$GG
  n / 2 * log(2 * f.func.full(GG, gg, beta)) + sparsity * log(n) +
    2 * penalty * (logfactorial(p^2 * d) - logfactorial(sparsity) - logfactorial(p ^
                                                                                        2 * d - sparsity))
}


#' @title Information criterion for factor-adjusted VAR estimation
#' @keywords internal
yw.ic <- function(xx,
                  method = c("lasso", "ds"),
                  lambda.max = NULL,
                  var.order = 1,
                  penalty = NULL,
                  path.length = 10,
                  q = 0,
                  kern.bw = NULL) {
  n <- ncol(xx)
  p <- nrow(xx)
  if(is.null(kern.bw))
    kern.bw <- 4 * floor((n / log(n))^(1 / 3))
  if(is.null(penalty))
    penalty <- 1 / (1 + exp(5 - p / n))
  if(is.null(lambda.max))
    lambda.max <- max(abs(xx %*% t(xx) / n)) * 1
  lambda.path <-
    round(exp(seq(
      log(lambda.max), log(lambda.max * .0001), length.out = path.length
    )), digits = 10)

  ic.err.mat <-
    matrix(0, nrow = path.length, ncol = length(var.order))
  dimnames(ic.err.mat)[[1]] <- lambda.path
  dimnames(ic.err.mat)[[2]] <- var.order


  acv <-
    dyn.pca(xx,
            q = q,
            kern.bw = kern.bw,
            mm = max(var.order))$acv$Gamma_i

  for (jj in 1:length(var.order)) {
    mg <- make.gg(acv, var.order[jj])
    gg <- mg$gg
    GG <- mg$GG
    test.gg <- mg$gg
    test.GG <- mg$GG
    for (ii in 1:path.length) {
      if(method == "ds")
        beta <- var.dantzig(GG, gg, lambda = lambda.path[ii])$beta
      if(method == "lasso")
        beta <- var.lasso(GG, gg, lambda = lambda.path[ii])$beta
      beta <- threshold(beta)$thr.mat
      sparsity <- sum(beta[, 1] != 0)
      ifelse(sparsity == 0,pen <- 0,
             pen <-2 * penalty * (logfactorial(var.order[jj] * p) - log(sparsity) - log(var.order[jj] * p - sparsity)))

      ic.err.mat[ii, jj] <-
        ic.err.mat[ii, jj] + n / 2 * log(2 * f.func.mat(GG, gg, beta)[1]) + sparsity * log(n) + pen
    }
  }

  #ic.err.mat[ic.err.mat < 0] <- Inf
  ic.err.mat[is.nan(ic.err.mat)] <- Inf
  lambda.min <-
    min(lambda.path[apply(ic.err.mat, 1, min) == min(apply(ic.err.mat, 1, min))])
  order.min <-
    min(var.order[apply(ic.err.mat, 2, min) == min(apply(ic.err.mat, 2, min))])


  out <-
    list(
      lambda = lambda.min,
      order.min = order.min,
      error = ic.err.mat,
      lambda.path = lambda.path,
      var.order = var.order
    )
  return(out)
}

#' @title Plotting output for tuning parameter selection in fnets
#' @description Tuning plots for S3 objects of class \code{fnets}.
#' Produces up to two plots visualising CV and IC procedures for selecting tuning parameters and the VAR order.
#' @details See Owens, Cho and Barigozzi (2024+) for further details.
#' @param x \code{fnets} object
#' @param ... additional arguments
#' @return CV/IC plot for the VAR component, and CV plot for the lrpc component (when \code{x$do.lrpc = TRUE}).
#' @seealso \link[fnets]{fnets}
#' @references Owens, D., Cho, H. & Barigozzi, M. (2024+) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. The R Journal (to appear).
#' @importFrom graphics par abline box axis legend matplot
#' @keywords internal
tuning_plot <- function(x, ...){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  data <- attr(x, "data")
  lambda.min <- data$lambda
  order.min <- data$order.min
  error <- data$error
  lambda.path <- data$lambda.path
  var.order <- data$var.order
  args <- attr(x, "args")
  par(mfrow = c(1,1+x$do.lrpc)) ## check

  if(args$tuning == "bic") {
    ylab = "IC"
    main = "IC for VAR parameter estimation"
  } else if(args$tuning == "cv") {
    ylab = "CV"
    main = "CV for VAR parameter estimation"
  }

    par(xpd=FALSE)
    error.plot <- error
    if(any(is.infinite(error.plot))) {
      error.plot[which(is.infinite(error.plot))] <- NA
      keep <- colSums(!is.na(error.plot)) > 0
      error.plot <- error.plot[,keep]
    } else keep <- rep(TRUE, length(var.order))
    matplot(
      lambda.path,
      error.plot,
      type = "b",
      col = (2:(max(var.order) + 1))[keep],
      pch = (2:(max(var.order) + 1))[keep],
      log = "x",
      xlab = "lambda (log scale)",
      ylab = ylab,
      main = main
    )
    abline(v = lambda.min)
    legend(
      "topleft",
      legend = var.order[keep],
      col = (2:(max(var.order) + 1))[keep],
      pch = (2:(max(var.order) + 1))[keep],
      lty = 1
    )

  if (x$do.lrpc) {
    data <- attr(x$lrpc, "data")
    plot(
      data$eta.path,
      data$cv.err,
      type = "b",
      col = 2,
      pch = 2,
      log = "x",
      xlab = "eta (log scale)",
      ylab = "CV error",
      main = "CV for (LR)PC matrix estimation"
    )
    abline(v = data$eta)
  }

}


#' @title Forecasting idiosyncratic VAR process
#' @description Produces forecasts of the idiosyncratic VAR process
#' for a given forecasting horizon by estimating the best linear predictors
#' @param object \code{fnets} object
#' @param x input time series, with each row representing a variable
#' @param cpre output of \link[fnets]{common.predict}
#' @param n.ahead forecast horizon
#' @return a list containing
#' \item{is}{ in-sample estimator of the idiosyncratic component (with each column representing a variable)}
#' \item{fc}{ forecasts of the idiosyncratic component for a given forecasting horizon \code{h} (with each column representing a variable)}
#' \item{n.ahead}{ forecast horizon}
#' @examples
#' \dontrun{
#' out <- fnets(data.unrestricted,
#' do.lrpc = FALSE, var.args = list(n.cores = 2))
#' cpre <- common.predict(out)
#' ipre <- idio.predict(out, cpre)
#' }
#' @importFrom stats as.ts
#' @keywords internal
idio.predict <- function(object, x, cpre, n.ahead = 1) {
  p <- dim(x)[1]
  n <- dim(x)[2]
  # if(attr(object, 'factor') == 'dynamic'){

  xx <- x - object$mean.x
  beta <- object$idio.var$beta
  if(is.null(beta)) beta <- object$beta
  d <- dim(beta)[1] / p
  A <- t(beta)

  is <- xx - t(cpre$is)
  if(n.ahead >= 1) {
    fc <- matrix(0, nrow = p, ncol = n.ahead)
    for (ii in 1:n.ahead) {
      for (ll in 1:d)
        fc[, ii] <-
          fc[, ii] + A[, p * (ll - 1) + 1:p] %*% is[, n + ii - ll]
      is <- cbind(is, fc[, ii])
    }
  } else {
    fc <- NA
  }

  out <- list(is = as.ts(t(is[, 1:n])), fc = as.ts(t(fc)), n.ahead = n.ahead)
  return(out)
}

#' @keywords internal
make.gg <- function(acv, d) {
  p <- dim(acv)[1]
  gg <- matrix(0, nrow = p * d, ncol = p)
  GG <- matrix(0, p * d, p * d)
  for (ll in 1:d) {
    gg[(ll - 1) * p + 1:p, ] <- acv[, , ll + 1]
    for (lll in ll:d) {
      GG[(ll - 1) * p + 1:p, (lll - 1) * p + 1:p] <-
        t(acv[, , 1 + lll - ll])
      GG[(lll - 1) * p + 1:p, (ll - 1) * p + 1:p] <-
        acv[, , 1 + lll - ll]
    }
  }
  out <- list(gg = gg, GG = GG)
  return(out)
}

#' @keywords internal
f.func <- function(GG, gg, A) {
  return(.5 * sum(diag(t(A) %*% GG %*% A - 2 * t(A) %*% gg)))
}

#' @keywords internal
gradf.func <- function(GG, gg, A) {
  return(GG %*% (A) - gg)
}

#' @keywords internal
Q.func <- function(A, A.up, L, GG, gg) {
  Adiff <- (A - A.up)
  return(f.func(GG, gg, A.up) + sum(Adiff * gradf.func(GG, gg, A.up)) + 0.5 * L * norm(Adiff, "F") ^
           2)
}

#' @keywords internal
prox.func <- function(B, lambda, L, GG, gg) {
  b <- B - (1 / L) * gradf.func(GG, gg, B)
  sgn <- sign(b)
  ab <- abs(b)
  sub <- ab - 2 * lambda / L
  sub[sub < 0] <- 0
  out <- sub * sgn
  return(as.matrix(out))
}

#' @title Threshold the entries of the input matrix at a data-driven level
#' @description Threshold the entries of the input matrix at a data-driven level.
#' This can be used to perform edge selection for VAR parameter, inverse innovation covariance, and long-run partial correlation networks.
#' @details See Owens, Cho & Barigozzi (2024+) for more information on the threshold selection process
#' @param mat input parameter matrix
#' @param path.length number of candidate thresholds
#' @return an S3 object of class \code{threshold}, which contains the following fields:
#' \item{threshold}{ data-driven threshold}
#' \item{thr.mat}{ thresholded input matrix}
#' @seealso \link[fnets]{plot.threshold}, \link[fnets]{print.threshold}
#' @example R/examples/thresh_ex.R
#' @importFrom graphics par
#' @references Owens, D., Cho, H. & Barigozzi, M. (2024+) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. The R Journal (to appear).
#' @export
threshold <- function(mat,
                      path.length = 500) {
  path.length <- posint(path.length)

    p <- dim(mat)[1]
    M <- max(abs(mat), 1e-3)
    m <- max(min(abs(mat)), M * .01, 1e-4)
    rseq <-
      round(exp(seq(log(M), log(m), length.out = path.length)), digits = 10)
    cusum <- ratio <- rseq * 0
    for (ii in 1:path.length) {
      A <- mat
      A[abs(A) < rseq[ii]] <- 0
      edges <- sum(A != 0)
      ratio[ii] <- edges / (prod(dim(mat)) - edges)
    }
    dif <- diff(ratio) / diff(rseq)
    for (ii in 2:(path.length - 1))
      cusum[ii] <-
      (mean(dif[2:ii - 1]) - mean(dif[ii:(path.length - 1)])) * (ii / sqrt(path.length)) * (1 - ii / path.length)

    thr <- rseq[which.max(abs(cusum))]

    A <- mat
    A[abs(A) < thr] <- 0


  out <- list(threshold = thr, thr.mat = A)
  attr(out, "seqs") <- list(rseq = rseq, ratio = ratio, dif = dif, cusum = cusum)
  attr(out, "class") <- "threshold"
  return(out)
}

#' @title Plotting the thresholding procedure
#' @method plot threshold
#' @description Plotting method for S3 objects of class \code{threshold}.
#' Produces a plot visualising three diagnostics for the thresholding procedure, with threshold values t_k (x axis) against
#' (i) Ratio_k, the ratio of the number of non-zero to zero entries in the matrix, as the threshold varies
#' (ii) Diff_k, the first difference of \code{Ratio_k}
#' (iii) |CUSUM_k|, the absolute scaled cumulative sums of \code{Diff_k}
#' @details See Owens, Cho and Barigozzi (2024+) for further details.
#' @param x \code{threshold} object
#' @param plots logical vector, which plots to use (Ratio, Diff, CUSUM respectively)
#' @param ... additional arguments
#' @return A network plot produced as per the input arguments
#' @references Owens, D., Cho, H. & Barigozzi, M. (2024+) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. The R Journal (to appear).
#' @seealso \link[fnets]{threshold}
#' @example R/examples/thresh_ex.R
#' @export
plot.threshold <- function(x, plots = c(TRUE, FALSE, TRUE), ...){
  seqs <- attr(x, "seqs")
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if(sum(plots) == 0 | length(plots) != 3){
    warning("setting all plots to TRUE")
    plots <- c(TRUE, TRUE, TRUE)
  }
  par(mfrow = c(1, sum(plots)))
  if(plots[1]){
    plot(seqs$rseq, seqs$ratio, type = "l", xlab = "Threshold", ylab = "Ratio")
    abline(v = x$threshold)
  }
  if(plots[2]){
    plot(seqs$rseq[-1], seqs$dif, type = "l", xlab = "Threshold", ylab = "Diff")
    abline(v = x$threshold)
  }
  if(plots[3]){
    plot(seqs$rseq, abs(seqs$cusum), type = "l", xlab = "Threshold", ylab = "|Cusum|")
    abline(v = x$threshold)
  }

}

#' @title Print threshold
#' @method print threshold
#' @description Prints a summary of a \code{threshold} object
#' @param x \code{threshold} object
#' @param ... not used
#' @return NULL, printed to console
#' @seealso \link[fnets]{threshold}
#' @example R/examples/thresh_ex.R
#' @export
print.threshold <- function(x,
                        ...){
  cat(paste("Thresholded matrix \n"))
  cat(paste("Threshold: ", x$threshold, "\n", sep = ""))
  cat(paste("Non-zero entries: ", sum(x$thr.mat != 0), "/", prod(dim(x$thr.mat)), "\n", sep = ""))
}
