#' @title \code{l1}-regularised Yule-Walker estimation for VAR processes
#' @description Estimates the VAR parameter matrices via \code{l1}-regularised Yule-Walker estimation
#' and innovation covariance matrix via constrained \code{l1}-minimisation.
#' @details Further information can be found in Barigozzi, Cho and Owens (2022).
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param method a string specifying the method to be adopted for VAR process estimation; possible values are:
#' \itemize{
#'    \item{\code{"lasso"}}{ Lasso-type \code{l1}-regularised \code{M}-estimation}
#'    \item{\code{"ds"}}{ Dantzig Selector-type constrained \code{l1}-minimisation}
#' }
#' @param var.order order of the VAR process; if a vector of integers is supplied, the order is chosen via \code{tuning}
#' @param lambda regularisation parameter; if \code{lambda = NULL}, \code{tuning} is employed to select the parameter
#' @param tuning.args a list specifying arguments for \code{tuning}
#' for selecting the regularisation parameter (and VAR order). It contains:
#' \itemize{
#'    \item{\code{tuning}} a string specifying the selection procedure for \code{var.order} and \code{lambda}; possible values are:
#'    \itemize{
#'       \item{\code{"cv"}}{ cross validation}
#'       \item{\code{"bic"}}{ information criterion}
#'    }
#'    \item{\code{n.folds}}{ if \code{tuning = "cv"}, positive integer number of folds}
#'    \item{\code{penalty}}{ if \code{tuning = "bic"}, penalty multiplier between 0 and 1; if \code{penalty = NULL}, it is set to \code{1/(1+exp(dim(x)[1])/dim(x)[2]))}} by default
#'    \item{\code{path.length}}{ positive integer number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{\code{do.plot}}{ whether to plot the output of the cross validation step}
#' }
#' @param do.threshold whether to perform adaptive thresholding of VAR parameter estimator with \link[fnets]{threshold}
#' @param n.iter maximum number of descent steps, by default depends on \code{var.order}; applicable when \code{method = "lasso"}
#' @param tol numerical tolerance for increases in the loss function; applicable when \code{method = "lasso"}
#' @param n.cores number of cores to use for parallel computing, see \link[parallel]{makePSOCKcluster}; applicable when \code{method = "ds"}
#' @return a list which contains the following fields:
#' \item{beta}{ estimate of VAR parameter matrix; each column contains parameter estimates for the regression model for a given variable}
#' \item{Gamma}{ estimate of the innovation covariance matrix}
#' \item{lambda}{ regularisation parameter}
#' \item{var.order}{ VAR order}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#' @example R/examples/var_ex.R
#' @importFrom parallel detectCores
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. arXiv preprint arXiv:2301.11675.
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
                        path.length = 10,
                        do.plot = FALSE
                      ),
                      do.threshold = FALSE,
                      n.iter = NULL,
                      tol = 0,
                      n.cores = min(parallel::detectCores() - 1, 3)) {
  p <- dim(x)[1]
  n <- dim(x)[2]

  tuning.args <- check.list.arg(tuning.args)

  method <- match.arg(method, c("lasso", "ds"))
  tuning <- match.arg(tuning.args$tuning, c("cv", "bic"))
  if(center)
    mean.x <- apply(x, 1, mean)
  else
    mean.x <- rep(0, p)
  xx <- x - mean.x
  dpca <- dyn.pca(xx, q = 0)
  acv <- dpca$acv


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
      do.plot = tuning.args$do.plot,
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
      kern.bw = NULL,
      do.plot = tuning.args$do.plot
    )
  }

  mg <- make.gg(acv$Gamma_i, icv$var.order)
  gg <- mg$gg
  GG <- mg$GG

  if(method == "lasso"){
    if(is.null(n.iter)) n.iter <- var.order*100
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
  ive$var.order <- icv$var.order
  ive$mean.x <- mean.x
  if(do.threshold)
    ive$beta <- threshold(ive$beta, do.plot = tuning.args$do.plot)

  attr(ive, "class") <- "fnets"
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
          if(f.func(GG, gg, prox) <= Q.func(prox, y, L.bar, GG, gg)) {
            found <- TRUE
          } else {
            L.bar <- L.bar * gamma
          }
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
           n.cores = min(parallel::detectCores() - 1, 3)) {
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
#' @importFrom graphics abline legend matplot
#' @keywords internal
yw.cv <- function(xx,
                  method = c("lasso", "ds"),
                  lambda.max = NULL,
                  var.order = 1,
                  n.folds = 1,
                  path.length = 10,
                  q = 0,
                  kern.bw = NULL,
                  do.plot = FALSE,
                  n.cores = min(parallel::detectCores() - 1, 3)) {
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

  if(do.plot) {
    par(xpd=FALSE)
    cv.err.mat.plot <- cv.err.mat
    if(any(is.infinite(cv.err.mat.plot))) {
      cv.err.mat.plot[which(is.infinite(cv.err.mat.plot))] <- NA
      keep <- colSums(!is.na(cv.err.mat.plot)) > 0
      cv.err.mat.plot <- cv.err.mat.plot[,keep]
    } else keep <- rep(1, length(var.order))
    keep <- as.logical(keep)
    matplot(
      lambda.path,
      cv.err.mat.plot,
      type = "b",
      col = (2:(max(var.order) + 1))[keep],
      pch = (2:(max(var.order) + 1))[keep],
      log = "x",
      xlab = "lambda (log scale)",
      ylab = "CV error",
      main = "CV for VAR parameter estimation"
    )
    abline(v = lambda.min)
    legend(
      "topleft",
      legend = var.order[keep],
      col = (2:(max(var.order) + 1))[keep],
      pch = (2:(max(var.order) + 1))[keep],
      lty = 1
    )
  }

  out <-
    list(
      lambda = lambda.min,
      var.order = order.min,
      cv.error = cv.err.mat,
      lambda.path = lambda.path
    )
  return(out)
}


# ebic

#' @title logarithmic factorial of `n`
#' @keywords internal
log.factorial <- function(n)
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
    2 * penalty * (log.factorial(p^2 * d) - log.factorial(sparsity) - log.factorial(p ^
                                                                                        2 * d - sparsity))
}


#' @title Information criterion for factor-adjusted VAR estimation
#' @importFrom graphics abline legend matplot
#' @keywords internal
yw.ic <- function(xx,
                  method = c("lasso", "ds"),
                  lambda.max = NULL,
                  var.order = 1,
                  penalty = NULL,
                  path.length = 10,
                  q = 0,
                  kern.bw = NULL,
                  do.plot = FALSE) {
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
      beta <- threshold(beta, do.plot = FALSE)$thr.mat
      sparsity <- sum(beta[, 1] != 0)
      if(sparsity == 0) {
        pen <- 0
      } else {
        pen <-
          2 * penalty * (log.factorial(var.order[jj] * p) - log(sparsity) - log(var.order[jj] * p - sparsity))
      }
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

  if(do.plot) {
    par(xpd=FALSE)
    ic.err.mat.plot <- ic.err.mat
    if(any(is.infinite(ic.err.mat.plot))) {
      ic.err.mat.plot[which(is.infinite(ic.err.mat.plot))] <- NA
      keep <- colSums(!is.na(ic.err.mat.plot)) > 0
      ic.err.mat.plot <- ic.err.mat.plot[,keep]
    } else keep <- rep(1, length(var.order))
    matplot(
      lambda.path,
      ic.err.mat.plot,
      type = "b",
      col = (2:(max(var.order) + 1))[keep],
      pch = (2:(max(var.order) + 1))[keep],
      log = "x",
      xlab = "lambda (log scale)",
      ylab = "IC",
      main = "IC for VAR parameter estimation"
    )
    abline(v = lambda.min)
    legend(
      "topleft",
      legend = var.order[keep],
      col = (2:(max(var.order) + 1))[keep],
      pch = (2:(max(var.order) + 1))[keep],
      lty = 1
    )
  }

  out <-
    list(
      lambda = lambda.min,
      var.order = order.min,
      ic.error = ic.err.mat,
      lambda.path = lambda.path
    )
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
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. arXiv preprint arXiv:2301.11675.
#' @examples
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.unrestricted(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#' out <- fnets(x, q = NULL, var.order = 1, var.method = "lasso",
#' do.lrpc = FALSE, var.args = list(n.cores = 2))
#' cpre <- common.predict(out, x, h = 1, r = NULL)
#' ipre <- idio.predict(out, x, cpre, h = 1)
#' @export
idio.predict <- function(object, x, cpre, h = 1) {
  p <- dim(x)[1]
  n <- dim(x)[2]
  # if(attr(object, 'factor') == 'dynamic'){

  xx <- x - object$mean.x
  beta <- object$idio.var$beta
  d <- dim(beta)[1] / p
  A <- t(beta)

  is <- xx - cpre$is
  if(h >= 1) {
    fc <- matrix(0, nrow = p, ncol = h)
    for (ii in 1:h) {
      for (ll in 1:d)
        fc[, ii] <-
          fc[, ii] + A[, p * (ll - 1) + 1:p] %*% is[, n + ii - ll]
      is <- cbind(is, fc[, ii])
    }
  } else {
    fc <- NA
  }
  # }
  # if(attr(object, 'factor') == 'static'){
  #   is <- matrix(0, p, n)
  #   fc <- 0
  # }


  out <- list(is = is[, 1:n], fc = fc, h = h)
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

#' @title Edge selection for VAR parameter, inverse innovation covariance, and long-run partial correlation matrices
#' @description Threshold the entries of the input matrix at a data-driven level to perform edge selection
#' @details See Liu, Zhang, and Liu (2021) for more information on the threshold selection process
#' @param mat input parameter matrix
#' @param path.length number of candidate thresholds
#' @param do.plot whether to plot thresholding output
#' @return a list which contains the following fields:
#' \item{threshold}{ data-driven threshold}
#' \item{thr.mat}{ thresholded input matrix}
#' @examples
#' \donttest{
#' set.seed(123)
#' A <- diag(.7, 50) + rnorm(50^2, 0, .1)
#' threshold.A <- threshold(A)
#' }
#' @importFrom graphics par
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network analysis for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Liu, B., Zhang, X. & Liu, Y. (2021) Simultaneous Change Point Inference and Structure Recovery for High Dimensional Gaussian Graphical Models. Journal of Machine Learning Research, 22(274), 1--62.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. arXiv preprint arXiv:2301.11675.
#' @export
threshold <- function(mat,
                      path.length = 500,
                      do.plot = FALSE) {
  if(!is.null(attr(mat, "thresholded"))) {
    warning("This matrix has already been thresholded. Returning input.")
    A <- mat
    thr <- 0
  } else {
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

    if(do.plot) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      par(mfrow = c(1, 2))
      plot(rseq, ratio, type = "l", xlab = "threshold")
      abline(v = thr)
      # plot(rseq[-1], dif, type = "l", xlab = "threshold")
      # abline(v = thr)
      plot(rseq, abs(cusum), type = "l", xlab = "threshold")
      abline(v = thr)
    }

    A <- mat
    A[abs(A) < thr] <- 0

    attr(A, "thresholded") <- TRUE
  }

  out <- list(threshold = thr, thr.mat = A)
  return(out)
}
