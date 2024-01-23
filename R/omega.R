#' @title Parametric estimation of long-run partial correlations of factor-adjusted VAR processes
#' @description Returns a parametric estimate of long-run partial correlations of the VAR process
#' from the VAR parameter estimates and the inverse of innovation covariance matrix obtained via constrained \code{l1}-minimisation.
#' @details See Barigozzi, Cho and Owens (2024+) for further details, and Cai, Liu and Zhou (2016) for further details on the adaptive estimation procedure.
#' @param object \code{fnets} object
#' @param eta \code{l1}-regularisation parameter; if \code{eta = NULL}, it is selected by cross validation
#' @param tuning.args a list specifying arguments for the cross validation procedure
#' for selecting the tuning parameter involved in long-run partial correlation matrix estimation. It contains:
#' \describe{
#'    \item{\code{n.folds}}{ positive integer number of folds}
#'    \item{\code{path.length}}{ positive integer number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#' }
#' @param lrpc.adaptive whether to use the adaptive estimation procedure
#' @param eta.adaptive \code{l1}-regularisation parameter for Step 1 of the adaptive estimation procedure; if \code{eta.adaptive = NULL}, the default choice is \code{2 * sqrt(log(dim(x)[1])/dim(x)[2])}
#' @param do.correct whether to correct for any negative entries in the diagonals of the inverse of long-run covariance matrix
#' @param do.threshold whether to perform adaptive thresholding of \code{Delta} and \code{Omega} parameter estimators with \link[fnets]{threshold}
#' @param n.cores number of cores to use for parallel computing, see \link[parallel]{makePSOCKcluster}
#' @return a list containing
#' \item{Delta}{ estimated inverse of the innovation covariance matrix}
#' \item{Omega}{ estimated inverse of the long-run covariance matrix}
#' \item{pc}{ estimated innovation partial correlation matrix}
#' \item{lrpc}{ estimated long-run partial correlation matrix}
#' \item{eta}{ \code{l1}-regularisation parameter}
#' \item{lrpc.adaptive}{ input argument }
#' @references Barigozzi, M., Cho, H. & Owens, D. (2024+) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. Journal of Business & Economic Statistics (to appear).
#' @references Cai, T. T., Liu, W., & Zhou, H. H. (2016) Estimating sparse precision matrix: Optimal rates of convergence and adaptive estimation. The Annals of Statistics, 44(2), 455-488.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2024+) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. The R Journal (to appear).
#' @examples
#' \donttest{
#' out <- fnets(data.unrestricted, do.lrpc = FALSE, var.args = list(n.cores = 2))
#' plrpc <- par.lrpc(out, n.cores = 2)
#' out$lrpc <- plrpc
#' out$do.lrpc <- TRUE
#' plot(out, type = "pc", display = "network")
#' plot(out, type = "lrpc", display = "heatmap")
#' }
#' @importFrom parallel detectCores
#' @export
par.lrpc <- function(object,
                     eta = NULL,
                     tuning.args = list(n.folds = 1,
                                        path.length = 10),
                     lrpc.adaptive = FALSE,
                     eta.adaptive = NULL,
                     do.correct = TRUE,
                     do.threshold = FALSE,
                     n.cores = 1) {
  is.var <- FALSE
  p <- dim(object$acv$Gamma_x)[1]
  if(is.null(p)) {
    p <- dim(object$beta)[2]
    is.var <- TRUE
  }
  x <- t(as.matrix(attr(object, "args")$x))
  xx <- x - object$mean.x
  n <- dim(x)[2]

  tuning.args <- check.list.arg(tuning.args)
  n.cores <- posint(n.cores)

  if(!is.var){
    GG <- object$idio.var$Gamma
    A <- t(object$idio.var$beta)
  } else {
    GG <- object$Gamma
    A <- t(object$beta)
  }
  d <- dim(A)[2] / p

  A1 <- diag(1, p)
  for (ll in 1:d)
    A1 <- A1 - A[, (ll - 1) * p + 1:p]

  if (is.null(eta)) {
    dcv <- direct.cv(
      object,
      xx,
      target = "acv",
      symmetric = "min",
      n.folds = tuning.args$n.folds,
      path.length = tuning.args$path.length,
      q = object$q,
      kern.bw = object$kern.bw,
      n.cores = n.cores,
      lrpc.adaptive = lrpc.adaptive,
      eta.adaptive = eta.adaptive,
      is.var = is.var
    )
    eta <- dcv$eta
  } else eta <- max(0, eta)
  if (lrpc.adaptive) {
    if (is.null(eta.adaptive)) {
      eta.adaptive <- 2 * sqrt(log(p) / n)
    }
    Delta <- adaptive.direct.inv.est(
      GG,
      n,
      eta = eta,
      eta.adaptive = eta.adaptive,
      symmetric = "min",
      do.correct = do.correct,
      n.cores = n.cores
    )$DD
  } else {
    Delta <- direct.inv.est(
      GG,
      eta = eta,
      symmetric = "min",
      do.correct = do.correct,
      n.cores = n.cores
    )$DD
  }
  if(do.threshold){
    dd <- diag(Delta)
    Delta <- threshold(Delta)$thr.mat
    diag(Delta) <- dd
  }
  Omega <- 2 * pi * t(A1) %*% Delta %*% A1
  if (do.correct)
    if(!is.null(object$spec$Sigma_i))
      Omega <- correct.diag(Re(object$spec$Sigma_i[, , 1]), Omega)
  if(do.threshold){
    dd <- diag(Omega)
    Omega <- threshold(Omega)$thr.mat
    diag(Omega) <- dd
  }
  pc <- -t(t(Delta) / sqrt(abs(diag(Delta)))) / sqrt(abs(diag(Delta)))
  lrpc <- -t(t(Omega) / sqrt(diag(Omega))) / sqrt(diag(Omega))
  out <-
    list(
      Delta = Delta,
      Omega = Omega,
      pc = pc,
      lrpc = lrpc,
      eta = eta,
      lrpc.adaptive = lrpc.adaptive
    )
  attr(out, "data") <- dcv
  return(out)
}

#' @keywords internal
#' @importFrom parallel detectCores
#' @importFrom graphics abline
direct.cv <-
  function(object,
           xx,
           target = c("spec", "acv"),
           symmetric = c("min", "max", "avg", "none"),
           n.folds = 1,
           path.length = 10,
           q = 0,
           kern.bw = NULL,
           n.cores = 1,
           lrpc.adaptive = FALSE,
           eta.adaptive = NULL,
           is.var = FALSE) {
    n <- ncol(xx)
    p <- nrow(xx)

    if (is.null(kern.bw))
      kern.bw <- 4 * floor((n / log(n)) ^ (1 / 3))
    target <- match.arg(target, c("spec", "acv"))
    if (target == "spec") {
      GG <- Re(object$spec$Sigma_i[, , 1])
      eta.max <- max(abs(GG))
      eta.path <-
        round(exp(seq(
          log(eta.max), log(eta.max * .01), length.out = path.length
        )), digits = 10)
    }
    if (target == "acv") {
      ifelse(is.var, A <- t(object$beta), A <- t(object$idio.var$beta))
      d <- dim(A)[2] / p
      ifelse(is.var, GG <- object$Gamma, GG <- object$idio.var$Gamma)
      eta.max <- max(abs(GG))
      ifelse(lrpc.adaptive, eta.max.2 <- 2 * sqrt(log(p) / n), eta.max.2 <- eta.max)
      eta.path <-
        round(exp(seq(
          log(eta.max.2), log(eta.max * .01), length.out = path.length
        )), digits = 10)
    }

    cv.err <- rep(0, length = path.length)
    ind.list <- split(1:n, ceiling(n.folds * (1:n) / n))
    for (fold in 1:n.folds) {
      train.ind <- 1:ceiling(length(ind.list[[fold]]) * .5)
      train.x <- xx[, ind.list[[fold]][train.ind]]
      test.x <- xx[, ind.list[[fold]][-train.ind]]
      if (target == "spec") {
        train.GG <-
          Re(dyn.pca(train.x, q = q, kern.bw = kern.bw)$spec$Sigma_i[, , 1])
        test.GG <-
          Re(dyn.pca(test.x, q = q, kern.bw = kern.bw)$spec$Sigma_i[, , 1])
      }
      if (target == "acv") {
        train.G0 <-
          dyn.pca(train.x,
                  q = q,
                  kern.bw = kern.bw,
                  mm = d)$acv$Gamma_i
        test.G0 <-
          dyn.pca(test.x,
                  q = q,
                  kern.bw = kern.bw,
                  mm = d)$acv$Gamma_i
        train.GG <- train.G0[, , 1]
        test.GG <- test.G0[, , 1]
        for (ll in 1:d) {
          train.GG <-
            train.GG - A[, (ll - 1) * p + 1:p] %*% train.G0[, , ll + 1]
          test.GG <-
            test.GG - A[, (ll - 1) * p + 1:p] %*% test.G0[, , ll + 1]
        }
      }

      for (ii in 1:path.length) {
        if (lrpc.adaptive) {
          DD <-
            adaptive.direct.inv.est(
              train.GG,
              n = n,
              eta = eta.path[ii],
              eta.adaptive = eta.adaptive,
              symmetric = symmetric,
              n.cores = n.cores
            )$DD
        } else {
          DD <-
            direct.inv.est(
              train.GG,
              eta = eta.path[ii],
              symmetric = symmetric,
              n.cores = n.cores
            )$DD
        }
        DG <- DD %*% test.GG
        sv <- svd(DG, nu = 0, nv = 0)
        cv.err[ii] <- cv.err[ii] + sum(sv$d) - sum(log(sv$d)) - p
      }
    }

    eta.min <- eta.path[which.min(cv.err)]

    out <- list(eta = eta.min,
                cv.error = cv.err,
                eta.path = eta.path)
    return(out)
  }

#' @keywords internal
#' @importFrom parallel makePSOCKcluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lpSolve lp
direct.inv.est <-
  function(GG,
           eta = NULL,
           symmetric = c("min", "max", "avg", "none"),
           do.correct = FALSE,
           n.cores = 1) {
    p <- dim(GG)[1]
    f.obj <- rep(1, 2 * p)
    f.con <- rbind(-GG, GG)
    f.con <- cbind(f.con, -f.con)
    f.dir <- rep("<=", 2 * p)

    cl <- parallel::makePSOCKcluster(n.cores)
    doParallel::registerDoParallel(cl)

    ii <- 1
    DD <-
      foreach::foreach(
        ii = 1:p,
        .combine = "cbind",
        .multicombine = TRUE,
        .export = c("lp")
      ) %dopar% {
        ee <- rep(0, p)
        ee[ii] <- 1
        b1 <- rep(eta, p) - ee
        b2 <- rep(eta, p) + ee
        f.rhs <- c(b1, b2)
        lpout <- lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs)
        lpout$solution[1:p] - lpout$solution[-(1:p)]
      }
    parallel::stopCluster(cl)

    DD <- make.symmetric(DD, symmetric)
    if (do.correct)
      DD <- correct.diag(GG, DD)

    out <- list(DD = DD,
                eta = eta,
                symmetric = symmetric)
    return(out)
  }

#' @keywords internal
#' @importFrom parallel makePSOCKcluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom lpSolve lp
adaptive.direct.inv.est <-
  function(GG,
           n,
           eta = NULL,
           eta.adaptive = NULL,
           symmetric = c("min", "max", "avg", "none"),
           do.correct = FALSE,
           n.cores = 1) {
    p <- dim(GG)[1]
    f.obj <- rep(1, 2 * p)
    GG.n <- GG + diag(1 / n, p) # add ridge
    dGG <- pmax(diag(GG), 0)
    f.dir <- rep("<=", 2 * p)

    f.con.0 <- rbind(-GG.n, GG.n) # initialise
    f.con.0 <- cbind(f.con.0, -f.con.0)
    ## Step 1 //
    cl <- parallel::makePSOCKcluster(n.cores)
    doParallel::registerDoParallel(cl)
    if (is.null(eta.adaptive))
      eta.adaptive <- 2 * sqrt(log(p) / n)
    ii <- 1
    step1.index <- which(dGG <= sqrt(n / log(p)))
    f.con.1 <- rbind(f.con.0, 0)
    f.dir.1 <- c(f.dir, "==") # diagonals are positive
    DD.1 <-
      foreach::foreach(
        ii = step1.index,
        .combine = "cbind",
        .multicombine = TRUE,
        .export = c("lp")
      ) %dopar% {
        f.con.ii <- f.con.1
        ii.replace <- eta.adaptive * pmax(dGG, dGG[ii])
        f.con.ii[1:(2 * p), ii] <-
          f.con.ii[1:(2 * p), ii] - ii.replace # mutate cols
        f.con.ii[1:(2 * p), ii + p] <-
          f.con.ii[1:(2 * p), ii + p] - ii.replace # mutate cols
        f.con.ii[2 * p + 1, ii + p] <- 1 # diagonals are positive
        ee <- rep(0, p)
        ee[ii] <- 1
        b1 <- -ee
        b2 <- ee
        f.rhs <- c(b1, b2, 0)
        lpout <- lpSolve::lp("min", f.obj, f.con.ii, f.dir.1, f.rhs)
        lpout$solution[1:p] - lpout$solution[p + (1:p)]
      }
    dDD.1 <- diag(DD.1)
    dDD.1[!step1.index] <- sqrt(log(p) / n)
    ## Step 2 //
    if (is.null(eta)) {
      eta <- 2 * sqrt(log(p) / n)
    }
    ii <- 1
    DD.2 <-
      foreach::foreach(
        ii = 1:p,
        .combine = "cbind",
        .multicombine = TRUE,
        .export = c("lp")
      ) %dopar% {
        ee <- rep(0, p)
        ee[ii] <- 1
        bb <- eta * sqrt(dGG) * sqrt(dDD.1[ii])
        b1 <- bb - ee
        b2 <- bb + ee
        f.rhs <- c(b1, b2)
        lpout <- lpSolve::lp("min", f.obj, f.con.0, f.dir, f.rhs)
        lpout$solution[1:p] - lpout$solution[-(1:p)]
      }
    parallel::stopCluster(cl)
    DD.2 <- make.symmetric(DD.2, symmetric)
    if (do.correct) {
      tmp <- gen.inverse(GG)
      ind <- which(diag(DD.2) == 0)
      diag(DD.2)[ind] <- tmp[ind]
    }
    out <- list(DD = DD.2,
                eta = eta,
                symmetric = symmetric)
    return(out)
  }


#' @keywords internal
gen.inverse <- function(GG) {
  p <- dim(GG)[1]
  sv <- svd(GG)
  L <- GG * 0
  diag(L)[sv$d > 0] <- sv$d[sv$d > 0]
  return(diag(sv$u %*% L %*% t(sv$u)))
}

#' @keywords internal
correct.diag <- function(GG, DD) {
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
