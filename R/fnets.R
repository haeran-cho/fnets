#' @title Factor-adjusted network estimation
#' @description Operating under factor-adjusted vector autoregressive (VAR) model,
#' the function estimates the spectral density and autocovariance matrices of the factor-driven common component and the idiosyncratic VAR process,
#' the impulse response functions and common shocks for the common component,
#' and VAR parameters, innovation covariance matrix and long-run partial correlations for the idiosyncratic component.
#' @details See Barigozzi, Cho and Owens (2021) for further details.
#' List arguments do not need to be specified with all list components; any missing entries will be filled in with the default argument.
#'
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param factor.model a string specifying the method to be adopted for factor model estimation; possible values are:
#' \itemize{
#'    \item{\code{"unrestricted"}}{ unrestricted factor model}
#'    \item{\code{"restricted"}}{ restricted factor model}
#' }
#' @param q number of unrestricted factors. If \code{q = NULL}, the factor number is estimated by an information criterion-based approach of Hallin and Liška (2007) or Bai and Ng (2002), or eigenvalue ratio, see \link[fnets]{hl.factor.number} and \link[fnets]{bn.factor.number} for further details
#' @param q.method a string specifying the factor number selection method when \code{factor.model = "restricted"}; possible values are:
#' \itemize{
#'    \item{\code{"bn"}}{ information criteria of Bai and Ng (2002)}
#'    \item{\code{"er"}}{ eigenvalue ratio}
#' }
#' @param ic.op choice of the information criterion, see \link[fnets]{hl.factor.number} and \link[fnets]{bn.factor.number} for further details
#' @param kern.const constant multiplied to \code{floor((dim(x)[2]/log(dim(x)[2]))^(1/3)))} which determines the kernel bandwidth for dynamic PCA
#' @param common.args a list specifying the tuning parameters required for estimating the impulse response functions and common shocks. It contains:
#' \itemize{
#'    \item{\code{var.order}}{ order of the blockwise VAR representation of the common component. If \code{var.order = NULL}, it is selected blockwise by Schwarz criterion}
#'    \item{\code{max.var.order}}{ maximum blockwise VAR order for the Schwarz criterion}
#'    \item{\code{trunc.lags}}{ truncation lag for impulse response function estimation}
#'    \item{\code{n.perm}}{ number of cross-sectional permutations involved in impulse response function estimation}
#' }
#' @param idio.var.order order of the idiosyncratic VAR process; if a vector of integers is supplied, the order is chosen via \code{tuning}
#' @param idio.method a string specifying the method to be adopted for idiosyncratic VAR process estimation; possible values are:
#' \itemize{
#'        \item{\code{"lasso"}}{ Lasso-type \code{l1}-regularised \code{M}-estimation}
#'        \item{\code{"ds"}}{ Dantzig Selector-type constrained \code{l1}-minimisation}
#' }
#' @param idio.args a list specifying the tuning parameters required for estimating the idiosyncratic VAR process. It contains:
#' \itemize{
#'    \item{\code{tuning}} a string specifying the selection procedure for \code{idio.var.order} and \code{lambda}; possible values are:
#'    \itemize{
#'       \item{\code{"cv"}}{ cross validation}
#'       \item{\code{"ic"}}{ information criterion}
#'    }
#'    \item{\code{n.iter}}{ maximum number of descent steps; applicable when \code{idio.method = "lasso"}}
#'    \item{\code{tol}}{ numerical tolerance for increases in the loss function; applicable when \code{idio.method = "lasso"}}
#'    \item{\code{n.cores}}{ number of cores to use for parallel computing, see \link[parallel]{makePSOCKcluster}; applicable when \code{idio.method = "ds"}}
#' }
#' @param idio.threshold whether to perform adaptive thresholding of \code{beta} with \link[fnets]{threshold}
#' @param lrpc.method a string specifying the type of estimator for long-run partial correlation matrix estimation; possible values are:
#' \itemize{
#'    \item{\code{"par"}}{ parametric estimator based on the VAR model assumption}
#'    \item{\code{"npar"}}{ nonparametric estimator from inverting the long-run covariance matrix of the idiosyncratic component via constrained \code{l1}-minimisation}
#'    \item{\code{"none"}}{ do not estimate the long-run partial correlation matrix}
#' }
#' @param lrpc.adaptive whether to use the adaptive estimation procedure
#' @param cv.args a list specifying arguments for \code{tuning}
#' for selecting the tuning parameters involved in VAR parameter and (long-run) partial correlation matrix estimation. It contains:
#' \itemize{
#'    \item{\code{n.folds}}{ if \code{tuning = "cv"}, number of folds}
#'    \item{\code{penalty}}{ if \code{tuning = "ic"}, penalty multiplier between 0 and 1; if \code{penalty = NULL}, defaults to \code{1/(1+exp(dim(x)[1])/dim(x)[2]))}}
#'    \item{\code{path.length}}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{\code{do.plot}}{ whether to plot the output of the cross validation step}
#' }
#' @return an S3 object of class \code{fnets}, which contains the following fields:
#' \item{q}{ number of factors}
#' \item{spec}{ if \code{factor.model = "unrestricted"} a list containing estimates of the spectral density matrices for \code{x}, common and idiosyncratic components}
#' \item{acv}{ a list containing estimates of the autocovariance matrices for \code{x}, common and idiosyncratic components}
#' \item{common.irf}{  if \code{factor.model = "unrestricted"} and \code{q >= 1}, a list containing estimators of the impulse response functions (as an array of dimension \code{(p, q, trunc.lags + 2)})
#' and common shocks (an array of dimension \code{(q, n)}) for the common component}
#' \item{lam}{ if \code{factor.model = "restricted"}, factor loadings}
#' \item{f}{ if \code{factor.model = "restricted"}, factor series}
#' \item{idio.var}{ a list containing the following fields:
#' \itemize{
#' \item{\code{beta}}{ estimate of VAR parameter matrix; each column contains parameter estimates for the regression model for a given variable}
#' \item{\code{Gamma}}{ estimate of the innovation covariance matrix}
#' \item{\code{lambda}}{ regularisation parameter}
#' \item{\code{convergence}}{ returned when \code{idio.method = "lasso"}; indicates whether a convergence criterion is met}
#' \item{\code{var.order}}{ VAR order}
#' }}
#' \item{lrpc}{ see the output of \link[fnets]{par.lrpc} if \code{lrpc.method = 'par'}
#' and that of \link[fnets]{npar.lrpc} if \code{lrpc.method = 'npar'}}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#' \item{idio.method}{ input parameter}
#' \item{lrpc.method}{ input parameter}
#' \item{kern.const}{ input parameter}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @references Bai, J. & Ng, S. (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221. \cr
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.unrestricted(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#' out <- fnets(x,
#'   q = NULL, idio.var.order = 1, idio.method = "lasso",
#'   lrpc.method = "par", cv.args = list(n.folds = 1, path.length = 10, do.plot = TRUE)
#' )
#' pre <- predict(out, x, h = 1, common.method = "unrestricted")
#' plot(out, type = "granger", display = "network", threshold = .05)
#' plot(out, type = "lrpc", display = "heatmap", threshold = .05)
#' }
#' @seealso \link[fnets]{predict.fnets}, \link[fnets]{plot.fnets}
#' @importFrom graphics par
#' @export
fnets <- function(x, center = TRUE, factor.model = c("unrestricted", "restricted"), q = NULL, q.method = c("bn","er"), ic.op = NULL, kern.const = 4,
                  common.args = list(var.order = NULL, max.var.order = NULL, trunc.lags = 20, n.perm = 10),
                  idio.var.order = 1, idio.method = c("lasso", "ds"),
                  idio.args = list(tuning = c("cv", "ic"), n.iter = 100, tol = 0, n.cores = min(parallel::detectCores() - 1, 3)),
                  idio.threshold = FALSE,
                  lrpc.method = c("par", "npar", "none"), lrpc.adaptive = FALSE,
                  cv.args = list(n.folds = 1, penalty = NULL, path.length = 10, do.plot = FALSE)) {
  p <- dim(x)[1]
  n <- dim(x)[2]

  idio.args <- check.list.arg(idio.args)
  cv.args <- check.list.arg(cv.args)

  q.method <- match.arg(q.method, c("bn", "er"))
  idio.method <- match.arg(idio.method, c("lasso", "ds"))
  tuning <- match.arg(idio.args$tuning, c("cv", "ic"))
  factor.model <- match.arg(factor.model, c("unrestricted", "restricted"))
  lrpc.method <- match.arg(lrpc.method, c("par", "npar", "none"))
  if (center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  if (factor.model == "restricted") {
    spca <- static.pca(xx, q.max = NULL, q = q, q.method = q.method, ic.op = ic.op)
    q <- spca$q
    lam <- spca$lam
    f <- spca$f
    spec <- NULL
    acv <- spca$acv
    cve <- NULL
  }
  if (factor.model == "unrestricted") {
    ## dynamic pca
    dpca <- dyn.pca(xx, q, ic.op, kern.const)
    q <- dpca$q
    spec <- dpca$spec
    acv <- dpca$acv
    lam <- NULL
    f <- NULL
    ## common VAR estimation
    cve <- common.irf.estimation(xx,
      Gamma_c = acv$Gamma_c, q = q,
      var.order = common.args$var.order, max.var.order = common.args$max.var.order,
      trunc.lags = common.args$trunc.lags, n.perm = common.args$n.perm
    )
  }

  ## idio estimation
  if (cv.args$do.plot) par(mfrow = c(1, 1 + lrpc.method %in% c("par", "npar")))
  if (tuning == "cv") {
    icv <- yw.cv(xx,
      method = idio.method,
      lambda.max = NULL, var.order = idio.var.order,
      n.folds = cv.args$n.folds, path.length = cv.args$path.length,
      q = q, kern.const = kern.const, do.plot = cv.args$do.plot
    )
  }

  if (tuning == "ic") {
    icv <- yw.ic(xx,
      method = idio.method,
      lambda.max = NULL, var.order = idio.var.order,
      penalty = cv.args$penalty, path.length = cv.args$path.length,
      q = q, kern.const = kern.const, do.plot = cv.args$do.plot
    )
  }


  mg <- make.gg(acv$Gamma_i, icv$var.order)
  gg <- mg$gg
  GG <- mg$GG
  if (idio.method == "lasso") ive <- var.lasso(GG, gg, lambda = icv$lambda, symmetric = "min", n.iter = idio.args$n.iter, tol = idio.args$tol)
  if (idio.method == "ds") ive <- var.dantzig(GG, gg, lambda = icv$lambda, symmetric = "min", n.cores = idio.args$n.cores)
  ive$var.order <- icv$var.order
  if (idio.threshold) ive$beta <- threshold(ive$beta, do.plot = cv.args$do.plot)$thr.mat

  out <- list(
    q = q, spec = spec, lam = lam, f = f, acv = acv,
    common.irf = cve, idio.var = ive, mean.x = mean.x,
    idio.method = idio.method, lrpc.method = lrpc.method, kern.const = kern.const
  )

  attr(out, "factor") <- factor.model

  ## lrpc estimation
  if (lrpc.method %in% c("par", "npar")) {
    if (lrpc.method == "par") lrpc <- par.lrpc(out, x, eta = NULL, cv.args = cv.args, lrpc.adaptive = lrpc.adaptive)
    if (lrpc.method == "npar") lrpc <- npar.lrpc(out, x, eta = NULL, cv.args = cv.args)
    out$lrpc <- lrpc
  } else {
    out$lrpc <- NA
  }

  attr(out, "class") <- "fnets"
  return(out)
}

#' @title Dynamic PCA
#' @description Performs principal components analysis in frequency domain for identifying common and idiosyncratic components.
#' @param xx centred input time series matrix, with each row representing a variable
#' @param q number of factors. If \code{q = NULL}, the factor number is estimated by an information criterion-based approach of Hallin and Liška (2007)
#' @param ic.op choice of the information criterion. Currently the three options from Hallin and Liška (2007) (\code{ic.op = 1, 2} or \code{3}) and
#' their variations with logarithm taken on the cost (\code{ic.op = 4, 5} or \code{6}) are implemented,
#' with \code{ic.op = 5} recommended as a default choice based on numerical experiments
#' @param kern.const constant multiplied to \code{floor((dim(x)[2]/log(dim(x)[2]))^(1/3)))} which determines the kernel bandwidth for dynamic PCA
#' @param mm bandwidth; if \code{mm = NULL}, it is chosen using \code{kern.const}
#' @return a list containing
#' \item{q}{ number of factors}
#' \item{hl}{ if \code{q = NULL}, the output from \link[fnets]{hl.factor.number}}
#' \item{spec}{ a list containing the estimates of the spectral density matrices for \code{x}, common and idiosyncratic components}
#' \item{acv}{ a list containing estimates of the autocovariance matrices for \code{x}, common and idiosyncratic components}
#' \item{kern.const}{ input parameter}
#' @importFrom stats fft
#' @keywords internal
dyn.pca <- function(xx, q = NULL, ic.op = 5, kern.const = 4, mm = NULL) {
  p <- dim(xx)[1]
  n <- dim(xx)[2]

  if (is.null(ic.op)) ic.op <- 5
  if (is.null(mm)) mm <- min(max(1, kern.const * floor((n / log(n))^(1 / 3))), floor(n / 4) - 1) else mm <- min(max(mm, 1, kern.const * floor((n / log(n))^(1 / 3))), floor(n / 4) - 1)
  len <- 2 * mm
  w <- Bartlett.weights(((-mm):mm) / mm)

  ## dynamic pca

  if (!is.null(q)) {
    q <- as.integer(q)
    hl <- NA
    Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * mm + 1))
    for (h in 0:(mm - 1)) {
      Gamma_x[, , h + 1] <- xx[, 1:(n - h)] %*% t(xx[, 1:(n - h) + h]) / n
      Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + mm + 1]
      if (h != 0) {
        Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_xw[, , 2 * mm + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
      }
    }
    Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)
    sv <- list(1:(mm + 1))
    if (q > 0) for (ii in 1:(mm + 1)) sv[[ii]] <- svd(Sigma_x[, , ii], nu = q, nv = 0)
  }
  if (is.null(q)) {
    q.max <- min(50, floor(sqrt(min(n - 1, p))))
    hl <- hl.factor.number(xx, q.max, mm, w, center = FALSE)
    q <- hl$q.hat[ic.op]
    Gamma_x <- hl$Gamma_x
    Sigma_x <- hl$Sigma_x
    sv <- hl$sv
  }

  Gamma_c <- Gamma_i <- Sigma_c <- Sigma_i <- Sigma_x * 0
  if (q >= 1) {
    for (ii in 1:(mm + 1)) {
      Sigma_c[, , ii] <- sv[[ii]]$u[, 1:q, drop = FALSE] %*% diag(sv[[ii]]$d[1:q], q) %*% Conj(t(sv[[ii]]$u[, 1:q, drop = FALSE]))
      if (ii > 1) {
        Sigma_c[, , 2 * mm + 1 - (ii - 1) + 1] <- Conj(Sigma_c[, , ii])
      }
    }
    Gamma_c <- aperm(apply(Sigma_c, c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) * (2 * pi) / (2 * mm + 1)
    Gamma_c <- Re(Gamma_c)
  }
  Sigma_i <- Sigma_x - Sigma_c
  Gamma_i <- Gamma_x - Gamma_c

  spec <- list(Sigma_x = Sigma_x, Sigma_c = Sigma_c, Sigma_i = Sigma_i)
  acv <- list(Gamma_x = Gamma_x, Gamma_c = Re(Gamma_c), Gamma_i = Re(Gamma_i))

  out <- list(q = q, hl = hl, spec = spec, acv = acv, kern.const = kern.const)
  return(out)
}

#' @title Factor number estimator of Hallin and Liška (2007)
#' @description Estimates the number of factors by minimising an information criterion over sub-samples of the data.
#' Currently the three information criteria proposed in Hallin and Liška (2007) (\code{ic.op = 1, 2} or \code{3})
#' and their variations with logarithm taken on the cost (\code{ic.op = 4, 5} or \code{6}) are implemented,
#' with \code{ic.op = 5} recommended as a default choice based on numerical experiments.
#' @details See Hallin and Liška (2007) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param q.max maximum number of factors; if \code{q.max = NULL}, a default value is selected as \code{min(50, floor(sqrt(min(dim(x)[2] - 1, dim(x)[1]))))}
#' @param mm integer representing the kernel bandwidth
#' @param w vector of length \code{2 * mm + 1} containing symmetric weights; if \code{w = NULL}, default weights are generated using the Bartlett kernel and \code{mm}
#' @param do.plot whether to plot the values of six information criteria
#' @param center whether to de-mean the input \code{x} row-wise
#' @return a list containing
#' \item{q.hat}{ a vector containing minimisers of the six information criteria}
#' \item{Gamma_x}{ an array containing the estimates of the autocovariance matrices of \code{x} at \code{2 * mm + 1} lags}
#' \item{Sigma_x}{ an array containing the estimates of the spectral density matrices of \code{x} at \code{2 * mm + 1} Fourier frequencies}
#' \item{sv}{ a list containing the singular value decomposition of \code{Sigma_x}}
#' @example R/examples/hl_ex.R
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @importFrom graphics par abline box axis legend
#' @importFrom stats var
#' @export
hl.factor.number <- function(x, q.max = NULL, mm, w = NULL, do.plot = FALSE, center = TRUE) {
  p <- dim(x)[1]
  n <- dim(x)[2]
  if (is.null(q.max)) q.max <- min(50, floor(sqrt(min(n - 1, p))))

  if (center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  if (is.null(w)) w <- Bartlett.weights(((-mm):mm) / mm)

  p.seq <- floor(3 * p / 4 + (1:10) * p / 40)
  n.seq <- n - (9:0) * floor(n / 20)
  const.seq <- seq(.001, 2, by = .01)
  IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 2 * 3))

  Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * mm + 1))

  for (kk in 1:10) {
    nn <- n.seq[kk]
    pp <- p.seq[kk]
    pen <- c(
      (1 / mm^2 + sqrt(mm / nn) + 1 / pp) * log(min(pp, mm^2, sqrt(nn / mm))),
      1 / sqrt(min(pp, mm^2, sqrt(nn / mm))),
      1 / min(pp, mm^2, sqrt(nn / mm)) * log(min(pp, mm^2, sqrt(nn / mm)))
    )

    for (h in 0:(mm - 1)) {
      Gamma_x[, , h + 1] <- xx[, 1:(nn - h)] %*% t(xx[, 1:(nn - h) + h]) / nn
      Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + mm + 1]
      if (h != 0) {
        Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_xw[, , 2 * mm + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
      }
    }
    Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)
    sv <- list(1:(mm + 1))

    tmp <- rep(0, q.max + 1)
    for (ii in 1:(mm + 1)) {
      if (kk == length(n.seq)) nu <- q.max else nu <- 0
      sv[[ii]] <- svd(Sigma_x[1:pp, 1:pp, ii], nu = nu, nv = 0)
      dd <- sum(sv[[ii]]$d)
      tmp[1] <- tmp[1] + dd / pp / (2 * mm + 1)
      for (jj in 1:q.max) {
        dd <- dd - sv[[ii]]$d[jj]
        tmp[jj + 1] <- tmp[jj + 1] + dd / pp / (2 * mm + 1)
      }
      for (jj in 1:length(const.seq)) {
        for (pen.op in 1:3) {
          IC[, jj, kk, 3 * 0 + pen.op] <- tmp + (0:q.max) * const.seq[jj] * pen[pen.op]
          IC[, jj, kk, 3 * 1 + pen.op] <- log(tmp) + (0:q.max) * const.seq[jj] * pen[pen.op]
        }
      }
    }
  }

  q.mat <- apply(IC, c(2, 3, 4), which.min)
  Sc <- apply(q.mat, c(1, 3), var)
  q.hat <- rep(0, 6)
  for (ii in 1:6) {
    ss <- Sc[, ii]
    if (min(ss) > 0) {
      q.hat[ii] <- min(q.mat[max(which(ss == min(ss))), , ii]) - 1
    } else {
      if (sum(ss[-length(const.seq)] != 0 & ss[-1] == 0)) {
        q.hat[ii] <- q.mat[which(ss[-length(const.seq)] != 0 & ss[-1] == 0)[1] + 1, 10, ii] - 1
      } else {
        q.hat[ii] <- min(q.mat[max(which(ss == 0)), , ii]) - 1
      }
    }
  }

  if (do.plot) {
    par(mfrow = c(2, 3))
    for (ii in 1:6) {
      plot(const.seq, q.mat[, 10, ii] - 1, type = "b", pch = 1, col = 2, bty = "n", axes = FALSE, xlab = "constant", ylab = "", main = paste("IC ", ii))
      box()
      axis(1, at = pretty(range(const.seq)))
      axis(2, at = pretty(range(q.mat[, 10, ii] - 1)), col = 2, col.ticks = 2, col.axis = 2)
      par(new = TRUE)
      plot(const.seq, Sc[, ii], col = 4, pch = 2, type = "b", bty = "n", axes = FALSE, xlab = "", ylab = "")
      axis(4, at = pretty(range(Sc[, ii])), col = 4, col.ticks = 4, col.axis = 4)
      legend("topright", legend = c("q", "Sc"), col = c(2, 4), lty = c(1, 1), pch = c(1, 2), bty = "n")
    }
  }

  ls <- list(q.hat = q.hat, Gamma_x = Gamma_x, Sigma_x = Sigma_x, sv = sv)
  return(ls)
}

#' @title Forecasting by fnets
#' @method predict fnets
#' @description Produces forecasts of the data for a given forecasting horizon by
#' separately estimating the best linear predictors of common and idiosyncratic components
#' @param object \code{fnets} object
#' @param x input time series matrix, with each row representing a variable
#' @param h forecasting horizon
#' @param common.method a string specifying the method for common component forecasting; possible values are:
#' \itemize{
#'    \item{\code{"restricted"}}{ performs forecasting under a restricted factor model}
#'    \item{\code{"unrestricted"}}{ performs forecasting under an unrestrictive, blockwise VAR representation of the common component}
#' }
#' @param r number of restricted factors; if \code{common.method = "restricted"} and \code{r = NULL},
#' it is estimated as the maximiser of the ratio of the successive eigenvalues of the estimate of the common component covariance matrix,
#' see Ahn and Horenstein (2013)
#' @param ... not used
#' @return a list containing
#' \item{forecast}{ forecasts for the given forecasting horizon}
#' \item{common.pred}{ a list containing forecasting results for the common component}
#' \item{idio.pred}{ a list containing forecasting results for the idiosyncratic component}
#' \item{mean.x}{ \code{mean.x} argument from \code{object}}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @seealso \link[fnets]{fnets}, \link[fnets]{common.predict}, \link[fnets]{idio.predict}
#' @export
predict.fnets <- function(object, x, h = 1, common.method = c("restricted", "unrestricted"), r = NULL, ...) {
  cpre <- common.predict(object, x, h, common.method, r)
  ipre <- idio.predict(object, x, cpre, h)

  out <- list(
    forecast = cpre$fc + ipre$fc,
    common.pred = cpre, idio.pred = ipre,
    mean.x = object$mean.x
  )
  return(out)
}

#' @title Plotting the networks estimated by fnets
#' @method plot fnets
#' @description Plotting method for S3 objects of class \code{fnets}.
#' Produces a plot visualising three networks underlying factor-adjusted VAR processes:
#' (i) directed network representing Granger causal linkages, as given by estimated VAR transition matrices summed across the lags,
#' (ii) undirected network representing contemporaneous linkages after accounting for lead-lag dependence, as given by partial correlations of VAR innovations,
#' (iii) undirected network summarising (i) and (ii) as given by long-run partial correlations of VAR processes.
#' @details See Barigozzi, Cho and Owens (2021) for further details.
#' @param x \code{fnets} object
#' @param type a string specifying which of the above three networks (i)--(iii) to visualise; possible values are
#' \itemize{
#'    \item{\code{"granger"}}{ directed network representing Granger causal linkages}
#'    \item{\code{"pc"}}{ undirected network representing contemporaneous linkages; available when \code{x$lrpc.method = "par"}}
#'    \item{\code{"lrpc"}}{ undirected network summarising Granger causal and contemporaneous linkages; available when \code{x$lrpc.method = "par"} or \code{x$lrpc.method = "npar"}}
#' }
#' @param display a string specifying how to visualise the network; possible values are:
#' \itemize{
#'    \item{\code{"network"}}{ as an \code{igraph} object, see \link[igraph]{plot.igraph}}
#'    \item{\code{"heatmap"}}{ as a heatmap, see \link[fields]{imagePlot}}
#' }
#' @param names a character vector containing the names of the vertices
#' @param groups an integer vector denoting any group structure of the vertices
#' @param threshold if \code{threshold > 0}, hard thresholding is performed on the matrix giving rise to the network of interest
#' @param ... additional arguments
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @seealso \link[fnets]{fnets}
#' @import igraph
#' @importFrom fields imagePlot
#' @importFrom grDevices rainbow
#' @importFrom graphics mtext axis
#' @importFrom RColorBrewer brewer.pal
#' @export
plot.fnets <- function(x, type = c("granger", "pc", "lrpc"), display = c("network", "heatmap"),
                       names = NA, groups = NA, threshold = 0, ...) {
  type <- match.arg(type, c("granger", "pc", "lrpc"))
  display <- match.arg(display, c("network", "heatmap"))

  p <- dim(x$acv$Gamma_x)[1]
  A <- matrix(0, nrow = p, ncol = p)

  if(is.null(x$idio.var)){
    warning(paste0("object contains no idiosyncratic component"))
  } else {
    if (type == "granger") {
      d <- dim(x$idio.var$beta)[1] / p
      for (ll in 1:d) A <- A + t(x$idio.var$beta)[, (ll - 1) * p + 1:p]
      nm <- "Granger causal"
    }

    if (type == "pc") {
      if (x$lrpc.method != "par") {
        stop(paste0("Partial correlation matrix is undetected"))
      } else {
        A <- x$lrpc$pc
        nm <- "Partial correlation"
      }
    }

    if (type == "lrpc") {
      if (!(x$lrpc.method %in% c("par", "npar"))) {
        stop(paste0("Long-run partial correlation matrix is undetected"))
      } else {
        A <- x$lrpc$lrpc
        nm <- "Long-run partial correlation"
      }
    }
    nm <- paste(nm, display, sep = " ")

    A[abs(A) < threshold] <- 0

    if (!is.na(groups[1])) {
      grps <- perm <- c()
      K <- length(unique(groups))
      for (ii in 1:K) {
        permii <- which(groups == unique(groups)[ii])
        perm <- c(perm, permii)
        grps <- c(grps, rep(ii, length(permii)))
      }
    } else {
      perm <- 1:p
      grps <- rep(1, p)
      K <- 1
    }
    grp.col <- rep(rainbow(K, alpha = 1), table(grps))
    A <- A[perm, perm]
    if (!is.na(names[1])) names <- names[perm]

    if (display == "network") {
      v.col <- rep(rainbow(K, alpha = .2), table(grps))
      if (type == "granger") g <- igraph::graph_from_adjacency_matrix(A, mode = "directed", weighted = TRUE, diag = FALSE, ...)
      if (type %in% c("pc", "lrpc")) g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE, ...)
      lg <- igraph::layout_in_circle(g)
      igraph::plot.igraph(g,
                          main = nm, layout = lg, vertex.label = names, vertex.label.font = 2,
                          vertex.shape = "circle", vertex.color = v.col,
                          vertex.label.color = grp.col, vertex.label.cex = 0.6,
                          edge.color = "gray40", edge.arrow.size = 0.5
      )
    } else if (display == "heatmap") {
      heat.cols <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
      if (type == "granger") mv <- max(1e-3, abs(A))
      if (type %in% c("pc", "lrpc")) {
        A[abs(A) > 1] <- sign(A[abs(A) > 1])
        diag(A) <- 0
        mv <- 1.01
      }
      breaks <- seq(-mv, mv, length.out = 12)

      fields::imagePlot(A,
                        axes = FALSE, col = heat.cols,
                        breaks = breaks, main = nm, ...
      )
      if (!is.na(names[1]) || !is.na(groups[1])) {
        if (is.na(names[1])) names <- groups[perm]
        for (ii in 1:p) mtext(text = names[ii], at = (ii - 1) / (p - 1), side = 1, las = 2, cex = .8, col = grp.col[ii])
        for (ii in 1:p) mtext(text = names[ii], at = (ii - 1) / (p - 1), side = 2, las = 2, cex = .8, col = grp.col[ii])
      }
    }
  }
}
