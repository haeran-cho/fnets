



#' @title Factor model estimation
#' @description Unrestricted and restricted factor model estimation
#' @details See Barigozzi, Cho and Owens (2021) for further details.
#'
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param factor.model a string specifying the method to be adopted for factor model estimation; possible values are:
#' \itemize{
#'    \item{\code{"unrestricted"}}{ unrestricted factor model}
#'    \item{\code{"restricted"}}{ restricted factor model}
#' }
#' @param q number of factors. If \code{q = NULL}, the factor number is estimated by an information criterion-based approach of Hallin and Liška (2007) or Bai and Ng (2002), see \link[fnets]{hl.factor.number} and \link[fnets]{bn.factor.number} for further details
#' @param q.method a string specifying the factor number selection method when \code{factor.model = "restricted"}; possible values are:
#' \itemize{
#'    \item{\code{"bn"}}{ information criteria of Bai and Ng (2002)}
#'    \item{\code{"er"}}{ eigenvalue ratio}
#' }
#' @param ic.op choice of the information criterion, see \link[fnets]{hl.factor.number} or \link[fnets]{bn.factor.number} for further details
#' @param kern.const constant multiplied to \code{floor((dim(x)[2]/log(dim(x)[2]))^(1/3)))} which determines the kernel bandwidth for dynamic PCA
#' @param common.args a list specifying the tuning parameters required for estimating the impulse response functions and common shocks. It contains:
#' \itemize{
#'    \item{\code{var.order}}{ order of the blockwise VAR representation of the common component. If \code{var.order = NULL}, it is selected blockwise by Schwarz criterion}
#'    \item{\code{max.var.order}}{ maximum blockwise VAR order for the Schwarz criterion}
#'    \item{\code{trunc.lags}}{ truncation lag for impulse response function estimation}
#'    \item{\code{n.perm}}{ number of cross-sectional permutations involved in impulse response function estimation}
#' }
#' @return an S3 object of class \code{fm}, which contains the following fields:
#' \item{q}{ number of factors}
#' \item{spec}{ if \code{factor.model = "unrestricted"} a list containing estimates of the spectral density matrices for \code{x}, common and idiosyncratic components}
#' \item{acv}{ a list containing estimates of the autocovariance matrices for \code{x}, common and idiosyncratic components}
#' \item{common.irf}{ if if \code{factor.model = "unrestricted"} and \code{q >= 1}, a list containing estimators of the impulse response functions (as an array of dimension \code{(p, q, trunc.lags + 2)})
#' and common shocks (an array of dimension \code{(q, n)}) for the common component}
#' \item{lam}{ if \code{factor.model = "restricted"} factor loadings}
#' \item{f}{ if \code{factor.model = "restricted"} factor series}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @references Bai, J. & Ng, S. (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221. \cr
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.restricted(n, p)
#' x <- common$data
#' out <- fnets.factor.model(x, factor.model = "restricted")
#' }
#' @export
fnets.factor.model <- function(x, center = TRUE, factor.model = c("unrestricted", "restricted"), q = NULL, q.method = c("bn","er"),
                               ic.op = NULL, kern.const = 4,
                               common.args = list(var.order = NULL, max.var.order = NULL, trunc.lags = 20, n.perm = 10)) {
  p <- dim(x)[1]
  n <- dim(x)[2]

  if (center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  common.args <- check.list.arg(common.args)

  q.method <- match.arg(q.method, c("bn", "er"))
  factor.model <- match.arg(factor.model, c("unrestricted", "restricted"))
  if (factor.model == "restricted") {
    spca <- static.pca(xx, q = q, q.method = q.method, ic.op = ic.op, kern.const = kern.const)
    q <- spca$q
    lam <- spca$lam
    f <- spca$f
    spec <- NULL
    acv <- spca$acv
    cve <- NULL
  }
  # dynamic pca
  if (factor.model == "unrestricted") {
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



  out <- list(
    q = q, spec = spec, lam = lam, f = f, acv = acv,
    common.irf = cve, mean.x = mean.x
  )
  attr(out, "factor") <- factor.model
  attr(out, "class") <- "fm" #"fnets"

  return(out)
}


#' @title Factor number estimator of Bai and Ng (2002)
#' @description Estimates the number of restricted factors by minimising an information criterion.
#' Currently the five information criteria proposed in Bai and Ng (2002) (\code{ic.op = 1,...,5}) are implemented,
#' with \code{ic.op = 2} recommended as a default choice based on numerical experiments.
#' @details See Bai and Ng (2002) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param lam,f loading and factor matrices; if \code{lam = NULL} or \code{f = NULL}, these are obtained with PCA
#' @param q.max maximum number of factors; if \code{q.max = NULL}, a default value is selected as \code{min(50, floor(sqrt(min(dim(x)[2] - 1, dim(x)[1]))))}
#' @param ic.op chosen information criterion
#' @param do.plot whether to plot the value of the information criterion
#' @param center whether to de-mean the input \code{x} row-wise
#' @return a list containing
#' \item{q.hat}{ the mimimiser of the chosen information criteria}
#' \item{lam}{ loading matrix}
#' \item{f}{ factor series}
#' \item{q.max}{ maximum number of factors}
#' \item{ic}{ vector of information criteria values}
#' @example R/examples/bn_ex.R
#' @references Bai, J. & Ng, S. (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221. \cr
#' @importFrom graphics abline matplot
#' @export
bn.factor.number <- function(x, lam = NULL, f = NULL, q.max = NULL, ic.op = 2, do.plot = FALSE, center = TRUE) {
  p <- dim(x)[1]
  n <- dim(x)[2]
  cnt <- min(n, p)
  if (center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  if (is.null(q.max)) q.max <- min(50, floor(sqrt(min(n - 1, p))))

  if (is.null(lam) | is.null(f)) {
    covx <- xx %*% t(xx) / n
    eig <- eigen(covx, symmetric = TRUE)
    lam <- eig$vectors[, 1:(cnt - 1), drop = FALSE] * sqrt(p)
    f <- t(xx) %*% (eig$vectors[, 1:(cnt - 1), drop = FALSE]) / sqrt(p)
  }

  ic <- rep(0, 1 + q.max)
  ic[1] <- (ic.op <= 4) * log(mean(xx^2)) + (ic.op == 5) * mean(xx^2)
  l <- 1
  while (l <= q.max) {
    hchi <- lam[, 1:l, drop = FALSE] %*% t(f[, 1:l, drop = FALSE]) # lam[, 1:l, drop=FALSE]%*%f[1:l, , drop=FALSE]
    ic[l + 1] <- (ic.op <= 4) * log(mean((xx - hchi)^2)) +
      (ic.op == 1) * l * (n + p) / (n * p) * log(n * p / (n + p)) +
      (ic.op == 2) * l * (n + p) / (n * p) * log(cnt) +
      (ic.op == 3) * l * log(cnt) / cnt +
      (ic.op == 4) * l * ((n + p - l) * log(n * p) / (n * p) + (n + p) / (n * p) * log(cnt)) / 2 +
      (ic.op == 5) * (mean((xx - hchi)^2) + l * mean((xx - hchi)^2) * (n + p - l) * log(n * p) / (n * p))
    l <- l + 1
  }
  q.hat <- which(ic == min(ic)) - 1
  if (do.plot) {
    par(mfrow = c(1, 1))
    matplot(0:q.max, ic,
      type = "p", col = 2, pch = 2,
      xlab = "Factor number", ylab = "IC", main = "IC for factor number selection"
    )
    abline(v = q.hat)
  }
  return(list(q.hat = q.hat, lam = lam, f = f, q.max = q.max, ic = ic))
}




#' @title Static PCA
#' @keywords internal
static.pca <- function(xx, q.max = NULL, q = NULL, q.method = c("bn","er"), ic.op = 2, kern.const = 4, mm = NULL) {
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  cnt <- min(n, p)
  if (is.null(mm)) mm <- min(max(1, kern.const * floor((n / log(n))^(1 / 3))), floor(n / 4) - 1) else mm <- min(max(mm, 1, kern.const * floor((n / log(n))^(1 / 3))), floor(n / 4) - 1)


  if (is.null(q.max)) q.max <- min(50, floor(sqrt(min(n - 1, p))))
  q.method <- match.arg(q.method, c("bn", "er"))
  if (is.null(ic.op)) ic.op <- 2

  covx <- xx %*% t(xx) / n
  eig <- eigen(covx, symmetric = TRUE)
  lam <- eig$vectors[, 1:(cnt - 1), drop = FALSE] * sqrt(p)
  f <- t(xx) %*% (eig$vectors[, 1:(cnt - 1), drop = FALSE]) / sqrt(p)

  bn <- NULL
  if (is.null(q)) {
    if(q.method == "er") q <- which.max(eig$values[1:q.max] / eig$values[1 + 1:q.max])
    if(q.method == "bn") {
      bn <- bn.factor.number(xx, lam = lam, f = f, q.max = q.max, ic.op = ic.op, do.plot = FALSE, center = FALSE)
      q <- bn$q.hat
    }
  }

  # q <- as.integer(q); hl <- NA
  proj <- eig$vectors[, 1:q, drop = FALSE] %*% t(eig$vectors[, 1:q, drop = FALSE])
  Gamma_c <- Gamma_x <- array(0, dim = c(p, p, 2 * mm + 1))
  for (h in 0:(mm - 1)) {
    Gamma_x[, , h + 1] <- xx[, 1:(mm - h)] %*% t(xx[, 1:(mm - h) + h]) / n
    Gamma_c[, , h + 1] <- proj %*% Gamma_x[, , h + 1] %*% proj
    if (h != 0) {
      Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
      Gamma_c[, , 2 * mm + 1 - h + 1] <- t(Gamma_c[, , h + 1])
    }
  }


  Gamma_i <- Gamma_x - Gamma_c
  acv <- list(Gamma_x = Gamma_x, Gamma_c = Re(Gamma_c), Gamma_i = Re(Gamma_i))
  return(list(q = q, lam = lam, f = f, acv = acv))
}


#' @title Forecasting for factor models
#' @method predict fm
#' @description Produces forecasts of the data for a given forecasting horizon by
#'estimating the best linear predictors of the common component
#' @param object \code{fm} object
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
#' @seealso \link[fnets]{fnets.factor.model}, \link[fnets]{common.predict}
#' @export
predict.fm <- function(object, x, h = 1, common.method = c("restricted", "unrestricted"), r = NULL, ...) {
  common.method <- match.arg(common.method, c("unrestricted", "restricted"))
  out <- common.predict(object = object, x = x, h = h, common.method = common.method, r = r)
  return(out)
}
