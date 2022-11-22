#' @title Factor number estimator of Bai and Ng (2002)
#' @description Estimates the number of restricted factors by minimising an information criterion.
#' Currently three information criteria proposed in Owens, Cho, and Barigozzi (2022) (\code{pen.op = 1,2,3}) are implemented,
#' with \code{pen.op = 2} recommended as a default choice based on numerical experiments.
#' @details See Bai and Ng (2002) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param covx covariance of \code{x}
#' @param q.max maximum number of factors; if \code{q.max = NULL}, a default value is selected as \code{min(50, floor(sqrt(min(dim(x)[2] - 1, dim(x)[1]))))}
#' @param do.plot whether to plot the value of the information criterion
#' @param center whether to de-mean the input \code{x} row-wise
#' @return a list containing
#' \item{q.hat}{ the mimimiser of the chosen information criteria}
#' \item{sv}{ svd of \code{covx}}
#' @example R/examples/bn_ex.R
#' @references preprint
#' @references Bai, J. & Ng, S. (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221. 
#' @references Alessi, L., Barigozzi, M.,  & Capasso, M. (2010) Improved penalization for determining the number of factors in approximate factor models. Statistics & Probability Letters, 80(23-24):1806–1813.
#' @importFrom graphics abline
#' @keywords internal
bn.factor.number <- function(x, covx = NULL, q.max = NULL,  do.plot = FALSE, center = TRUE) {
  p <- dim(x)[1]
  n <- dim(x)[2]
  cnt <- min(n, p)
  if (center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  if (is.null(q.max)) q.max <- min(50, floor(sqrt(min(n - 1, p))))

  if(is.null(covx)) covx <- xx %*% t(xx) / n

  p.seq <- floor(3 * p / 4 + (1:10) * p / 40)
  n.seq <- n - (9:0) * floor(n / 20)
  const.seq <- seq(.001, 2, by = .01)
  IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 3))


  for (kk in 1:10) {
    nn <- n.seq[kk]
    pp <- p.seq[kk]
    pen <- c(
      (n + p) / (n * p) * log(n * p / (n + p)),
      (n + p) / (n * p) * log(cnt),
      log(cnt) / cnt
    )

    mm <- floor(n / 4) - 1
    #sv <- list(1:(mm + 1))
    tmp <- rep(0, q.max + 1)
    if (kk == length(n.seq)) nu <- q.max else nu <- 0
    sv <- svd(covx[1:pp, 1:pp], nu = nu, nv = 0)
    dd <- sum(sv$d)
    tmp[1] <- tmp[1] + dd / pp #/ (2 * mm + 1)
    for (jj in 1:q.max) {
      dd <- dd - sv$d[jj]
      tmp[jj + 1] <- tmp[jj + 1] + dd / pp #/ (2 * mm + 1)
    }
    for (jj in 1:length(const.seq)) {
      for (pen.op in 1:3) {
        IC[, jj, kk, pen.op] <- tmp + (0:q.max) * const.seq[jj] * pen[pen.op]
      }
    }

  }
  #

  q.mat <- apply(IC, c(2, 3, 4), which.min)
  Sc <- apply(q.mat, c(1, 3), var)
  q.hat <- rep(0, 3)
  for (ii in 1:3) {
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
    par(mfrow = c(1, 3))
    for (ii in 1:3) {
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
  return(list(q.hat = q.hat, sv = sv))
}




#' @title Static PCA
#' @keywords internal
static.pca <- function(xx, q = NULL, q.method = c("ic","er"), q.max = NULL, pen.op = 2, kern.bw = NULL, mm = NULL) {
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  cnt <- min(n, p)
  if(is.null(kern.bw)) kern.bw <- 4 * floor((n/log(n))^(1/3))
  if (is.null(mm)) mm <- min(max(1, kern.bw), floor(n / 4) - 1) else mm <- min(max(mm, 1, kern.bw), floor(n / 4) - 1)


  if (is.null(q.max)) q.max <- min(50, floor(sqrt(min(n - 1, p))))
  q.method <- match.arg(q.method, c("ic", "er"))
  if (is.null(pen.op)) pen.op <- 2

  covx <- xx %*% t(xx) / n
  eig <- eigen(covx, symmetric = TRUE)
  lam <- eig$vectors[, 1:(cnt - 1), drop = FALSE] * sqrt(p)
  f <- t(xx) %*% (eig$vectors[, 1:(cnt - 1), drop = FALSE]) / sqrt(p)


  if(is.null(q)){
    if(q.method == "er") q <- which.max(eig$values[1:q.max] / eig$values[1 + 1:q.max])
    if(q.method == "ic") {
      bn <- bn.factor.number(xx, covx = covx, q.max = q.max, do.plot = FALSE, center = FALSE)
      q <- bn$q.hat[pen.op]
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
#' @param forecast.restricted whether to forecast using a restricted or unrestricted, blockwise VAR representation of the common component
#' @param r number of restricted factors, or a string specifying the factor number selection method when \code{forecast.restricted = TRUE};
#'  possible values are:
#' \itemize{
#'    \item{\code{"ic"}}{ information criteria of Bai and Ng (2002)}
#'    \item{\code{"er"}}{ eigenvalue ratio}
#' }
#' @param ... not used
#' @return a list containing
#' \item{is}{ in-sample predictions}
#' \item{forecast}{ forecasts for the given forecasting horizon}
#' \item{r}{ factor number}
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @seealso \link[fnets]{fnets.factor.model}, \link[fnets]{common.predict}
#' @export
predict.fm <- function(object, x, h = 1, forecast.restricted = TRUE, r = c("ic","er"), ...) {
  out <- common.predict(object = object, x = x, h = h, forecast.restricted = forecast.restricted, r = r)
  return(out)
}
