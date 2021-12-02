library(igraph)

#' @title Factor-adjusted network analysis
#' @description This function estimates the spectral density and autocovariance matrices of the common and the idiosyncratic components, impulse response function and common shocks, and (sparse) VAR transition matrix and innovation covariance matrix.
#' @details
#'  Further information can be found in Barigozzi, Cho and Owens (2021).
#'
#' @param x input time series matrix, with each row representing a time series
#' @param q the number of factors, if q=NULL this is selected by the information criterion-based estimator of Hallin and Liska (2007)
#' @param ic.op an index number for the information criterion
#' @param common.var.args A list specifying the estimator for the common component. This contains:
#' \itemize{
#'    \item{\code{'var.order'}}{the order of the VAR model, if NULL then selected blockwise by BIC}
#'    \item{\code{'max.var.order'}}{the maximum order of the VAR model for the BIC to consider}
#'    \item{\code{'trunc.lags'}}{the order of the MA representation}
#'    \item{\code{'n.perm'}}{number of cross-sectional permutations}
#' }
#' @param idio.method A string specifying the type of l1-regularised estimator for the idiosyncratic VAR matrix, possible values are:
#' \itemize{
#'    \item{\code{'lasso'}}{Lasso estimator}
#'    \item{\code{'ds'}}{Dantzig Selector}
#' }
#' @param idio.cv.args  A list specifying arguments to the cross-validation (CV) procedure for the idiosyncratic VAR. This contains:
#' \itemize{
#'    \item{\code{'n.folds'}}{number of folds}
#'    \item{\code{'path.length'}}{number of lambda values to consider}
#'    \item{\code{'symmetric'}}{symmetrisation method for Gamma matrix}
#'    \item{\code{'cv.plot'}}{Boolean selecting whether to plot the CV curve}
#' }
#' @param center demean the input `x`
#' @return An S3 object of class \code{fnets}, which contains the following fields:
#' \itemize{
#' \item{q: see Arguments}
#' \item{spec: Spectral density matrices}
#' \item{acv: Autocovariance matrices}
#' \item{common.var: Estimated common component}
#' \item{idio.var: Estimated idiosyncratic component}
#' \item{mean.x: Removed means of x}
#' \item{kern.bandwidth.const: Constant to determine bandwidth size}
#' }
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @examples fnets(sample.data, q=2, idio.method = "lasso")
#' @export
fnets <- function(x, q = NULL, ic.op = 4, kern.bandwidth.const = 4,
                           common.var.args = list(var.order = 1, max.var.order = NULL, trunc.lags = 20, n.perm = 10),
                           idio.var.order = 1, idio.method = c('ds', 'lasso'),
                           idio.cv.args = list(n.folds = 1, path.length = 10, symmetric = 'min', cv.plot = TRUE),
                           center = TRUE){
  p <- dim(x)[1]
  n <- dim(x)[2]

  idio.method <- match.arg(idio.method, c('lasso', 'ds'))
  if(center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  # dynamic pca
  dpca <- dyn.pca(xx, q, ic.op, kern.bandwidth.const)
  q <- dpca$q
  spec <- dpca$spec
  acv <- dpca$acv

  ## common estimation
  cve <- common.var.estimation(xx, Gamma_c = acv$Gamma_c, q = q,
                               var.order = common.var.args$var.order, max.var.order = common.var.args$max.var.order,
                               trunc.lags = common.var.args$trunc.lags, n.perm = common.var.args$n.perm)

  ## idio estimation
  mg <- make.gg(acv$Gamma_i, idio.var.order)
  gg <- mg$gg; GG <- mg$GG

  icv <- idio.cv(xx, lambda.max = max(abs(GG)), var.order = idio.var.order, idio.method = idio.method,
          path.length = idio.cv.args$path.length, n.folds = idio.cv.args$n.folds,
          q = q, kern.bandwidth.const = kern.bandwidth.const, cv.plot = idio.cv.args$n.folds)
  if(idio.method == 'lasso') ive <- var.lasso(GG, gg, lambda = icv$lambda, symmetric = idio.cv.args$symmetric)
  if(idio.method == 'ds') ive <- var.dantzig(GG, gg, lambda = icv$lambda, symmetric = idio.cv.args$symmetric)

  out <- list(q = q, spec = spec, acv = acv,
              common.var = cve, idio.var = ive, mean.x = mean.x,
              kern.bandwidth.const = kern.bandwidth.const)
  attr(out, 'class') <- 'fnets'
  return(out)

}

#' @title Dynamic PCA
#' @description Performs principal components analysis of the autocovariance matrices.
#' @param xx centred input time series matrix, with each row representing a time series
#' @param q the number of factors, if q=NULL this is selected by the information criterion-based estimator of Hallin and Liska (2007)
#' @param ic.op an index number for the information criterion
#' @param kern.bandwidth.const constant to determine bandwidth size
#' @return A list containing
#' \itemize{
#' \item{q: see Arguments}
#' \item{spec: Spectral density matrices}
#' \item{acv: Autocovariance matrices}
#' \item{kern.bandwidth.const: Constant to determine bandwidth size}
#' }
#' @example dyn.pca(sample.data, q=2)
#' @export
dyn.pca <- function(xx, q = NULL, ic.op = 4, kern.bandwidth.const = 4){

  p <- dim(xx)[1]
  n <- dim(xx)[2]

  mm <- min(max(1, kern.bandwidth.const * floor((n/log(n))^(1/3))), floor(n/4) - 1)
  len <- 2 * mm
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- weights.Bartlett(((-mm):mm)/mm)

  # dynamic pca

  if(!is.null(q)){
    q <- as.integer(q)
    Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * mm + 1))
    for(h in 0:(mm - 1)){
      Gamma_x[, , h + 1] <- xx[, 1:(n - h)] %*% t(xx[, 1:(n - h) + h])/n
      Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + mm + 1]
      if(h != 0){
        Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_xw[, , 2 * mm + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
      }
    }
    Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)
    sv <- list(1:(mm + 1))
    if(q > 0) for(ii in 1:(mm + 1)) sv[[ii]] <- svd(Sigma_x[, , ii], nu = q, nv = 0)
  }
  if(is.null(q)){
    q.max <- min(50, floor(sqrt(min(n - 1, p))))
    qq <- hl.factor.number(xx, q.max, mm, w)
    q <- qq$q.hat[ic.op]
    Gamma_x <- qq$Gamma_x
    Sigma_x <- qq$Sigma_x
    sv <- qq$sv
  }

  Gamma_c <- Gamma_i <- Sigma_c <- Sigma_i <- Sigma_x * 0
  if(q >= 1){
    for(ii in 1:(mm + 1)){
      Sigma_c[,, ii] <- sv[[ii]]$u[, 1:q, drop = FALSE] %*% diag(sv[[ii]]$d[1:q], q) %*% Conj(t(sv[[ii]]$u[, 1:q, drop = FALSE]))
      if(ii > 1){
        Sigma_c[,, 2 * mm + 1 - (ii - 1) + 1] <- Conj(Sigma_c[,, ii])
      }
    }
    Gamma_c <-  (aperm(apply(Sigma_c, c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) ) * (2 * pi) / (2 * mm + 1)
    Gamma_c <- Re(Gamma_c)
  }
  # max(abs(Im(Gamma_c)))
  # Gamma_i <-  (aperm(apply(Sigma_i, c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) ) * (2 * pi) / (2 * mm + 1)
  Sigma_i <- Sigma_x - Sigma_c
  Gamma_i <- Gamma_x - Gamma_c

  spec <- list(Sigma_x = Sigma_x, Sigma_c = Sigma_c, Sigma_i = Sigma_i)
  acv <- list(Gamma_x = Gamma_x, Gamma_c = Re(Gamma_c), Gamma_i = Re(Gamma_i))

  out <- list(q = q, spec = spec, acv = acv, kern.bandwidth.const = kern.bandwidth.const)
  return(out)

}

#' @title Factor number estimator of Hallin and Liska (2011)
#' @description Selects the factor number `q` based on 6 information criteria.
#' @param xx centred input time series matrix, with each row representing a time series
#' @param q.max the maximum number of factors to consider
#' @param mm bandwidth
#' @param w weights
#' @return A list containing
#' \itemize{
#' \item{q.hat: Estimated factor numbers corresponding to each criterion}
#' \item{Gamma_x: Autocovariance of x}
#' \item{Sigma_x: Spectral density of x}
#' \item{sv: singular value decomposition of Sigma_x}
#' }
#' @references Hallin, M., & Liška, R. (2007). Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @export
hl.factor.number <- function(xx, q.max, mm, w){

  p <- dim(xx)[1]; n <- dim(xx)[2]

  p.seq <- floor(3*p/4 + (1:10) * p/40)
  n.seq <- n - (9:0) * floor(n/20)
  const.seq <- seq(.001, 2, by = .01)
  IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 2 * 3))

  Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * mm + 1))

  for(kk in 1:10){
    nn <- n.seq[kk]
    pp <- p.seq[kk]
    pen <- c((1/mm^2 + sqrt(mm/nn) + 1/pp) * log(min(pp, mm^2, sqrt(nn/mm))),
             1/sqrt(min(pp, mm^2, sqrt(nn/mm))),
             1/min(pp, mm^2, sqrt(nn/mm)) * log(min(pp, mm^2, sqrt(nn/mm))))

    for(h in 0:(mm - 1)){
      Gamma_x[, , h + 1] <- xx[, 1:(nn - h)] %*% t(xx[, 1:(nn - h) + h])/nn
      Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + mm + 1]
      if(h != 0){
        Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_xw[, , 2 * mm + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
      }
    }
    Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)
    sv <- list(1:(mm + 1))

    tmp <- rep(0, q.max + 1)
    for(ii in 1:(mm + 1)){
      if(kk == length(n.seq)) nu <- q.max else nu <- 0
      sv[[ii]] <- svd(Sigma_x[1:pp, 1:pp, ii], nu = nu, nv = 0)
      dd <- sum(sv[[ii]]$d)
      tmp[1] <- tmp[1] + dd/pp/(2 * mm + 1)
      for(jj in 1:q.max) {
        dd <- dd - sv[[ii]]$d[jj]
        tmp[jj + 1] <- tmp[jj + 1] + dd/pp/(2 * mm + 1)
      }
      for(jj in 1:length(const.seq)){
        for(pen.op in 1:3){
          IC[, jj, kk, 3 * 0 + pen.op] <- tmp + (0:q.max) * const.seq[jj] * pen[pen.op]
          IC[, jj, kk, 3 * 1 + pen.op] <- log(tmp) + (0:q.max) * const.seq[jj] * pen[pen.op]
        }
      }
    }
  }

  q.mat <- apply(IC, c(2, 3, 4), which.min)
  Sc <- apply(q.mat, c(1, 3), var)
  q.hat <- rep(0, 6)
  for(ii in 1:6){
    ss <- Sc[, ii]
    if(min(ss) > 0){
      q.hat[ii] <- q.mat[which.min(ss), ii] - 1
    } else{
      q.hat[ii] <- q.mat[which(ss[-length(const.seq)] != 0 & ss[-1] == 0)[1] + 1, 10, ii] - 1
    }
  }

  par(mfrow = c(2, 3))
  for(ii in 1:6){
    plot(const.seq, q.mat[, 10, ii] - 1, type = 'b', pch = 1, col = 2, bty = 'n', axes = FALSE, xlab = 'constant', ylab = '', main = paste('IC ', ii))
    box()
    axis(1, at = pretty(range(const.seq)))
    axis(2, at = pretty(range(q.mat[, 10, ii] - 1)), col = 2, col.ticks = 2, col.axis = 2)
    par(new = TRUE)
    plot(const.seq, Sc[, ii], col = 4, pch = 2, type = 'b', bty = 'n', axes = FALSE, xlab = '', ylab = '')
    axis(4, at = pretty(range(Sc[, ii])), col = 4, col.ticks = 4, col.axis = 4)
    legend('topright', legend = c('q', 'Sc'), col = c(2, 4), lty = c(1, 1), pch = c(1, 2), bty = 'n')
  }

  ls <- list(q.hat = q.hat, Gamma_x = Gamma_x, Sigma_x = Sigma_x, sv = sv)
  return(ls)

}

#' @title Prediction by fnets
#' @method predict fnets
#' @description Predicts common and idiosyncratic components from a fnets object for new data
#' @param object fnets object
#' @param x input time series matrix, with each row representing a time series
#' @param h forecast horizon
#' @param common.method which of "static" or "var" to forecast the common component with
#' @param r factor number, if r=NULL this is selected using the maximal eigenratio
#' @return A list containing
#' \itemize{
#' \item{fitted: x in-sample estimation}
#' \item{forecast: x forecast}
#' \item{common.pred: Prediction for the factor-driven common component}
#' \item{idio.pred: Prediction for the idiosyncratic component}
#' \item{x.mean: removed mean of x}
#' }
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @export
predict.fnets <- function(object, x, h = 1, common.method = c('static', 'var'), r = NULL){

  cpre <- common.predict(object, x, h, common.method, r)
  ipre <- idio.predict(object, x, cpre, h)

  out <- list(fitted = cpre$is + ipre$is, forecast = cpre$fc + ipre$fc,
              common.pred = cpre, idio.pred = ipre, x.mean = object$mean.x)
  return(out)

}

#' @title Plot fnets object
#' @method plot fnets
#' @references Barigozzi, M., Cho, H., & Owens, D. (2021) Factor-adjusted network analysis for high-dimensional time series.
#' @export
plot.fnets <- function(object, names = NULL, groups = NULL, threshold = 0, size = NULL, ...){
  A <- abs(t(object$idio.var$beta)) ## absolute value graph
  A[A < threshold] <- 0
  p <- dim(A)[1]
  d <- dim(A)[2]/dim(A)[1]

  mark.groups <- list()
  perm <- c()
  if(!is.null(groups)){
    for (ii in unique(groups)) {
      permii <- which(groups == ii)
      mark.groups[[ii]] <- permii
      perm <- c(perm, permii)
    }
   # perm <- Matrix::invPerm(perm)
  } else perm <- 1:p
  # Granger causal network



  granger.mat <- matrix(0, p, p)
  for(ll in 1:d) granger.mat <- granger.mat + A[, (ll - 1) * p + 1:p]
  granger.mat <- granger.mat[perm,perm]#granger <- permute(granger, perm) #
  granger <- igraph::graph_from_adjacency_matrix(granger.mat, mode = "directed", weighted=TRUE, diag = FALSE)
  if(!is.null(size)) {degree <- igraph::degree(granger, mode = size, normalized = F)
  degree <- degree/max(degree)*15 }else degree <- rep(15, p)
  if(!is.null(names)) V(granger)$name <- names[perm]
  #if(length(perm) == p) granger <- permute(granger, perm) #granger.mat <- granger.mat[perm,perm]
  l_granger <- igraph::layout_in_circle(granger)



  par(mfrow = c(1,1))
  plot.igraph(granger, main = "Granger", vertex.label = V(granger)$name, layout = l_granger, vertex.label.font = 2, #mark.groups = mark.groups,
              vertex.label.color = "black", edge.color = "gray40", edge.arrow.size = 0.1, vertex.size = degree, #arrow.width = 5,
              vertex.shape ="circle", vertex.color = groups[perm], vertex.label.cex = 0.8) #"light blue"
}



##

#' @export
#' @keywords internal
weights.Bartlett <- function(x) 1 - abs(x)


