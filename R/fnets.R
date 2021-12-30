#' @title Factor-adjusted network estimation
#' @description Operating under factor-adjusted vector autoregressive (VAR) model, 
#' the function estimates the spectral density and autocovariance matrices of the factor-driven common component and the idiosyncratic VAR process,
#' the impulse response functions and common shocks for the common component, 
#' and VAR parameters and innovation covariance matrix for the idiosyncratic component.
#' @details See Barigozzi, Cho and Owens (2021) for further details.
#'
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param q number of factors. If \code{q = NULL}, the factor number is estimated by an information criterion-based approach of Hallin and Liška (2007), see \code{\link[fnets]{hl.factor.number}} for further details
#' @param ic.op choice of the information criterion, see \code{\link[fnets]{hl.factor.number}} for further details
#' @param kern.const constant multiplied to \code{floor((dim(x)[2]/log(dim(x)[2]))^(1/3)))} which determines the kernel bandwidth for dynamic PCA
#' @param common.var.args a list specifying the tuning parameters required for estimating the impulse response functions and common shocks. It contains:
#' \itemize{
#'    \item{var.order}{ order of the blockwise VAR representation of the common component. If \code{var.order = NULL}, it is selected blockwise by Schwarz criterion}
#'    \item{max.var.order}{ maximum blockwise VAR order for the Schwarz criterion}
#'    \item{trunc.lags}{ truncation lag for impulse response function estimation}
#'    \item{n.perm}{ number of cross-sectional permutations involved in impluse response function estimation}
#' }
#' @param idio.var.order order of the idiosyncratic VAR process; if a vector of integers is supplied, the order is chosen via cross validation
#' @param idio.method a string specifying the method to be adopted for idiosyncratic VAR process estimation; possible values are:
#' \itemize{
#'    \item{"lasso"}{ Lasso-type \code{l1}-regularised \code{M}-estimation}
#'    \item{"ds"}{ Dantzig Selector-type constrained \code{l1}-minimisation}
#' }
#' @param lrpc.method a string specifying the type of estimator for long-run partial correlation matrix estimation; possible values are:
#' \itemize{
#'    \item{"param"}{ parametric estimator based on the VAR model assumption}
#'    \item{"nonpar"}{ nonparametric estimator from inverting the long-run covariance matrix of the idiosyncratic component via constrained \code{l1}-minimisation}
#'    \item{"none"}{ do not estimate the long-run partial correlation matrix}
#' }
#' @param cv.args a list specifying arguments for the cross validation procedures 
#' for selecting the tuning parameters involved in VAR parameter and (long-run) partial correlation matrix estimation. It contains:
#' \itemize{
#'    \item{n.folds}{ number of folds}
#'    \item{path.length}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{do.plot}{ whether to plot the output of the cross validation step}
#' }
#' @return an S3 object of class \code{fnets}, which contains the following fields:
#' \itemize{
#' \item{q}{ number of factors}
#' \item{spec}{ a list containing estimates of the spectral density matrices for \code{x}, common and idiosyncratic components}
#' \item{acv}{ a list containing estimates of the autocovariance matrices for \code{x}, common and idiosyncratic components}
#' \item{common.var}{ if \code{q >= 1}, a list containing estimators of the impulse response functions (as an array of dimension \code{(p, q, trunc.lags + 2)}) 
#' and common shocks (an array of dimension \code{(q, n)}) for the common component}
#' \item{idio.var}{ a list containing the following fields:
#' \itemize{
#' \item{beta}{ estimate of VAR parameter matrix; each column contains parameter estimates for the regression model for a given variable}
#' \item{Gamma}{ estimate of the innovation covariance matrix}
#' \item{lambda}{ regularisation parameter}
#' \item{var.order}{ VAR order}
#' }}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#' \item{kern.const}{ input parameter}
#' }
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @example R/examples/fnets.R
#' @seealso \code{\link[fnets]{predict.fnets}}, \code{\link[fnets]{plot.fnets}}
#' @export
fnets <- function(x, center = TRUE, q = NULL, ic.op = 4, kern.const = 4,
                  common.var.args = list(var.order = 1, max.var.order = NULL, trunc.lags = 20, n.perm = 10),
                  idio.var.order = 1, idio.method = c('lasso', 'ds'),
                  lrpc.method = c('param', 'nonpar', 'none'),
                  cv.args = list(n.folds = 1, path.length = 10, do.plot = FALSE)){
  p <- dim(x)[1]
  n <- dim(x)[2]

  idio.method <- match.arg(idio.method, c('lasso', 'ds'))
  lrpc.method <- match.arg(lrpc.method, c('param', 'nonpar', 'none'))
  if(center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  # dynamic pca
  dpca <- dyn.pca(xx, q, ic.op, kern.const)
  q <- dpca$q
  spec <- dpca$spec
  acv <- dpca$acv

  ## common estimation
  cve <- common.var.estimation(xx, Gamma_c = acv$Gamma_c, q = q,
                               var.order = common.var.args$var.order, max.var.order = common.var.args$max.var.order,
                               trunc.lags = common.var.args$trunc.lags, n.perm = common.var.args$n.perm)

  ## idio estimation

  if(cv.args$do.plot) par(mfrow = c(1, 1 + lrpc.method %in% c('param', 'nonpar')))
    
  icv <- yw.cv(xx, lambda.max = NULL, var.order = idio.var.order, method = idio.method,
               path.length = cv.args$path.length, n.folds = cv.args$n.folds,
               q = q, kern.const = kern.const, do.plot = cv.args$do.plot)
  mg <- make.gg(acv$Gamma_i, icv$var.order)
  gg <- mg$gg; GG <- mg$GG
  if(idio.method == 'lasso') ive <- var.lasso(GG, gg, lambda = icv$lambda, symmetric = 'min')
  if(idio.method == 'ds') ive <- var.dantzig(GG, gg, lambda = icv$lambda, symmetric = 'min')
  ive$var.order <- icv$var.order
  
  out <- list(q = q, spec = spec, acv = acv,
              common.var = cve, idio.var = ive, mean.x = mean.x,
              kern.const = kern.const)
  attr(out, 'class') <- 'fnets'
  
  ## lrpc estimation
  if(lrpc.method %in% c('param', 'nonpar')){
    
    
    out$lrpc <- lrpc
  }
  
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
#' @return a list containing
#' \itemize{
#' \item{q}{ number of factors}
#' \item{spec}{ a list containing the estimates of the spectral density matrices for \code{x}, common and idiosyncratic components}
#' \item{acv}{ a list containing estimates of the autocovariance matrices for \code{x}, common and idiosyncratic components}
#' \item{kern.const}{ input parameter}
#' }
#' @keywords internal
dyn.pca <- function(xx, q = NULL, ic.op = 4, kern.const = 4){

  p <- dim(xx)[1]
  n <- dim(xx)[2]

  mm <- min(max(1, kern.const * floor((n/log(n))^(1/3))), floor(n/4) - 1)
  len <- 2 * mm
  w <- Bartlett.weights(((-mm):mm)/mm)

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
    qq <- hl.factor.number(xx, q.max, mm, w, center = FALSE)
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
  Sigma_i <- Sigma_x - Sigma_c
  Gamma_i <- Gamma_x - Gamma_c

  spec <- list(Sigma_x = Sigma_x, Sigma_c = Sigma_c, Sigma_i = Sigma_i)
  acv <- list(Gamma_x = Gamma_x, Gamma_c = Re(Gamma_c), Gamma_i = Re(Gamma_i))

  out <- list(q = q, spec = spec, acv = acv, kern.const = kern.const)
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
#' @param do.plot whether to produce a plot of six information criteria values
#' @param center whether to de-mean the input \code{x} row-wise
#' @return a list containing
#' \itemize{
#' \item{q.hat}{ a vector containing minimisers of the six information criteria}
#' \item{Gamma_x}{ an array containing the estimates of the autocovariance matrices of \code{x} at \code{2 * mm + 1} lags}
#' \item{Sigma_x}{ an array containing the estimates of the spectral density matrices of \code{x} at \code{2 * mm + 1} Fourier frequencies}
#' \item{sv}{ a list containing the singular value decomposition of \code{Sigma_x}}
#' }
#' @example R/examples/hlfactornumber.R
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @export
hl.factor.number <- function(x, q.max = NULL, mm, w = NULL, do.plot = TRUE, center = TRUE){
  p <- dim(x)[1]; n <- dim(x)[2]
  q.max <- min(50, floor(sqrt(min(n - 1, p))))
  
  if(center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  if(is.null(w)) w <- Bartlett.weights(((-mm):mm)/mm)

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
  
  if(do.plot){
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
#'    \item{"restricted"}{ performs forecasting under a restrictive static factor model}
#'    \item{"unrestricted"}{ performs forecasting under an unrestrictive, blockwise VAR representation of the common component}
#' @param r number of static factors; if \code{common.method = "restricted"} and \code{r = NULL}, 
#' it is estimated as the maximiser of the ratio of the successive eigenvalues of the estimate of the common component covariance matrix,
#' see Ahn and Horenstein (2013)
#' @param ... further arguments
#' @return a list containing
#' \itemize{
#' \item{forecast}{ forecasts for the given forecasting horizon}
#' \item{common.pred}{ a list containing forecasting results for the common component}
#' \item{idio.pred}{ a list containing forecasting results for the idiosyncratic component}
#' \item{mean.x}{ \code{mean.x} argument from \code{object}}
#' }
#' @example R/examples/predict.R
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @export
predict.fnets <- function(object, x, h = 1, common.method = c('restricted', 'unrestricted'), r = NULL, ...){

  cpre <- common.predict(object, x, h, common.method, r)
  ipre <- idio.predict(object, x, cpre, h)

  out <- list(forecast = cpre$fc + ipre$fc,
              common.pred = cpre, idio.pred = ipre, 
              mean.x = object$mean.x)
  return(out)

}

#' @title Plotting the output from network estimation via fnets
#' @method plot fnets
#' @description Plotting method for S3 objects of class \code{fnets}. 
#' Produces a plot visualising the Granger causal network which is determined by
#' aggregating the estimated VAR transition matrices across the lags.
#' @param x \code{fnets} object
#' @param display a string specifying which to be plotted as visualisation of the Granger network underpinning the idiosyncratic VAR process; possible values are:
#' \itemize{
#'    \item{"network"}{ an \code{igraph} object}
#'    \item{"heatmap"}{ performs forecasting under an unrestrictive, blockwise VAR representation of the common component}
#' @param names a character vector containing the names of the vertices
#' @param groups an integer vector denoting any group structure of the vertices
#' @param threshold if \code{threshold > 0} hard thresholding is applied to the aggregated VAR transition matrix before plotting
#' @param size a string specifying the type of degree to be used when \code{display = "network"}; possible values are \code{"all"}, \code{"out"}, \code{"in"} or \code{"total"}
#' @param ... additional arguments
#' @references Barigozzi, M., Cho, H. & Owens, D. (2021) FNETS: Factor-adjusted network analysis for high-dimensional time series.
#' @example R/examples/plot.R
#' @export
plot.fnets <- function(x, display = "network", names = NULL, groups = NULL, threshold = 0, size = NULL, ...){
  
  A <- abs(t(x$idio.var$beta)) 
  A[A < threshold] <- 0
  p <- dim(A)[1]
  d <- dim(A)[2]/dim(A)[1]

  if(!is.null(groups)){
    perm <- c()
    for(ii in unique(groups)){
      permii <- which(groups == ii)
      perm <- c(perm, permii)
    }
  } else{
    perm <- 1:p
    groups <- rep(1, p)
  }
  grp.col <- rep(rainbow(max(groups), alpha = 0.2), table(groups))

  # Granger causal network
  granger.mat <- matrix(0, p, p)
  for(ll in 1:d) granger.mat <- granger.mat + A[, (ll - 1) * p + 1:p]
  granger.mat <- granger.mat[perm, perm]

  if(display == "network"){
    granger <- igraph::graph_from_adjacency_matrix(granger.mat, 
                                                   mode = "directed", weighted = TRUE, diag = FALSE)
    if(!is.null(size)) {
      degree <- igraph::degree(granger, mode = size, normalized = FALSE)
      degree <- degree/max(degree) * 15 
    } else degree <- rep(15, p)
    if(!is.null(names)) V(granger)$name <- names[perm]
    l_granger <- igraph::layout_in_circle(granger)
    plot.igraph(granger, main = "Granger network", vertex.label = V(granger)$name, layout = l_granger, vertex.label.font = 2,
              vertex.label.color = "black", edge.color = "gray40", edge.arrow.size = 0.1, vertex.size = degree,
              vertex.shape = "circle", vertex.color = grp.col, vertex.label.cex = 0.8)
  }
  else if(display == "heatmap"){
    mv <- max(abs(granger.mat)) * 1.01
    mv <- max(mv, 1e-4)
    fields::imagePlot(granger.mat, axes = FALSE,
                      col = (RColorBrewer::brewer.pal(9, 'Reds') ),
                      breaks = seq(0, mv, length.out = 10),
                      main = "Granger network")
    if(is.null(names)) names <- 1:p
    for(ii in 1:p) mtext(text = names[perm[ii]], at = (ii - 1)/(p - 1), side = 1, las = 2, cex = .6, col = grp.col[ii])
    for(ii in 1:p) mtext(text = names[perm[ii]], at = (ii - 1)/(p - 1), side = 2, las = 2, cex = .6, col = grp.col[ii])
  }
  else warning(paste("display must be either 'network' or 'heatmap'"))
  
}



##

# #' @export
#' @title Bartlett weights
#' @description internal function
#' @keywords internal
Bartlett.weights <- function(x) 1 - abs(x)


