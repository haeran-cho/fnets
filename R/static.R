




#' @title Static PCA
#' @keywords internal
static.pca <-
  function(xx,
           q = NULL,
           q.method = c("ic", "er"),
           q.max = NULL,
           pen.op = 2,
           kern.bw = NULL,
           mm = NULL) {
    p <- dim(xx)[1]
    n <- dim(xx)[2]
    cnt <- min(n, p)
    if (is.null(kern.bw))
      kern.bw <- 4 * floor((n / log(n)) ^ (1 / 3))
    if (is.null(mm))
      mm <-
      min(max(1, kern.bw), floor(n / 4) - 1)
    else
      mm <- min(max(mm, 1, kern.bw), floor(n / 4) - 1)


    if (is.null(q.max))
      q.max <- min(50, floor(sqrt(min(n - 1, p))))
    q.method <- match.arg(q.method, c("ic", "er"))
    if (is.null(pen.op))
      pen.op <- 2

    covx <- xx %*% t(xx) / n
    eig <- eigen(covx, symmetric = TRUE)
    lam <- eig$vectors[, 1:(cnt - 1), drop = FALSE] * sqrt(p)
    f <-
      t(xx) %*% (eig$vectors[, 1:(cnt - 1), drop = FALSE]) / sqrt(p)


    if (is.null(q)) {
      if (q.method == "er"){
        q.method.out <- eig$values[1:q.max] / eig$values[1 + 1:q.max]
        q <- which.max(q.method.out)
      }

      if (q.method == "ic") {
        q.method.out <-
          abc.factor.number(
            xx,
            covx = covx,
            q.max = q.max,
            do.plot = FALSE,
            center = FALSE
          )
        q <- q.method.out$q.hat[pen.op]
      }
    }



    # q <- as.integer(q); hl <- NA
    proj <-
      eig$vectors[, 1:q, drop = FALSE] %*% t(eig$vectors[, 1:q, drop = FALSE])
    Gamma_c <- Gamma_x <- array(0, dim = c(p, p, 2 * mm + 1))
    for (h in 0:(mm - 1)) {
      Gamma_x[, , h + 1] <-
        xx[, 1:(mm - h)] %*% t(xx[, 1:(mm - h) + h]) / n
      Gamma_c[, , h + 1] <- proj %*% Gamma_x[, , h + 1] %*% proj
      if (h != 0) {
        Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_c[, , 2 * mm + 1 - h + 1] <- t(Gamma_c[, , h + 1])
      }
    }


    Gamma_i <- Gamma_x - Gamma_c
    acv <-
      list(
        Gamma_x = Gamma_x,
        Gamma_c = Re(Gamma_c),
        Gamma_i = Re(Gamma_i)
      )
    return(list(
      q = q,
      lam = lam,
      f = f,
      acv = acv,
      q.method.out = q.method.out
    ))
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
#'    \item{\code{"ic"}}{ information criteria of Alessi, Barigozzi & Capasso (2010)}
#'    \item{\code{"er"}}{ eigenvalue ratio}
#' }
#' @param ... not used
#' @return a list containing
#' \item{is}{ in-sample predictions}
#' \item{forecast}{ forecasts for the given forecasting horizon}
#' \item{r}{ factor number}
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @references Alessi, L., Barigozzi, M.,  & Capasso, M. (2010) Improved penalization for determining the number of factors in approximate factor models. Statistics & Probability Letters, 80(23-24):1806–1813.
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022)
#' @seealso \link[fnets]{fnets.factor.model}, \link[fnets]{common.predict}
#' @export
predict.fm <-
  function(object,
           x,
           h = 1,
           forecast.restricted = TRUE,
           r = c("ic", "er"),
           ...) {
    out <-
      common.predict(
        object = object,
        x = x,
        h = h,
        forecast.restricted = forecast.restricted,
        r = r
      )
    return(out)
  }
