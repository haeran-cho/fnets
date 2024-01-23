#' @title Factor model estimation
#' @description Performs factor modelling under either restricted (static) or unrestricted (dynamic) factor models
#' @details See Barigozzi, Cho and Owens (2024+) for further details.
#'
#' @param x input time series each column representing a time series variable; it is coerced into a \link[stats]{ts} object
#' @param center whether to de-mean the input \code{x}
#' @param fm.restricted whether to estimate a restricted factor model using static PCA
#' @param q Either a string specifying the factor number selection method when \code{fm.restricted = TRUE}; possible values are:
#' \describe{
#'    \item{\code{"ic"}}{ information criteria-based methods of Alessi, Barigozzi & Capasso (2010) when \code{fm.restricted = TRUE} or Hallin and Liška (2007) when \code{fm.restricted = FALSE}}
#'    \item{\code{"er"}}{ eigenvalue ratio of Ahn and Horenstein (2013) when \code{fm.restricted = TRUE} or Avarucci et al. (2022) when \code{fm.restricted = FALSE}}
#' }
#' or the number of unrestricted factors, see \link[fnets]{factor.number}
#' @param ic.op choice of the information criterion penalty, see \link[fnets]{hl.factor.number} or \link[fnets]{abc.factor.number} for further details
#' @param kern.bw a positive integer specifying the kernel bandwidth for dynamic PCA;
#' by default, it is set to \code{floor(4 *(dim(x)[2]/log(dim(x)[2]))^(1/3)))}.
#' When \code{fm.restricted = TRUE}, it is used to compute the number of lags for which autocovariance matrices are estimated
#' @param common.args a list specifying the tuning parameters required for estimating the impulse response functions and common shocks. It contains:
#' \describe{
#'    \item{\code{factor.var.order}}{ order of the blockwise VAR representation of the common component. If \code{factor.var.order = NULL}, it is selected blockwise by Schwarz criterion}
#'    \item{\code{max.var.order}}{ maximum blockwise VAR order for the Schwarz criterion}
#'    \item{\code{trunc.lags}}{ truncation lag for impulse response function estimation}
#'    \item{\code{n.perm}}{ number of cross-sectional permutations involved in impulse response function estimation}
#' }
#' @return an S3 object of class \code{fm}, which contains the following fields:
#' \item{q}{ number of factors}
#' \item{spec}{ if \code{fm.restricted = FALSE} a list containing estimates of the spectral density matrices for \code{x}, common and idiosyncratic components}
#' \item{acv}{ a list containing estimates of the autocovariance matrices for \code{x}, common and idiosyncratic components}
#' \item{loadings}{ if \code{fm.restricted = TRUE}, factor loadings; if \code{fm.restricted = FALSE} and \code{q >= 1},
#' a list containing estimators of the impulse response functions (as an array of dimension \code{(p, q, trunc.lags + 2)})}
#' \item{factors}{ if \code{fm.restricted = TRUE}, factor series; else, common shocks (an array of dimension \code{(q, n)})}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @references Alessi, L., Barigozzi, M.,  & Capasso, M. (2010) Improved penalization for determining the number of factors in approximate factor models. Statistics & Probability Letters, 80(23-24):1806–1813.
#' @references Avarucci, M., Cavicchioli, M., Forni, M., & Zaffaroni, P. (2022) The main business cycle shock(s): Frequency-band estimation of the number of dynamic factors.
#' @references Barigozzi, M., Cho, H. & Owens, D. (2024+) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. Journal of Business & Economic Statistics (to appear).
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2024+) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. The R Journal (to appear).
#' @seealso \link[fnets]{print.fm}, \link[fnets]{predict.fm}
#' @examples
#' \donttest{
#' out <- fnets.factor.model(data.restricted, fm.restricted = TRUE)
#' }
#' @importFrom stats as.ts
#' @export
fnets.factor.model <-
  function(x,
           center = TRUE,
           fm.restricted = FALSE,
           q = c("ic", "er"),
           ic.op = NULL,
           kern.bw = NULL,
           common.args = list(
             factor.var.order = NULL,
             max.var.order = NULL,
             trunc.lags = 20,
             n.perm = 10
           )) {
    x <- t(as.ts(x))
    p <- dim(x)[1]
    n <- dim(x)[2]

    if(!is.null(ic.op)) ic.op <- max(1, min(as.integer(ic.op), 6))

    ifelse(center, mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
    xx <- x - mean.x

    if(!is.numeric(q)) {
      q.method <- match.arg(q, c("ic", "er"))
      q <- NULL
    } else {
      q <- posint(q, 0)
      q.method <- NULL
    }

    args <- as.list(environment())
    args$x <- t(args$x)

    common.args <- check.list.arg(common.args)
    ifelse(is.null(kern.bw),
           kern.bw <- 4 * floor((n / log(n))^(1/3)),
           kern.bw <- max(0, kern.bw))
    mm <- min(max(1, kern.bw), floor(n / 4) - 1)

    #if(!is.null(q.method)) q.method <- match.arg(q.method, c("ic", "er"))
    if(fm.restricted) {
      if(isTRUE(q==0)){
        spca <- list()
      } else {
        spca <-
          static.pca(
            xx,
            q = q,
            q.method = q.method,
            ic.op = ic.op,
            mm = mm
          )
        q <- spca$q
      }
      loadings <- spca$lam
      factors <- spca$f
      spec <- NULL
      acv <- spca$acv
    } else {
      dpca <- dyn.pca(xx, q, q.method, ic.op, kern.bw, mm)
      q <- dpca$q
      spec <- dpca$spec
      acv <- dpca$acv
      ## common VAR estimation
      cve <- common.irf.estimation(
        xx,
        Gamma_c = acv$Gamma_c,
        q = q,
        factor.var.order = common.args$factor.var.order,
        max.var.order = common.args$max.var.order,
        trunc.lags = common.args$trunc.lags,
        n.perm = common.args$n.perm
      )
      if (q >= 1) {
        loadings <- cve$irf.est
        factors <-  cve$u.est
      } else loadings <- factors <- NULL
    }
    out <- list(
      q = q,
      spec = spec,
      loadings = loadings,
      factors = factors,
      acv = acv,
      mean.x = mean.x
    )
    ifelse(fm.restricted,attr(out, "factor") <-"restricted",
           attr(out, "factor") <- "unrestricted")
    attr(out, "class") <- "fm"
    attr(out, "args") <- args
    return(out)
}

#' @title Dynamic PCA
#' @description Performs principal components analysis in frequency domain for identifying common and idiosyncratic components.
#' @param xx centred input time series matrix, with each row representing a variable
#' @param q number of factors. If \code{q = NULL}, the factor number is estimated by an information criterion-based approach of Hallin and Liška (2007)
#' @param q.method A string specifying the factor number selection method; possible values are:
#' \describe{
#'    \item{\code{"ic"}}{ information criteria-based methods of Alessi, Barigozzi & Capasso (2010) when \code{fm.restricted = TRUE} or Hallin and Liška (2007) when \code{fm.restricted = FALSE}}
#'    \item{\code{"er"}}{ eigenvalue ratio of Ahn and Horenstein (2013)}
#' }
#' @param ic.op choice of the information criterion penalty. Currently the three options from Hallin and Liška (2007) (\code{ic.op = 1, 2} or \code{3}) and
#' their variations with logarithm taken on the cost (\code{ic.op = 4, 5} or \code{6}) are implemented,
#' with \code{ic.op = 5} recommended as a default choice based on numerical experiments
#' @param kern.bw a positive integer specifying the kernel bandwidth for dynamic PCA; by default, it is set to \code{floor(4 *(dim(x)[2]/log(dim(x)[2]))^(1/3)))}
#' @param mm bandwidth
#' @return a list containing
#' \item{q}{ number of factors}
#' \item{q.method.out}{ if \code{q = NULL}, the output from the chosen \code{q.method}, either a vector of eigenvalue ratios or \link[fnets]{hl.factor.number}}
#' \item{spec}{ a list containing the estimates of the spectral density matrices for \code{x}, common and idiosyncratic components}
#' \item{acv}{ a list containing estimates of the autocovariance matrices for \code{x}, common and idiosyncratic components}
#' \item{kern.bw}{ input parameter}
#' @importFrom stats fft
#' @keywords internal
dyn.pca <-
  function(xx,
           q = NULL,
           q.method = c("ic", "er"),
           ic.op = 5,
           kern.bw = NULL,
           mm = NULL) {
    p <- dim(xx)[1]
    n <- dim(xx)[2]

    q.method <- match.arg(q.method, c("ic", "er"))
    if(is.null(ic.op))
      ic.op <- 5
    if(is.null(kern.bw)) kern.bw <-  floor(4 * (n / log(n))^(1/3))
    mm <- min(max(1, kern.bw), floor(n / 4) - 1)

    len <- 2 * mm
    w <- Bartlett.weights(((-mm):mm) / mm)

    ## dynamic pca
    if(is.null(q)) {
      q <- min(50, floor(sqrt(min(n - 1, p))))
      flag <- TRUE
    } else flag <- FALSE
    q <- as.integer(q)
    Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * mm + 1))
    for (h in 0:(mm - 1)) {
      Gamma_x[, , h + 1] <- xx[, 1:(n - h)] %*% t(xx[, 1:(n - h) + h]) / n
      Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + mm + 1]
      if(h != 0) {
        Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_xw[, , 2 * mm + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
      }
    }
    Sigma_x <-
      aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)
    sv <- list(1:(mm + 1))
    if(q > 0)
      for (ii in 1:(mm + 1))
        sv[[ii]] <- svd(Sigma_x[, , ii], nu = q, nv = 0)


    if(flag){
      if(q.method == "er") {
        eigs <- rep(0, q + 1)
        for (ii in 1:(mm + 1)){
          eigs <- eigs + sv[[ii]]$d[0:q + 1]
        }
        q.method.out <- eigs[1:q] / eigs[1 + 1:q]
        q <- which.max(q.method.out)
      } else if(q.method == "ic") {
        q.max <- min(50, floor(sqrt(min(n - 1, p))))
        q.method.out <-
          hl.factor.number(xx, q.max, mm, center = FALSE)
        q <- q.method.out$q.hat[ic.op]
      }
    }

    Gamma_c <- Gamma_i <- Sigma_c <- Sigma_i <- Sigma_x * 0
    if(q >= 1) {
      for (ii in 1:(mm + 1)) {
        Sigma_c[, , ii] <-
          sv[[ii]]$u[, 1:q, drop = FALSE] %*% diag(sv[[ii]]$d[1:q], q) %*% Conj(t(sv[[ii]]$u[, 1:q, drop = FALSE]))
        if(ii > 1) {
          Sigma_c[, , 2 * mm + 1 - (ii - 1) + 1] <- Conj(Sigma_c[, , ii])
        }
      }
      Gamma_c <-
        aperm(apply(Sigma_c, c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) * (2 * pi) / (2 * mm + 1)
      Gamma_c <- Re(Gamma_c)
    }
    Sigma_i <- Sigma_x - Sigma_c
    Gamma_i <- Gamma_x - Gamma_c

    spec <-
      list(Sigma_x = Sigma_x,
           Sigma_c = Sigma_c,
           Sigma_i = Sigma_i)
    acv <-
      list(
        Gamma_x = Gamma_x,
        Gamma_c = Re(Gamma_c),
        Gamma_i = Re(Gamma_i)
      )

    out <-
      list(
        q = q,
        spec = spec,
        acv = acv,
        kern.bw = kern.bw
      )
    if(q.method == "er" & flag) out$q.method.out <- q.method.out
    return(out)
  }

#' @title Static PCA
#' @keywords internal
static.pca <-
  function(xx,
           q = NULL,
           q.method = c("ic", "er"),
           q.max = NULL,
           ic.op = 2,
           mm = NULL){

    p <- dim(xx)[1]
    n <- dim(xx)[2]

    if(is.null(q.max)) q.max <- min(50, floor(sqrt(min(n - 1, p))))
    q.method <- match.arg(q.method, c("ic", "er"))

    if(is.null(ic.op))
      ic.op <- 5

    if(is.null(mm)) mm <- floor(min(n, p) / log(max(n, p)))

    covx <- xx %*% t(xx) / n
    eig <- eigen(covx, symmetric = TRUE)

    q.method.out <- 0
    if(is.null(q)) {
      if(q.method == "er") {
        q.method.out <- eig$values[1:q.max] / eig$values[1 + 1:q.max]
        q <- which.max(q.method.out)
      }

      if(q.method == "ic") {
        q.method.out <-
          abc.factor.number(
            xx,
            covx = covx,
            q.max = q.max,
            center = FALSE
          )
        q <- q.method.out$q.hat[ic.op]
      }
    }

    lam <- eig$vectors[, 1:q, drop = FALSE] * sqrt(p)
    f <- t(xx) %*% (eig$vectors[, 1:q, drop = FALSE]) / sqrt(p)

    proj <-
      eig$vectors[, 1:q, drop = FALSE] %*% t(eig$vectors[, 1:q, drop = FALSE])
    Gamma_c <- Gamma_x <- array(0, dim = c(p, p, 2 * mm + 1))
    for (h in 0:(mm - 1)) {
      Gamma_x[, , h + 1] <-
        xx[, 1:(mm - h)] %*% t(xx[, 1:(mm - h) + h]) / n
      Gamma_c[, , h + 1] <- proj %*% Gamma_x[, , h + 1] %*% proj
      if(h != 0) {
        Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_c[, , 2 * mm + 1 - h + 1] <- t(Gamma_c[, , h + 1])
      }
    }

    Gamma_i <- Gamma_x - Gamma_c
    acv <-
      list(
        Gamma_x = Gamma_x,
        Gamma_c = Gamma_c,
        Gamma_i = Gamma_i
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
#' @description Produces forecasts of the data input to \code{object} for a given forecasting horizon by
#'estimating the best linear predictors of the common component
#' @param object \code{fm} object
#' @param n.ahead forecasting horizon
#' @param fc.restricted if \code{fc.restricted = TRUE}, the forecast is generated under a restricted factor model
#' @param r number of static factors, or a string specifying the factor number selection method when \code{fc.restricted = TRUE};
#'possible values are:
#' \describe{
#'    \item{\code{"ic"}}{ information criteria of Alessi, Barigozzi & Capasso (2010)}
#'    \item{\code{"er"}}{ eigenvalue ratio of Ahn & Horenstein (2013) }
#' }
#' @param ... not used
#' @return a list containing
#' \item{is}{ in-sample predictions}
#' \item{forecast}{ forecasts for the given forecasting horizon}
#' \item{r}{ factor number}
#' @seealso \link[fnets]{fnets.factor.model}
#' @examples
#' out <- fnets.factor.model(data.restricted, fm.restricted = TRUE)
#' pre <- predict(out)
#' @export
predict.fm <-
  function(object,
           n.ahead = 1,
           fc.restricted = TRUE,
           r = c("ic", "er"),
           ...) {
    x <- t(attr(object, "args")$x)
    n.ahead <- posint(n.ahead)
    out <-
      common.predict(
        object = object,
        x = x,
        n.ahead = n.ahead,
        fc.restricted = fc.restricted,
        r = r
      )
    return(out)
  }

#' @title Print factor model
#' @method print fm
#' @description Prints a summary of a \code{fm} object
#' @param x \code{fm} object
#' @param ... not used
#' @return NULL, printed to console
#' @seealso \link[fnets]{fnets.factor.model}
#' @examples
#' out <- fnets.factor.model(data.restricted, q = "ic")
#' print(out)
#' @export
print.fm <- function(x,
                        ...){
  args <- attr(x, "args")
  cat(paste("Factor model with \n"))
  cat(paste("n: ", args$n, ", p: ", args$p,  "\n", sep = ""))
  cat(paste("Factor model: ", attr(x, "factor"), "\n", sep = ""))
  cat(paste("Number of factors: ", x$q, "\n", sep = ""))
  cat(paste("Factor number selection method: ", ifelse(is.null(args$q.method), "NA", args$q.method), "\n", sep = ""))
  if(!is.null(args$q.method)) if(args$q.method == "ic")
    cat(paste("Information criterion: ", ifelse(is.null(args$ic.op), "IC5", paste("IC", args$ic.op, sep = "")), "\n", sep = ""))
}
