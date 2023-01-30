#' @title Factor-adjusted network estimation
#' @description Operating under a factor-adjusted vector autoregressive (VAR) model,
#' the function estimates the spectral density and autocovariance matrices of the factor-driven common component and the idiosyncratic VAR process,
#' the impulse response functions and common shocks for the common component,
#' and VAR parameters, innovation covariance matrix and long-run partial correlations for the idiosyncratic component.
#' @details See Barigozzi, Cho and Owens (2022) and Owens, Cho and Barigozzi (2022) for further details.
#' List arguments do not need to be specified with all list components; any missing entries will be filled in with the default argument.
#'
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param fm.restricted whether to estimate a restricted factor model using static PCA
#' @param q Either the number of factors or a string specifying the factor number selection method; possible values are:
#' \itemize{
#'    \item{\code{"ic"}}{ information criteria-based methods of Alessi, Barigozzi & Capasso (2010) when \code{fm.restricted = TRUE} or Hallin and Liška (2007) when \code{fm.restricted = FALSE} modifying Bai and Ng (2002)}
#'    \item{\code{"er"}}{ eigenvalue ratio of Ahn and Horenstein (2013)}
#' }
#' see \link[fnets]{factor.number}.
#' @param ic.op choice of the information criterion penalty, see \link[fnets]{factor.number} for further details
#' @param kern.bw a positive integer specifying the kernel bandwidth for dynamic PCA; by default, it is set to \code{floor(4 *(dim(x)[2]/log(dim(x)[2]))^(1/3)))}.  When \code{fm.restricted = TRUE}, it is used to compute the number of lags for which autocovariance matrices are estimated
#' @param common.args a list specifying the tuning parameters required for estimating the impulse response functions and common shocks. It contains:
#' \itemize{
#'    \item{\code{factor.var.order}}{ order of the blockwise VAR representation of the common component. If \code{factor.var.order = NULL}, it is selected blockwise by Schwarz criterion}
#'    \item{\code{max.var.order}}{ maximum blockwise VAR order for the Schwarz criterion}
#'    \item{\code{trunc.lags}}{ truncation lag for impulse response function estimation}
#'    \item{\code{n.perm}}{ number of cross-sectional permutations involved in impulse response function estimation}
#' }
#' @param var.order order of the idiosyncratic VAR process; if a vector of integers is supplied, the order is chosen via \code{tuning}
#' @param var.method a string specifying the method to be adopted for idiosyncratic VAR process estimation; possible values are:
#' \itemize{
#'        \item{\code{"lasso"}}{ Lasso-type \code{l1}-regularised \code{M}-estimation}
#'        \item{\code{"ds"}}{ Dantzig Selector-type constrained \code{l1}-minimisation}
#' }
#' @param var.args a list specifying the tuning parameters required for estimating the idiosyncratic VAR process. It contains:
#' \itemize{
#'    \item{\code{n.iter}}{ maximum number of descent steps, by default depends on \code{var.order}; applicable when \code{var.method = "lasso"}}
#'    \item{\code{n.cores}}{ number of cores to use for parallel computing, see \link[parallel]{makePSOCKcluster}; applicable when \code{var.method = "ds"}}
#' }
#' @param do.threshold whether to perform adaptive thresholding of all parameter estimators with \link[fnets]{threshold}
#' @param do.lrpc whether to estimate the long-run partial correlation
#' @param lrpc.adaptive whether to use the adaptive estimation procedure
#' @param tuning.args a list specifying arguments for \code{tuning}
#' for selecting the tuning parameters involved in VAR parameter and (long-run) partial correlation matrix estimation. It contains:
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
#' @return an S3 object of class \code{fnets}, which contains the following fields:
#' \item{q}{ number of factors}
#' \item{spec}{ if \code{fm.restricted = FALSE} a list containing estimates of the spectral density matrices for \code{x}, common and idiosyncratic components}
#' \item{acv}{ a list containing estimates of the autocovariance matrices for \code{x}, common and idiosyncratic components}
#' \item{loadings}{ if \code{fm.restricted = TRUE}, factor loadings; if \code{fm.restricted = FALSE} and \code{q >= 1},
#' a list containing estimators of the impulse response functions (as an array of dimension \code{(p, q, trunc.lags + 2)})}
#' \item{factors}{ if \code{fm.restricted = TRUE}, factor series; else, common shocks (an array of dimension \code{(q, n)})}
#' \item{idio.var}{ a list containing the following fields:
#' \itemize{
#' \item{\code{beta}}{ estimate of VAR parameter matrix; each column contains parameter estimates for the regression model for a given variable}
#' \item{\code{Gamma}}{ estimate of the innovation covariance matrix}
#' \item{\code{lambda}}{ regularisation parameter}
#' \item{\code{var.order}}{ VAR order}
#' }}
#' \item{lrpc}{ see the output of \link[fnets]{par.lrpc}}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#' \item{var.method}{ input parameter}
#' \item{do.lrpc}{ input parameter}
#' \item{kern.bw}{ input parameter}
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @references Alessi, L., Barigozzi, M.,  & Capasso, M. (2010) Improved penalization for determining the number of factors in approximate factor models. Statistics & Probability Letters, 80(23-24):1806–1813.
#' @references Bai, J. & Ng, S. (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221.
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. arXiv preprint arXiv:2301.11675.
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.unrestricted(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#' out <- fnets(x,
#'   q = NULL, var.order = 1, var.method = "lasso", do.threshold = TRUE,
#'   do.lrpc = TRUE, tuning.args = list(tuning = "cv", n.folds = 1, path.length = 10, do.plot = TRUE),
#'   var.args = list(n.cores = 2)
#' )
#' pre <- predict(out, x, h = 1, common.method = "unrestricted")
#' plot(out, type = "granger", display = "network")
#' plot(out, type = "lrpc", display = "heatmap")
#' }
#' @seealso \link[fnets]{predict.fnets}, \link[fnets]{plot.fnets}
#' @importFrom graphics par
#' @export
fnets <-
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
           ),
           var.order = 1,
           var.method = c("lasso", "ds"),
           var.args = list(
             n.iter = NULL,
             n.cores = min(parallel::detectCores() - 1, 3)
           ),
           do.threshold = FALSE,
           do.lrpc = TRUE,
           lrpc.adaptive = FALSE,
           tuning.args = list(
             tuning = c("cv", "bic"),
             n.folds = 1,
             penalty = NULL,
             path.length = 10,
             do.plot = FALSE
           )) {
    p <- dim(x)[1]
    n <- dim(x)[2]

    var.args <- check.list.arg(var.args)
    common.args <- check.list.arg(common.args)
    tuning.args <- check.list.arg(tuning.args)

    if(!is.numeric(q)) {
      q.method <- match.arg(q, c("ic", "er"))
      q <- NULL
    } else
      q.method <- NULL

    var.method <- match.arg(var.method, c("lasso", "ds"))
    tuning <- match.arg(tuning.args$tuning, c("cv", "bic"))
    if(center)
      mean.x <- apply(x, 1, mean)
    else
      mean.x <- rep(0, p)
    xx <- x - mean.x

    if(fm.restricted) {
      spca <-
        static.pca(
          xx,
          q.max = NULL,
          q = q,
          q.method = q.method,
          ic.op = ic.op
        )
      q <- spca$q
      loadings <- spca$lam
      factors <- spca$f
      spec <- NA
      acv <- spca$acv
      kern.bw <- NA
    } else {
      if(is.null(kern.bw))
        kern.bw <-  floor(4 * (n / log(n))^(1/3))

      ## dynamic pca
      dpca <- dyn.pca(xx, q, q.method, ic.op, kern.bw)
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
      loadings <- cve$irf.est
      factors <- cve$u.est
    }

    ## idio estimation
    if(tuning.args$do.plot){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      par(mfrow = c(1, 1 + do.lrpc))
    }
    if(tuning == "cv") {
      icv <- yw.cv(
        xx,
        method = var.method,
        lambda.max = NULL,
        var.order = var.order,
        n.folds = tuning.args$n.folds,
        path.length = tuning.args$path.length,
        q = q,
        kern.bw = kern.bw,
        do.plot = tuning.args$do.plot
      )
    }

    if(tuning == "bic") {
      icv <- yw.ic(
        xx,
        method = var.method,
        lambda.max = NULL,
        var.order = var.order,
        penalty = tuning.args$penalty,
        path.length = tuning.args$path.length,
        q = q,
        kern.bw = kern.bw,
        do.plot = tuning.args$do.plot
      )
    }

    mg <- make.gg(acv$Gamma_i, icv$var.order)
    gg <- mg$gg
    GG <- mg$GG
    if(var.method == "lasso"){
      if(is.null(var.args$n.iter)) var.args$n.iter <- var.order*100
      ive <-
      var.lasso(
        GG,
        gg,
        lambda = icv$lambda,
        symmetric = "min",
        n.iter = var.args$n.iter
      )
    }
    if(var.method == "ds")
      ive <-
      var.dantzig(
        GG,
        gg,
        lambda = icv$lambda,
        symmetric = "min",
        n.cores = var.args$n.cores
      )
    ive$var.order <- icv$var.order
    if(do.threshold)
      ive$beta <-
      threshold(ive$beta, do.plot = tuning.args$do.plot)$thr.mat

    out <- list(
      q = q,
      spec = spec,
      acv = acv,
      loadings = loadings,
      factors = factors,
      idio.var = ive,
      mean.x = mean.x,
      var.method = var.method,
      do.lrpc = do.lrpc,
      kern.bw = kern.bw
    )

    if(fm.restricted)
      attr(out, "factor") <-
      "restricted"
    else
      attr(out, "factor") <- "unrestricted"

    ## lrpc estimation
    if(do.lrpc) {
      out$lrpc <-
        par.lrpc(
          out,
          x,
          eta = NULL,
          tuning.args = tuning.args,
          do.threshold = do.threshold,
          lrpc.adaptive = lrpc.adaptive,
          n.cores = var.args$n.cores
        )
    } else {
      out$lrpc <- NA
    }

    attr(out, "class") <- "fnets"
    return(out)
  }

#' @title Plotting the networks estimated by fnets
#' @method plot fnets
#' @description Plotting method for S3 objects of class \code{fnets}.
#' Produces a plot visualising three networks underlying factor-adjusted VAR processes:
#' (i) directed network representing Granger causal linkages, as given by estimated VAR transition matrices summed across the lags,
#' (ii) undirected network representing contemporaneous linkages after accounting for lead-lag dependence, as given by partial correlations of VAR innovations,
#' (iii) undirected network summarising (i) and (ii) as given by long-run partial correlations of VAR processes.
#' Edge widths are determined by edge weights.
#' @details See Barigozzi, Cho and Owens (2022) for further details.
#' @param x \code{fnets} object
#' @param type a string specifying which of the above three networks (i)--(iii) to visualise; possible values are
#' \itemize{
#'    \item{\code{"granger"}}{ directed network representing Granger causal linkages}
#'    \item{\code{"pc"}}{ undirected network representing contemporaneous linkages; available when \code{x$do.lrpc = TRUE}}
#'    \item{\code{"lrpc"}}{ undirected network summarising Granger causal and contemporaneous linkages; available when \code{x$do.lrpc = TRUE}}
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
#' @return A network plot produced as per the input arguments
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. arXiv preprint arXiv:2301.11675.
#' @seealso \link[fnets]{fnets}
#' @import igraph
#' @importFrom fields imagePlot
#' @importFrom grDevices rainbow
#' @importFrom graphics mtext axis
#' @importFrom RColorBrewer brewer.pal
#' @export
plot.fnets <-
  function(x,
           type = c("granger", "pc", "lrpc"),
           display = c("network", "heatmap"),
           names = NA,
           groups = NA,
           threshold = 0,
           ...) {
    type <- match.arg(type, c("granger", "pc", "lrpc"))
    display <- match.arg(display, c("network", "heatmap"))

    p <- dim(x$acv$Gamma_x)[1]
    A <- matrix(0, nrow = p, ncol = p)

    if(is.null(x$idio.var)) {
      warning(paste0("object contains no idiosyncratic component"))
    } else {
      if(type == "granger") {
        d <- dim(x$idio.var$beta)[1] / p
        for (ll in 1:d)
          A <- A + t(x$idio.var$beta)[, (ll - 1) * p + 1:p]
        nm <- "Granger causal"
      }

      if(type == "pc") {
        if(!x$do.lrpc & is.null(x$lrpc$pc) ){
          stop(paste0("Partial correlation matrix is undetected"))
        } else {
          A <- x$lrpc$pc
          nm <- "Partial correlation"
        }
      }

      if(type == "lrpc") {
        if(!x$do.lrpc & is.null(x$lrpc$lrpc)) {
          stop(paste0("Long-run partial correlation matrix is undetected"))
        } else {
          A <- x$lrpc$lrpc
          nm <- "Long-run partial correlation"
        }
      }
      nm <- paste(nm, display, sep = " ")

      A[abs(A) < threshold] <- 0

      if(!is.na(groups[1])) {
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
      if(!is.na(names[1]))
        names <- names[perm]

      if(display == "network") {
        v.col <- rep(rainbow(K, alpha = .2), table(grps))
        if(type == "granger")
          g <-
            igraph::graph_from_adjacency_matrix(A,
                                                mode = "directed",
                                                weighted = TRUE,
                                                diag = FALSE,
                                                ...)
        if(type %in% c("pc", "lrpc"))
          g <-
            igraph::graph_from_adjacency_matrix(A,
                                                mode = "undirected",
                                                weighted = TRUE,
                                                diag = FALSE,
                                                ...)
        lg <- igraph::layout_in_circle(g)
        igraph::plot.igraph(
          g,
          main = nm,
          layout = lg,
          vertex.label = names,
          vertex.label.font = 2,
          vertex.shape = "circle",
          vertex.color = v.col,
          vertex.label.color = grp.col,
          vertex.label.cex = 0.6,
          edge.color = "gray40",
          edge.arrow.size = 0.5,
          edge.width = .5 + 3 * igraph::E(g)$weight
        )
      } else if(display == "heatmap") {
        heat.cols <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
        if(type == "granger")
          mv <- max(1e-3, abs(A))
        if(type %in% c("pc", "lrpc")) {
          A[abs(A) > 1] <- sign(A[abs(A) > 1])
          diag(A) <- 0
          mv <- 1.01
        }
        breaks <- seq(-mv, mv, length.out = 12)

        fields::imagePlot(
          A,
          axes = FALSE,
          col = heat.cols,
          breaks = breaks,
          main = nm,
          ...
        )
        if(!is.na(names[1]) || !is.na(groups[1])) {
          if(is.na(names[1]))
            names <- groups[perm]
          for (ii in 1:p)
            mtext(
              text = names[ii],
              at = (ii - 1) / (p - 1),
              side = 1,
              las = 2,
              cex = .8,
              col = grp.col[ii]
            )
          for (ii in 1:p)
            mtext(
              text = names[ii],
              at = (ii - 1) / (p - 1),
              side = 2,
              las = 2,
              cex = .8,
              col = grp.col[ii]
            )
        }
      }
    }
  }

#' @title Forecasting by fnets
#' @method predict fnets
#' @description Produces forecasts of the data for a given forecasting horizon by
#' separately estimating the best linear predictors of common and idiosyncratic components
#' @param object \code{fnets} object
#' @param x input time series matrix, with each row representing a variable
#' @param h forecasting horizon
#' @param fc.restricted whether to forecast using a restricted or unrestricted, blockwise VAR representation of the common component
#' @param r number of restricted factors, or a string specifying the factor number selection method when \code{fc.restricted = TRUE};
#'  possible values are:
#' \itemize{
#'    \item{\code{"ic"}}{ information criteria of Bai and Ng (2002)}
#'    \item{\code{"er"}}{ eigenvalue ratio}
#' }
#' @param ... not used
#' @return a list containing
#' \item{forecast}{ forecasts for the given forecasting horizon}
#' \item{common.pred}{ a list containing forecasting results for the common component}
#' \item{idio.pred}{ a list containing forecasting results for the idiosyncratic component}
#' \item{mean.x}{ \code{mean.x} argument from \code{object}}
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling
#' @seealso \link[fnets]{fnets}, \link[fnets]{common.predict}, \link[fnets]{idio.predict}
#' @examples
#' set.seed(123)
#' n <- 500
#' p <- 50
#' common <- sim.restricted(n, p)
#' idio <- sim.var(n, p)
#' x <- common$data + idio$data
#' out <- fnets(x, q = 2, var.order = 1, var.method = "lasso",
#' do.lrpc = FALSE, var.args = list(n.cores = 2))
#' cpre.unr <- common.predict(out, x, h = 1, fc.restricted = FALSE, r = NULL)
#' cpre.res <- common.predict(out, x, h = 1, fc.restricted = TRUE, r = NULL)
#' ipre <- idio.predict(out, x, cpre.res, h = 1)
#' @export
predict.fnets <-
  function(object,
           x,
           h = 1,
           fc.restricted = TRUE,
           r = c("ic", "er"),
           ...) {
    cpre <- common.predict(object, x, h, fc.restricted, r)
    ipre <- idio.predict(object, x, cpre, h)

    out <- list(
      forecast = cpre$fc + ipre$fc,
      common.pred = cpre,
      idio.pred = ipre,
      mean.x = object$mean.x
    )
    return(out)
  }
