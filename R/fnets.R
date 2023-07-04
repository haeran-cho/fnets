#' @title Factor-adjusted network estimation
#' @description Under a factor-adjusted vector autoregressive (VAR) model,
#' the function estimates the spectral density and autocovariance matrices of the factor-driven common component and the idiosyncratic VAR process,
#' the impulse response functions and common shocks for the common component,
#' and VAR parameters, innovation covariance matrix and long-run partial correlations for the idiosyncratic component.
#' @details See Barigozzi, Cho and Owens (2022) and Owens, Cho and Barigozzi (2022) for further details.
#' List arguments do not need to be specified with all list components; any missing entries will be filled in with the default argument.
#'
#' @param x input time series each column representing a time series variable; it is coerced into a \link[stats]{ts} object
#' @param center whether to de-mean the input \code{x}
#' @param fm.restricted whether to estimate a restricted factor model using static PCA
#' @param q Either the number of factors or a string specifying the factor number selection method; possible values are:
#' \itemize{
#'    \item{\code{"ic"}}{ information criteria-based methods of Alessi, Barigozzi & Capasso (2010) when \code{fm.restricted = TRUE} or Hallin and Liška (2007) when \code{fm.restricted = FALSE}}
#'    \item{\code{"er"}}{ eigenvalue ratio of Ahn and Horenstein (2013) when \code{fm.restricted = TRUE} or Avarucci et al. (2022) when \code{fm.restricted = FALSE}}
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
#' @references Avarucci, M., Cavicchioli, M., Forni, M., & Zaffaroni, P. (2022) The main business cycle shock(s): Frequency-band estimation of the number of dynamic factors.
#' @references Barigozzi, M., Cho, H. & Owens, D. (2022) FNETS: Factor-adjusted network estimation and forecasting for high-dimensional time series. arXiv preprint arXiv:2201.06110.
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. arXiv preprint arXiv:2301.11675.
#' @examples
#' \donttest{
#' out <- fnets(data.unrestricted,
#'   do.threshold = TRUE,
#'   var.args = list(n.cores = 2)
#' )
#' pre <- predict(out, common.method = "unrestricted")
#' plot(out, type = "granger", display = "network")
#' plot(out, type = "lrpc", display = "heatmap")
#' }
#' @seealso \link[fnets]{predict.fnets}, \link[fnets]{plot.fnets}, \link[fnets]{print.fnets}
#' @importFrom stats as.ts
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
             n.cores = 1
           ),
           do.threshold = FALSE,
           do.lrpc = TRUE,
           lrpc.adaptive = FALSE,
           tuning.args = list(
             tuning = c("cv", "bic"),
             n.folds = 1,
             penalty = NULL,
             path.length = 10
           )) {
  x <- t(as.ts(x))
  p <- dim(x)[1]
  n <- dim(x)[2]


  var.args <- check.list.arg(var.args)
  common.args <- check.list.arg(common.args)
  tuning.args <- check.list.arg(tuning.args)

  if(!is.numeric(q)) {
    q.method <- match.arg(q, c("ic", "er"))
    q <- NULL
  } else {
    q <- posint(q, 0)
    q.method <- NULL
  }

  var.method <- match.arg(var.method, c("lasso", "ds"))
  tuning <- match.arg(tuning.args$tuning, c("cv", "bic"))

  args <- as.list(environment())
  args$x <- t(args$x)

  ifelse(center, mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
  xx <- x - mean.x

  if(!fm.restricted & is.null(kern.bw))
    kern.bw <-  floor(4 * (n / log(n))^(1/3))

  fm <- fnets.factor.model(t(xx),
                           center = FALSE,
                           fm.restricted = fm.restricted,
                           q = q,
                           ic.op = ic.op,
                           kern.bw = kern.bw,
                           common.args = common.args)

  if(fm.restricted) {
    spec <- NA
    kern.bw <- NA
  } else {
    spec <- fm$spec
  }
  q <- fm$q
  loadings <- fm$loadings
  factors <- fm$factors
  acv <- fm$acv

  ## idio estimation
  ive <- fnets.var.internal(xx, acv, method = c("lasso", "ds"),
                            lambda = NULL,
                            var.order = var.order,
                            tuning.args = tuning.args,
                            do.threshold = do.threshold,
                            n.iter = var.args$n.iter,
                            tol = var.args$tol,
                            n.cores = var.args$n.cores)
  ive$mean.x <- mean.x

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

  attr(out, "factor") <- ifelse(fm.restricted, "restricted", "unrestricted")
  attr(out, "args") <- args

  ## lrpc estimation
  ifelse(do.lrpc,
         out$lrpc <-
           par.lrpc(
             out,
             eta = NULL,
             tuning.args = tuning.args,
             do.threshold = do.threshold,
             lrpc.adaptive = lrpc.adaptive,
             n.cores = var.args$n.cores
           ),
         out$lrpc <- NA)

  attr(out, "class") <- "fnets"
  #tuning data
  attr(out, "data") <- attr(ive, "data")
  return(out)
}

#' @title internal function for \code{plot.fnets} and \code{network}
#' @keywords internal
plot_internal <- function(x,
                          type = c("granger", "pc", "lrpc"),
                          display = c("network", "heatmap"),
                          names = NA,
                          groups = NA,
                          group.colours = NA,
                          ...){
  is.var <- FALSE
  p <- dim(x$acv$Gamma_x)[1]
  if(is.null(p)) {
    p <- dim(x$beta)[2]
    is.var <- TRUE
  }
  A <- matrix(0, nrow = p, ncol = p)

  if(is.null(x$idio.var) & !is.var) {
    warning("object contains no idiosyncratic component")
  } else {
    if(type == "granger") {
      ifelse(is.var, beta <- x$beta, beta <- x$idio.var$beta)
      d <- dim(beta)[1] /p
      for (ll in 1:d)
        A <- A + t(beta)[, (ll - 1) * p + 1:p]
    }

    if(type == "pc") {
      if(!x$do.lrpc & is.null(x$lrpc$pc) ){
        stop("Partial correlation matrix is undetected")
      } else {
        A <- x$lrpc$pc
      }
    }

    if(type == "lrpc") {
      if(!x$do.lrpc & is.null(x$lrpc$lrpc)) {
        stop("Long-run partial correlation matrix is undetected")
      } else {
        A <- x$lrpc$lrpc
      }
    }

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

    if(is.na(group.colours[1])) group.colours <- grDevices::rainbow(K, alpha = .2 + .8*(display == "heatmap"))
    if(length(group.colours) != K){
      warning("Length of group.colours must be equal to number of groups; default colours will be used")
      group.colours <- grDevices::rainbow(K, alpha = .2 + .8*(display == "heatmap"))
    }
    grp.col <- rep(group.colours, table(grps))

    A <- A[perm, perm]
    if(!is.na(names[1]))
      names <- names[perm]

    return(list(A = A,
                names = names,
                grps = grps,
                grp.col = grp.col,
                K = K,
                perm = perm))
  }
}


#' @title Convert networks into igraph objects
#' @param object object
#' @param ... additional arguments
#' @seealso \link[fnets]{network.fnets}
#' @export
network <- function (object, ...) UseMethod("network", object)


#' @title Convert networks estimated by fnets into igraph objects
#' @method network fnets
#' @exportS3Method fnets::network
#' @description Converts S3 objects of class \code{fnets} into a network.
#' Produces an igraph object for the three networks underlying factor-adjusted VAR processes:
#' (i) directed network representing Granger causal linkages, as given by estimated VAR transition matrices summed across the lags,
#' (ii) undirected network representing contemporaneous linkages after accounting for lead-lag dependence, as given by partial correlations of VAR innovations,
#' (iii) undirected network summarising (i) and (ii) as given by long-run partial correlations of VAR processes.
#' When plotting the network, note that the edge weights may be negative since they correspond to the entries of the estimators of VAR parameters and (long-run) partial correlations.
#' @details See Barigozzi, Cho and Owens (2022) for further details.
#' @param object \code{fnets} object
#' @param type a string specifying which of the above three networks (i)--(iii) to visualise; possible values are
#' \itemize{
#'    \item{\code{"granger"}}{ directed network representing Granger causal linkages}
#'    \item{\code{"pc"}}{ undirected network representing contemporaneous linkages; available when \code{object$do.lrpc = TRUE}}
#'    \item{\code{"lrpc"}}{ undirected network summarising Granger causal and contemporaneous linkages; available when \code{x$do.lrpc = TRUE}}
#' }
#' @param names a character vector containing the names of the vertices
#' @param groups an integer vector denoting any group structure of the vertices
#' @param group.colours a vector denoting colours corresponding to \code{groups}
#' @param ... additional arguments to \code{igraph::graph_from_adjacency_matrix}
#' @return a list containing
#' \item{network}{ \code{igraph} object}
#' \item{names}{ input argument}
#' \item{groups}{ input argument}
#' \item{grp.col}{ vector of colours corresponding to each node}
#' \item{...}{ additional arguments to \code{igraph::graph_from_adjacency_matrix}}
#' @seealso \link[fnets]{fnets}, \link[fnets]{plot.fnets}
#' @examples
#' \donttest{
#' out <- fnets(data.unrestricted,
#'   do.threshold = TRUE,
#'   var.args = list(n.cores = 2)
#' )
#' net <- network(out, type = "granger")$network
#' plot(net, layout = igraph::layout_in_circle(net))
#' network(out, type = "pc")
#' network(out, type = "lrpc")
#' }
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
network.fnets <- function(object,
                          type = c("granger", "pc", "lrpc"),
                          names = NA,
                          groups = NA,
                          group.colours = NA,
                          ...) {

  type <- match.arg(type, c("granger", "pc", "lrpc"))
  int <- plot_internal(object, type, display = "network", names, groups, group.colours, ...)
  A <- int$A

  if(type == "granger")
    g <-
    igraph::graph_from_adjacency_matrix(A,
                                        mode = "directed",
                                        weighted = TRUE,
                                        diag = FALSE,
                                        ...)
  else if(type %in% c("pc", "lrpc"))
    g <-
    igraph::graph_from_adjacency_matrix(A,
                                        mode = "undirected",
                                        weighted = TRUE,
                                        diag = FALSE,
                                        ...)
  if(min(A) < 0)
    warning("Negative edge weights present. Take care when plotting.")
  return(list(network = g,
              names = int$names,
              groups = int$grps,
              grp.col = int$grp.col))
}

#' @title Plotting the networks estimated by fnets
#' @method plot fnets
#' @description Plotting method for S3 objects of class \code{fnets}.
#' When \code{display = "network"} or {"heatmap"}, it produces a plot visualising three networks underlying factor-adjusted VAR processes:
#' (i) directed network representing Granger causal linkages, as given by estimated VAR transition matrices summed across the lags,
#' (ii) undirected network representing contemporaneous linkages after accounting for lead-lag dependence, as given by partial correlations of VAR innovations,
#' (iii) undirected network summarising (i) and (ii) as given by long-run partial correlations of VAR processes.
#' Edge widths are determined by edge weights.
#' When \code{display = "tuning"}, it produces up to two plots (when \code{do.larpc = TRUE}) visualising
#' the outcome of CV or IC adopted for selecting the \code{l1}-regularisation parameters and the VAR order.
#' @details See Barigozzi, Cho and Owens (2022) for further details.
#' @param x \code{fnets} object
#' @param type a string specifying which of the above three networks (i)--(iii) to visualise
#' when \code{display = "network"} or \code{display = "heatmap"}; possible values are
#' \itemize{
#'    \item{\code{"granger"}}{ directed network representing Granger causal linkages}
#'    \item{\code{"pc"}}{ undirected network representing contemporaneous linkages; available when \code{x$do.lrpc = TRUE}}
#'    \item{\code{"lrpc"}}{ undirected network summarising Granger causal and contemporaneous linkages; available when \code{x$do.lrpc = TRUE}}
#' }
#' @param display a string specifying which plot to produce; possible values are
#' \itemize{
#'    \item{\code{"network"}}{ visualise the network as an \code{igraph} object, see \link[igraph]{plot.igraph}}
#'    \item{\code{"heatmap"}}{ visualise the network as a heatmap, see \link[fields]{imagePlot}}
#'    \item{\code{"tuning"}}{ visualise the outcome from CV or IC (specified by \code{tuning.args$tuning} of \link[fnets]{fnets})
#'    for selecting \code{l1}-regularisation parameters and the VAR order}
#' }
#' @param names a character vector containing the names of the network vertices
#' @param groups an integer vector denoting any group structure of the network vertices
#' @param group.colours a vector denoting colours corresponding to \code{groups}
#' @param ... additional arguments
#' @return A plot produced as per the input arguments
#' @seealso \link[fnets]{fnets}
#' @examples
#' \donttest{
#' out <- fnets(data.unrestricted,
#'   do.threshold = TRUE,
#'   var.args = list(n.cores = 2)
#' )
#' plot(out, type = "granger", display = "network",
#' groups = rep(c(1,2), 50/2), group.colours = c("orange","blue"))
#' plot(out, type = "lrpc", display = "heatmap")
#' plot(out, display = "tuning")
#' }
#' @importFrom igraph layout_in_circle plot.igraph E
#' @importFrom fields imagePlot
#' @importFrom grDevices rainbow
#' @importFrom graphics mtext axis
#' @importFrom RColorBrewer brewer.pal
#' @export
plot.fnets <-
    function(x,
             type = c("granger", "pc", "lrpc"),
             display = c("network", "heatmap", "tuning"),
             names = NA,
             groups = NA,
             group.colours = NA,
             ...) {
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))
      type <- match.arg(type, c("granger", "pc", "lrpc"))
      display <- match.arg(display, c("network", "heatmap", "tuning"))

      if(type == "granger"){
        nm <- "Granger causal"
      } else if(type == "pc"){
        nm <- "Partial correlation"
      } else if(type == "lrpc"){
        nm <- "Long-run partial correlation"
      }
      nm <- paste(nm, display, sep = " ")

      if(display == "network") {
        suppressWarnings({
        net <- network(x,
                       type = type,
                       names = names,
                       groups = groups,
                       group.colours = group.colours,
                       ...
                       ) })
        g <- net$network
        names <- net$names
        grp.col <- net$grp.col
        lg <- igraph::layout_in_circle(g)
        if(is.null(igraph::E(g)$weight)) igraph::E(g)$weight <- 0
        igraph::plot.igraph(
          g,
          main = nm,
          layout = lg,
          vertex.label = names,
          vertex.label.font = 2,
          vertex.shape = "circle",
          vertex.color = grp.col,
          vertex.label.color = "black",
          vertex.label.cex = 0.6,
          edge.color = "gray40",
          edge.arrow.size = 0.5,
          edge.width = .5 + 3 * abs(igraph::E(g)$weight)
        )
      } else if(display == "heatmap") {
        p <- attr(x, "args")$p
        int <- plot_internal(x, type, display, names, groups, group.colours, ...)
        A <- int$A
        grp.col <- int$grp.col
        names <- int$names
        groups <- int$grps
        perm <- int$perm

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
        if(!is.na(names[1]) || (!is.na(groups[1]) & !all(groups==1)) ) {
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
      } else if(display == "tuning") {
        tuning_plot(x)
      }
}

#' @title Forecasting by fnets
#' @method predict fnets
#' @description Produces forecasts of the data for a given forecasting horizon by
#' separately estimating the best linear predictors of common and idiosyncratic components
#' @param object \code{fnets} object
#' @param newdata input time series matrix; by default, uses input to \code{object}.
#' Valid only for the case where \code{newdata} is modelled as a VAR process without any factors
#' @param n.ahead forecasting horizon
#' @param fc.restricted whether to forecast using a restricted or unrestricted, blockwise VAR representation of the common component
#' @param r number of static factors, or a string specifying the factor number selection method when \code{fc.restricted = TRUE};
#'  possible values are:
#' \itemize{
#'    \item{\code{"ic"}}{ information criteria of Alessi, Barigozzi & Capasso (2010)}
#'    \item{\code{"er"}}{ eigenvalue ratio of Ahn & Horenstein (2013)}
#' }
#' @param ... not used
#' @return a list containing
#' \item{forecast}{ forecasts for the given forecasting horizon}
#' \item{common.pred}{ a list containing forecasting results for the common component}
#' \item{idio.pred}{ a list containing forecasting results for the idiosyncratic component}
#' \item{mean.x}{ \code{mean.x} argument from \code{object}}
#' @seealso \link[fnets]{fnets}
#' @examples
#' out <- fnets(data.restricted, q = 2, do.lrpc = FALSE, var.args = list(n.cores = 2))
#' pre.unr <- predict(out, fc.restricted = FALSE)
#' pre.res <- predict(out, fc.restricted = TRUE)
#' @importFrom stats as.ts
#' @export
predict.fnets <-
  function(object,
           newdata = NULL,
           n.ahead = 1,
           fc.restricted = TRUE,
           r = c("ic", "er"),
           ...) {
  if(is.null(newdata)) {
    newdata <- attr(object, "args")$x
  } else if (object$q >= 1) {
    stop("To produce forecasts when a common component is present, estimate a model on the new data. \n")
  }
  newdata <- as.ts(newdata)

  n.ahead <- posint(n.ahead)
  if(nrow(newdata) < n.ahead){
    n.ahead <- nrow(newdata)
    warning("Forecast horizon restricted by number of observations")
  }
  cpre <- common.predict(object, t(newdata), n.ahead, fc.restricted, r)
  ipre <- idio.predict(object, t(newdata), cpre, n.ahead)

  out <- list(
    forecast = cpre$fc + ipre$fc,
    common.pred = cpre,
    idio.pred = ipre,
    mean.x = object$mean.x
  )
  return(out)
}

#' @title Print fnets
#' @method print fnets
#' @description Prints a summary of a \code{fnets} object
#' @param x \code{fnets} object
#' @param ... not used
#' @return NULL, printed to console
#' @seealso \link[fnets]{fnets}
#' @examples \donttest{
#' out <- fnets(data.restricted, q = 2,
#' do.lrpc = FALSE, var.args = list(n.cores = 2))
#' print(out)
#' x <- sim.var(500, 50)$data
#' out <- fnets.var(x,
#' n.cores = 2)
#' print(out)
#' }
#' @export
print.fnets <- function(x,
                        ...){
  args <- attr(x, "args")
  if(!is.null(x$acv)) { #fnets.var
    method <- args$var.method
    beta <- x$idio.var$beta
    do.lrpc <- args$do.lrpc
    q <- x$q
  } else {
    method <- args$method
    beta <- x$beta
    do.lrpc <- FALSE
    q <- 0
  }
  cat(paste("Factor-adjusted vector autoregressive model with \n"))
  cat(paste("n: ", args$n, ", p: ", args$p,  "\n", sep = ""))
  cat(paste("Factor-driven common component ---------"), "\n")
  cat(paste("Factor model: ", attr(x, "factor"), "\n", sep = ""))
  cat(paste("Factor number: ", q, "\n", sep = ""))
  cat(paste("Factor number selection method: ", ifelse(is.null(args$q.method), "NA", args$q.method), "\n", sep = ""))
  if(!is.null(args$q.method)) if(args$q.method == "ic")
    cat(paste("Information criterion: ", ifelse(is.null(args$ic.op), "IC5", paste("IC", args$ic.op, sep = "")), "\n", sep = ""))
  cat(paste("Idiosyncratic VAR component ---------"), "\n")
  cat(paste("VAR order: ", args$var.order, "\n", sep = ""))
  cat(paste("VAR estimation method: ", method, "\n", sep = ""))
  cat(paste("Tuning method: ", args$tuning, "\n", sep = ""))
  cat(paste("Threshold: ", args$do.threshold, "\n", sep = ""))
  cat(paste("Non-zero entries: ", sum(beta != 0), "/", prod(dim(beta)), "\n", sep = ""))
  cat(paste("Long-run partial correlations ---------"), "\n")
  cat(paste("LRPC: ", do.lrpc, "\n", sep = ""))
}

