#' @title Factor number selection methods
#' @description Methods to estimate the number of factor.
#' When \code{method = 'er'}, the factor number is estimated by maximising the ration of successive eigenvalues.
#' When \code{method = 'ic'}, the information criterion-methods discussed in Hallin and Liška (2007) (when \code{fm.restricted = FALSE})
#' and Alessi, Barigozzi and Capasso (2010) (when \code{fm.restricted = TRUE}) are implemented.
#' The information criterion called by \code{ic.op = 5} (as an argument to \code{fnets} or \code{fnets.factor.model}) is recommended by default.
#' @details For further details, see references.
#' @param x input time series each column representing a time series variable; it is coerced into a \link[stats]{ts} object
#' @param fm.restricted whether to estimate the number of restricted or unrestricted factors
#' @param method A string specifying the factor number selection method; possible values are:
#' \itemize{
#'    \item{\code{"ic"}}{ information criteria-based methods of Alessi, Barigozzi & Capasso (2010) when \code{fm.restricted = TRUE} or Hallin and Liška (2007) when \code{fm.restricted = FALSE}}
#'    \item{\code{"er"}}{ eigenvalue ratio of Ahn and Horenstein (2013) when \code{fm.restricted = TRUE} or Avarucci et al. (2022) when \code{fm.restricted = FALSE}}
#' }
#' @param q.max maximum number of factors; if \code{q.max = NULL}, a default value is selected as \code{min(50, floor(sqrt(min(dim(x)[2] - 1, dim(x)[1]))))}
#' @param center whether to de-mean the input \code{x}
#' @return S3 object of class \code{factor.number}.
#' If \code{method = "ic"}, a vector containing minimisers of the six information criteria, otherwise, the maximiser of the eigenvalue ratio
#' @seealso \link[fnets]{plot.factor.number}, \link[fnets]{print.factor.number}
#' @example R/examples/factor_number_ex.R
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @references Alessi, L., Barigozzi, M., and Capasso, M. (2010) Improved penalization for determining the number of factors in approximate factor models. Statistics & Probability Letters, 80(23-24):1806–1813.
#' @references Avarucci, M., Cavicchioli, M., Forni, M., & Zaffaroni, P. (2022) The main business cycle shock(s): Frequency-band estimation of the number of dynamic factors.
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022) fnets: An R Package for Network Estimation and Forecasting via Factor-Adjusted VAR Modelling. arXiv preprint arXiv:2301.11675.
#' @importFrom stats var as.ts
#' @export
factor.number <-
  function(x,
           fm.restricted = FALSE,
           method = c("ic","er"),
           q.max = NULL,
           center = TRUE) {
    x <- t(as.ts(x))
    covx <- NULL
    p <- dim(x)[1]

    if(!is.null(q.max)) q.max <- posint(q.max)
    method <- match.arg(method, c("ic", "er"))
    args <- as.list(environment())
    if(method == "er"){
      ifelse(center, mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
      xx <- x - mean.x
      if(!fm.restricted){
        pca <- dyn.pca(xx, q.method = "er")
        q <- pca$q
      } else {
        pca <- static.pca(xx, q.method = "er")
        q <- pca$q
      }
      out <- q
    }

    if(method == "ic"){
      if(!fm.restricted) {
       temp <- hl.factor.number(
          x = x,
          q.max = q.max,
          mm = NULL,
          center = center
        )
      } else {
        temp <- abc.factor.number(
          x = x,
          covx = covx,
          q.max = q.max,
          center = center
        )
      }
      out <- temp$q.hat
    }
    attr(out, "factor") <- ifelse(fm.restricted, "restricted", "unrestricted")
    attr(out, "class") <- "factor.number"
    attr(out, "args") <- args
    if(method == "er") {
      attr(out, "data") <- pca$q.method.out
    } else {
      data <- attr(temp, "data")
      attr(out, "data") <- data
    }
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
#' @param mm a positive integer specifying the kernel bandwidth for dynamic PCA; by default, it is set to \code{floor(4 *(dim(x)[2]/log(dim(x)[2]))^(1/3)))}
#' @param w vector of length \code{2 * mm + 1} containing symmetric weights; if \code{w = NULL}, default weights are generated using the Bartlett kernel and \code{mm}
#' @param center whether to de-mean the input \code{x} row-wise
#' @return a list containing
#' \item{q.hat}{ a vector containing minimisers of the six information criteria}
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @importFrom graphics par abline box axis legend
#' @importFrom stats var
#' @keywords internal
hl.factor.number <-
  function(x,
           q.max = NULL,
           mm = NULL,
           center = TRUE) {
    p <- dim(x)[1]
    n <- dim(x)[2]
    if(is.null(q.max))
      q.max <- min(50, floor(sqrt(min(n - 1, p))))

    ifelse(center, mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
    xx <- x - mean.x

    if(is.null(mm))
      mm <-  floor(4 * (n / log(n))^(1 / 3))
    w <- Bartlett.weights(((-mm):mm) / mm)

    p.seq <- floor(3 * p / 4 + (1:10) * p / 40)
    n.seq <- n - (9:0) * floor(n / 20)
    const.seq <- seq(.001, 2, by = .01)
    IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 2 * 3))

    for (kk in 1:10) {
      nn <- n.seq[kk]
      pp <- p.seq[kk]
      pen <- c((1 / mm^2 + sqrt(mm / nn) + 1 / pp) * log(min(pp, mm^2, sqrt(nn / mm))),
               1 / sqrt(min(pp, mm^2, sqrt(nn / mm))),
               1 / min(pp, mm^2, sqrt(nn / mm)) * log(min(pp, mm^2, sqrt(nn / mm))))

      Gamma_x <- array(0, dim = c(pp, pp, 2 * mm + 1))
      for (h in 0:(mm - 1)) {
        Gamma_x[, , h + 1] <-
          xx[1:pp, 1:(nn - h)] %*% t(xx[1:pp, 1:(nn - h) + h]) / nn * w[h + mm + 1]
        if(h != 0) Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
      }
      Sigma_x <-
        aperm(apply(Gamma_x, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)

      tmp <- rep(0, q.max + 1)
      for (ii in 1:(mm + 1)) {
        ifelse(kk == length(n.seq),nu <- q.max, nu <- 0)
        sv <- svd(Sigma_x[,, ii], nu = nu, nv = 0)
        dd <- sum(sv$d)
        tmp[1] <- tmp[1] + dd / pp / (2 * mm + 1)
        for (jj in 1:q.max) {
          dd <- dd - sv$d[jj]
          tmp[jj + 1] <- tmp[jj + 1] + dd / pp / (2 * mm + 1)
        }
        for (jj in 1:length(const.seq)) {
          for (ic.op in 1:3) {
            IC[, jj, kk, 3 * 0 + ic.op] <-
              tmp + (0:q.max) * const.seq[jj] * pen[ic.op]
            IC[, jj, kk, 3 * 1 + ic.op] <-
              log(tmp) + (0:q.max) * const.seq[jj] * pen[ic.op]
          }
        }
      }
    }

    q.mat <- apply(IC, c(2, 3, 4), which.min)
    Sc <- apply(q.mat, c(1, 3), var)
    q.hat <- rep(0, 6)
    for (ii in 1:6) {
      ss <- Sc[, ii]
      if(min(ss) > 0) {
        q.hat[ii] <- min(q.mat[max(which(ss == min(ss))), , ii]) - 1
      } else {
        if(sum(ss[-length(const.seq)] != 0 & ss[-1] == 0)) {
          q.hat[ii] <-
            q.mat[which(ss[-length(const.seq)] != 0 &
                          ss[-1] == 0)[1] + 1, 10, ii] - 1
        } else {
          q.hat[ii] <- min(q.mat[max(which(ss == 0)), , ii]) - 1
        }
      }
    }



    out <-
      list(
        q.hat = q.hat
      )
    attr(out, "data") <- list(Sc=Sc, const.seq=const.seq, q.mat=q.mat)
    return(out)
  }

#' @title Factor number estimator of Alessi, Barigozzi and Capasso (2010)
#' @description Estimates the number of factors by minimising an information criterion over sub-samples of the data.
#' Currently the three information criteria proposed in Alessi, Barigozzi and Capasso (2010) (\code{ic.op = 1, 2, 3})
#' and their variations with logarithm taken on the cost (\code{ic.op = 4, 5, 6}) are implemented,
#' with \code{ic.op = 5} recommended as a default choice based on numerical experiments.
#' @details See Alessi, Barigozzi and Capasso (2010) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param covx covariance of \code{x}
#' @param q.max maximum number of factors; if \code{q.max = NULL}, a default value is selected as \code{min(50, floor(sqrt(min(dim(x)[2] - 1, dim(x)[1]))))}
#' @param center whether to de-mean the input \code{x} row-wise
#' @return a list containing
#' \item{q.hat}{ the mimimiser of the chosen information criteria}
#' @references Alessi, L., Barigozzi, M.,  & Capasso, M. (2010) Improved penalization for determining the number of factors in approximate factor models. Statistics & Probability Letters, 80(23-24):1806–1813.
#' @importFrom graphics abline
#' @keywords internal
abc.factor.number <-
  function(x,
           covx = NULL,
           q.max = NULL,
           center = TRUE) {
    p <- dim(x)[1]
    n <- dim(x)[2]
    ifelse(center, mean.x <- apply(x, 1, mean), mean.x <- rep(0, p))
    xx <- x - mean.x

    if(is.null(q.max))
      q.max <- min(50, floor(sqrt(min(n - 1, p))))

    if(is.null(covx))
      covx <- xx %*% t(xx) / n

    p.seq <- floor(3 * p / 4 + (1:10) * p / 40)
    n.seq <- n - (9:0) * floor(n / 20)
    const.seq <- seq(.001, 2, by = .01)
    IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 6))

    for (kk in 1:10) {
      nn <- n.seq[kk]
      pp <- p.seq[kk]
      pen <- c((nn + pp) / (nn * pp) * log(nn * pp / (nn + pp)),
               (nn + pp) / (nn * pp) * log(min(nn, pp)),
               log(min(nn, pp)) / min(nn, pp))

      tmp <- rep(0, q.max + 1)
      ifelse(kk == length(n.seq), nu <- q.max, nu <- 0)
      sv <- svd(covx[1:pp, 1:pp], nu = nu, nv = 0)
      dd <- sum(sv$d)
      tmp[1] <- tmp[1] + dd / pp
      for (jj in 1:q.max) {
        dd <- dd - sv$d[jj]
        tmp[jj + 1] <- tmp[jj + 1] + dd / pp
      }
      for (jj in 1:length(const.seq)) {
        for (ic.op in 1:3) {
          IC[, jj, kk, ic.op] <-
            tmp + (0:q.max) * const.seq[jj] * pen[ic.op]
          IC[, jj, kk, 3 * 1 + ic.op] <-
            log(tmp) + (0:q.max) * const.seq[jj] * pen[ic.op]
        }
      }
    }

    q.mat <- apply(IC, c(2, 3, 4), which.min)
    Sc <- apply(q.mat, c(1, 3), var)
    q.hat <- rep(0, 6)
    for (ii in 1:6) {
      ss <- Sc[, ii]
      if(min(ss) > 0) {
        q.hat[ii] <- min(q.mat[max(which(ss == min(ss))), , ii]) - 1
      } else {
        if(sum(ss[-length(const.seq)] != 0 & ss[-1] == 0)) {
          q.hat[ii] <-
            q.mat[which(ss[-length(const.seq)] != 0 &
                          ss[-1] == 0)[1] + 1, 10, ii] - 1
        } else {
          q.hat[ii] <- min(q.mat[max(which(ss == 0)), , ii]) - 1
        }
      }
    }


    out <- list(q.hat = q.hat, sv = sv)
    attr(out, "data") <- list(Sc=Sc, const.seq=const.seq, q.mat=q.mat)
    return(out)
  }

#' @title Plot factor number
#' @method plot factor.number
#' @description Plots the eigenvalue ratio or information criteria from a \code{factor.number} object
#' @param x \code{factor.number} object
#' @param ... not used
#' @return NULL, printed to console
#' @seealso \link[fnets]{factor.number}
#' @example R/examples/factor_number_ex.R
#' @importFrom graphics par abline box axis legend
#' @export
plot.factor.number <- function(x, ...){
  data <- attr(x, "data")
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if(attr(x, "args")$method == "er"){
    par(xpd = FALSE)
    plot(data, xlab = "q", ylab = "Eigenvalue Ratio")
    abline(v = x)
  } else {
    par(mfrow = c(2, 3))
    Sc <- data$Sc
    const.seq <- data$const.seq
    q.mat <- data$q.mat
      for (ii in 1:6) {
        plot(
          const.seq,
          q.mat[, 10, ii] - 1,
          type = "b",
          pch = 1,
          col = 2,
          bty = "n",
          axes = FALSE,
          xlab = "constant",
          ylab = "",
          main = paste("IC ", ii)
        )
        box()
        axis(1, at = pretty(range(const.seq)))
        axis(
          2,
          at = pretty(range(q.mat[, 10, ii] - 1)),
          col = 2,
          col.ticks = 2,
          col.axis = 2
        )
        par(new = TRUE)
        plot(
          const.seq,
          Sc[, ii],
          col = 4,
          pch = 2,
          type = "b",
          bty = "n",
          axes = FALSE,
          xlab = "",
          ylab = ""
        )
        axis(
          4,
          at = pretty(range(Sc[, ii])),
          col = 4,
          col.ticks = 4,
          col.axis = 4
        )
        legend(
          "topright",
          legend = c("q", "Sc"),
          col = c(2, 4),
          lty = c(1, 1),
          pch = c(1, 2),
          bty = "n"
        )
      }
  }
}

#' @title Print factor number
#' @method print factor.number
#' @description Prints a summary of a \code{factor.number} object
#' @param x \code{factor.number} object
#' @param ... not used
#' @return NULL, printed to console
#' @seealso \link[fnets]{factor.number}
#' @example R/examples/factor_number_ex.R
#' @export
print.factor.number <- function(x, ...){
  args <- attr(x, "args")
  cat(paste("Factor number selection \n"))
  cat(paste("Factor model: ", attr(x, "factor"), "\n", sep = ""))
  met <- ifelse(args$method == "er", "Eigenvalue ratio", "Information criterion")
  cat(paste("Method: ", met, "\n", sep = ""))
  if(args$method == "er"){
    cat(paste("Number of factors: ", x, "\n"), sep = "")
  }
  else{
    cat(paste("Number of factors:", "\n", sep = ""))
    for(ii in 1:6) cat(paste("IC", ii, ": ", x[ii], "\n", sep = ""))
  }
}
