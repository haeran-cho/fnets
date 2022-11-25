#' @title Factor number selection methods
#' @description Methods to estimate the number of factor.
#' Either maximises the ratio of successive eigenvalues, or minimises an information criterion over sub-samples of the data.
#' For \code{restricted = FALSE}, the three information criteria proposed in Hallin and Liška (2007) (\code{pen.op = 1, 2} or \code{3})
#' and their variations with logarithm taken on the cost (\code{pen.op = 4, 5} or \code{6}) are implemented,
#' with \code{pen.op = 5} recommended as a default choice based on numerical experiments.
#' For \code{restricted = TRUE}, the three information criteria in Owens, Cho, and Barigozzi (2022) are implemented,
#' with \code{pen.op = 2} recommended by default.
#' @details For further details, see references.
#' @param x input time series matrix, with each row representing a variable
#' @param fm.restricted whether to estimate the number of restricted or unrestricted factors
#' @param method A string specifying the factor number selection method; possible values are:
#' \itemize{
#'    \item{\code{"ic"}}{ information criteria-based methods of Alessi, Barigozzi & Capasso (2010) when \code{fm.restricted = TRUE} or Hallin and Liška (2007) when \code{fm.restricted = FALSE} modifying Bai and Ng (2002)}
#'    \item{\code{"er"}}{ eigenvalue ratio of Ahn and Horenstein (2013)}
#' };
#' @param q.max maximum number of factors; if \code{q.max = NULL}, a default value is selected as \code{min(50, floor(sqrt(min(dim(x)[2] - 1, dim(x)[1]))))}
#' @param mm bandwidth; defaults to \code{floor(4 * (dim(x)[2]/log(dim(x)[2]))^(1/3)))}
#' @param w vector of length \code{2 * mm + 1} containing symmetric weights; if \code{w = NULL}, default weights are generated using the Bartlett kernel and \code{mm}
#' @param covx covariance matrix of x
#' @param do.plot whether to plot the information criteria values
#' @param center whether to de-mean the input \code{x} row-wise
#' @return if \code{method = "ic"}, a list containing
#' \item{q.hat}{ a vector containing minimisers of the six information criteria}
#' \item{sv}{ a list containing the singular value decomposition of \code{Sigma_x}}
#' and if \code{restricted = FALSE}
#' \item{Gamma_x}{ an array containing the estimates of the autocovariance matrices of \code{x} at \code{2 * mm + 1} lags}
#' \item{Sigma_x}{ an array containing the estimates of the spectral density matrices of \code{x} at \code{2 * mm + 1} Fourier frequencies}
#' otherwise
#' \item{q.hat}{ the maximiser of the eigenvalue ratio}
#' \item{pca}{ dynamic or static pca output}
#'
#' @example R/examples/hl_ex.R
#' @example R/examples/abc_ex.R
#' @references Ahn, S. C. & Horenstein, A. R. (2013) Eigenvalue ratio test for the number of factors. Econometrica, 81(3), 1203--1227.
#' @references Alessi, L., Barigozzi, M., and Capasso, M. (2010) Improved penalization for determining the number of factors in approximate factor models. Statistics & Probability Letters, 80(23-24):1806–1813.
#' @references Bai, J. & Ng, S. (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221.
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022)
#' @importFrom graphics par abline box axis legend
#' @importFrom stats var
#' @export
factor.number <-
  function(x,
           fm.restricted = FALSE,
           method = c("ic","er"),
           q.max = NULL,
           mm = NULL,
           w = NULL,
           covx = NULL,
           do.plot = FALSE,
           center = TRUE) {
    p <- dim(x)[1]
    method <- match.arg(method, c("ic", "er"))
    if(method == "er"){
      if (center)
        mean.x <- apply(x, 1, mean)
      else
        mean.x <- rep(0, p)
      xx <- x - mean.x
      if (!fm.restricted){
        pca <- dyn.pca(xx, q.method = "er")
        q <- pca$q
      } else {
        pca <- static.pca(xx, q.method = "er")
        q <- pca$q
      }
      if(do.plot){
        plot(pca$q.method.out, xlab = "q", ylab = "Eigenvalue Ratio")
        abline(v = q)
      }
      out <- list(q.hat = q, pca = pca)
    }

    if(method == "ic"){
    if (!fm.restricted) {
       out <- hl.factor.number(
          x = x,
          q.max = q.max,
          mm = mm,
          w = w,
          do.plot = do.plot,
          center = center
        )

    } else {
      out <-
        abc.factor.number(
          x = x,
          covx = covx,
          q.max = q.max,
          do.plot = do.plot,
          center = center
        )
    }
    }
    return(out)
  }






#' @title Factor number estimator of Hallin and Liška (2007)
#' @description Estimates the number of factors by minimising an information criterion over sub-samples of the data.
#' Currently the three information criteria proposed in Hallin and Liška (2007) (\code{pen.op = 1, 2} or \code{3})
#' and their variations with logarithm taken on the cost (\code{pen.op = 4, 5} or \code{6}) are implemented,
#' with \code{pen.op = 5} recommended as a default choice based on numerical experiments.
#' @details See Hallin and Liška (2007) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param q.max maximum number of factors; if \code{q.max = NULL}, a default value is selected as \code{min(50, floor(sqrt(min(dim(x)[2] - 1, dim(x)[1]))))}
#' @param mm bandwidth; defaults to \code{floor(4 * (dim(x)[2]/log(dim(x)[2]))^(1/3)))}
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
#' @keywords internal
hl.factor.number <-
  function(x,
           q.max = NULL,
           mm = NULL,
           w = NULL,
           do.plot = FALSE,
           center = TRUE) {
    p <- dim(x)[1]
    n <- dim(x)[2]
    if (is.null(q.max))
      q.max <- min(50, floor(sqrt(min(n - 1, p))))

    if (center)
      mean.x <- apply(x, 1, mean)
    else
      mean.x <- rep(0, p)
    xx <- x - mean.x

    if (is.null(mm))
      mm <-  floor(4 * (n / log(n)) ^ (1 / 3))
    if (is.null(w))
      w <- Bartlett.weights(((-mm):mm) / mm)

    p.seq <- floor(3 * p / 4 + (1:10) * p / 40)
    n.seq <- n - (9:0) * floor(n / 20)
    const.seq <- seq(.001, 2, by = .01)
    IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 2 * 3))

    Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * mm + 1))

    for (kk in 1:10) {
      nn <- n.seq[kk]
      pp <- p.seq[kk]
      pen <- c((1 / mm ^ 2 + sqrt(mm / nn) + 1 / pp) * log(min(pp, mm ^ 2, sqrt(nn / mm))),
               1 / sqrt(min(pp, mm ^ 2, sqrt(nn / mm))),
               1 / min(pp, mm ^ 2, sqrt(nn / mm)) * log(min(pp, mm ^ 2, sqrt(nn / mm))))

      for (h in 0:(mm - 1)) {
        Gamma_x[, , h + 1] <-
          xx[, 1:(nn - h)] %*% t(xx[, 1:(nn - h) + h]) / nn
        Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + mm + 1]
        if (h != 0) {
          Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
          Gamma_xw[, , 2 * mm + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
        }
      }
      Sigma_x <-
        aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)
      sv <- list(1:(mm + 1))

      tmp <- rep(0, q.max + 1)
      for (ii in 1:(mm + 1)) {
        if (kk == length(n.seq))
          nu <- q.max
        else
          nu <- 0
        sv[[ii]] <- svd(Sigma_x[1:pp, 1:pp, ii], nu = nu, nv = 0)
        dd <- sum(sv[[ii]]$d)
        tmp[1] <- tmp[1] + dd / pp / (2 * mm + 1)
        for (jj in 1:q.max) {
          dd <- dd - sv[[ii]]$d[jj]
          tmp[jj + 1] <- tmp[jj + 1] + dd / pp / (2 * mm + 1)
        }
        for (jj in 1:length(const.seq)) {
          for (pen.op in 1:3) {
            IC[, jj, kk, 3 * 0 + pen.op] <-
              tmp + (0:q.max) * const.seq[jj] * pen[pen.op]
            IC[, jj, kk, 3 * 1 + pen.op] <-
              log(tmp) + (0:q.max) * const.seq[jj] * pen[pen.op]
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
          q.hat[ii] <-
            q.mat[which(ss[-length(const.seq)] != 0 &
                          ss[-1] == 0)[1] + 1, 10, ii] - 1
        } else {
          q.hat[ii] <- min(q.mat[max(which(ss == 0)), , ii]) - 1
        }
      }
    }

    if (do.plot) {
      par(mfrow = c(2, 3))
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

    ls <-
      list(
        q.hat = q.hat,
        Gamma_x = Gamma_x,
        Sigma_x = Sigma_x,
        sv = sv
      )
    return(ls)
  }









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
#' @example R/examples/abc_ex.R
#' @references preprint
#' @references Alessi, L., Barigozzi, M.,  & Capasso, M. (2010) Improved penalization for determining the number of factors in approximate factor models. Statistics & Probability Letters, 80(23-24):1806–1813.
#' @references Bai, J. & Ng, S. (2002) Determining the number of factors in approximate factor models. Econometrica. 70: 191-221.
#' @references Owens, D., Cho, H. & Barigozzi, M. (2022)
#' @importFrom graphics abline
#' @keywords internal
abc.factor.number <-
  function(x,
           covx = NULL,
           q.max = NULL,
           do.plot = FALSE,
           center = TRUE) {
    p <- dim(x)[1]
    n <- dim(x)[2]
    cnt <- min(n, p)
    if (center)
      mean.x <- apply(x, 1, mean)
    else
      mean.x <- rep(0, p)
    xx <- x - mean.x

    if (is.null(q.max))
      q.max <- min(50, floor(sqrt(min(n - 1, p))))

    if (is.null(covx))
      covx <- xx %*% t(xx) / n

    p.seq <- floor(3 * p / 4 + (1:10) * p / 40)
    n.seq <- n - (9:0) * floor(n / 20)
    const.seq <- seq(.001, 2, by = .01)
    IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 3))


    for (kk in 1:10) {
      nn <- n.seq[kk]
      pp <- p.seq[kk]
      pen <- c((n + p) / (n * p) * log(n * p / (n + p)),
               (n + p) / (n * p) * log(cnt),
               log(cnt) / cnt)

      mm <- floor(n / 4) - 1
      #sv <- list(1:(mm + 1))
      tmp <- rep(0, q.max + 1)
      if (kk == length(n.seq))
        nu <- q.max
      else
        nu <- 0
      sv <- svd(covx[1:pp, 1:pp], nu = nu, nv = 0)
      dd <- sum(sv$d)
      tmp[1] <- tmp[1] + dd / pp #/ (2 * mm + 1)
      for (jj in 1:q.max) {
        dd <- dd - sv$d[jj]
        tmp[jj + 1] <- tmp[jj + 1] + dd / pp #/ (2 * mm + 1)
      }
      for (jj in 1:length(const.seq)) {
        for (pen.op in 1:3) {
          IC[, jj, kk, pen.op] <-
            tmp + (0:q.max) * const.seq[jj] * pen[pen.op]
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
          q.hat[ii] <-
            q.mat[which(ss[-length(const.seq)] != 0 &
                          ss[-1] == 0)[1] + 1, 10, ii] - 1
        } else {
          q.hat[ii] <- min(q.mat[max(which(ss == 0)), , ii]) - 1
        }
      }
    }

    if (do.plot) {
      par(mfrow = c(1, 3))
      for (ii in 1:3) {
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
    return(list(q.hat = q.hat, sv = sv))
  }
