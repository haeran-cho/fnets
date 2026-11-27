#' @title Bartlett weights
#' @keywords internal
Bartlett.weights <- function(x)
  1 - abs(x)

#' @keywords internal
make.symmetric <- function(DD, symmetric) {
  symmetric <- match.arg(symmetric, c("min", "max", "avg", "none"))
  if(symmetric != "none") {
    if(symmetric == "min") {
      DD[abs(DD) > t(abs(DD))] <- 0
      DD <- DD + t(DD)
      diag(DD) <- diag(DD) / 2
    }
    if(symmetric == "max") {
      DD[abs(DD) < t(abs(DD))] <- 0
      DD <- DD + t(DD)
      diag(DD) <- diag(DD) / 2
    }
    if(symmetric == "avg")
      DD <- (DD + t(DD)) / 2
  }
  DD
}


#' @keywords internal
var.to.vma <- function(A, trunc.lags) {
  d <- dim(A)[1]
  s <- dim(A)[2]
  l <- s / d
  B <- array(0, dim = c(d, d, trunc.lags + 1))
  B[, , 1] <- diag(1, d)
  for (ii in 1:trunc.lags) {
    for (jj in 1:min(ii, l))
      B[, , ii + 1] <-
        B[, , ii + 1] + B[, , ii - jj + 1] %*% A[, (jj - 1) * d + 1:d]
  }
  B
}

#' @keywords internal
check.list.arg <- function(arg) {
  arg.name <- deparse(substitute(arg))
  if(arg.name == "var.args")
    default <-
      list(
        n.iter = 100,
        tol = 0,
        n.cores = min(parallel::detectCores() - 1, 3)
      )
  if (arg.name == "tuning.args")
    default <- list(
      tuning = c("cv", "bic"),
      cv.method = c("standard", "rolling"),
      n.folds = 10,
      n.blocks = 5,
      penalty = NULL,
      path.length = 10
    )
  if(arg.name == "common.args")
    default <-
      list(
        var.order = NULL,
        max.var.order = NULL,
        trunc.lags = 20,
        n.perm = 10
      )
  arg <- check.list.valid(arg)
  out <- check.list.match(arg, default)
  return(out)
}

#' @keywords internal
check.list.match <- function(arg, default) {
  for (ii in names(default)) {
    if(is.null(arg[[ii]]))
      arg[[ii]] <- default[[ii]]
  }
  return(arg)
}

#' @keywords internal
check.list.valid <- function(arg) {
  arg.name <- deparse(substitute(arg))
  if(arg.name == "var.args"){
    arg$n.iter <- posint(arg$n.iter)
    arg$tol <- max(0, arg$tol)
    arg$n.cores <- posint(arg$n.cores)
  }
  if(arg.name == "tuning.args"){
    arg$n.folds <- posint(arg$n.folds)
    arg$path.length <- posint(arg$path.length)
  }
  if(arg.name == "common.args"){
    arg$var.order <- posint(arg$var.order)
    arg$max.var.order <- posint(arg$max.var.order)
    arg$trunc.lags <- posint(arg$trunc.lags)
    arg$n.perm <- posint(arg$n.perm)
  }
  return(arg)
}

#' @keywords internal
posint <- function(x, min.x = 1){
  x.name <- deparse(substitute(x))
  if(x < min.x | x != round(x)){
    xx <- max(min.x, as.integer(x))
    if(xx == round(x))
      warning(cat(x.name,"=",x, "coerced to ", xx,"\n"))
    else
      stop(paste(x.name,"=",x, "cannot be coerced to correct input format. Must be integer greater than", min.x,"\n"))
    x <- xx
  }
  return(x)
}




#' @title Rolling cross-validation for glmnet VAR estimation
#'
#' @description
#' Selects the glmnet regularisation parameter using an expanding-window
#' rolling cross-validation scheme that preserves temporal ordering.
#'
#' The observations are first divided into consecutive blocks according to
#' their original time indices. The first block is used as the initial
#' training sample. For each subsequent block, all earlier blocks are used
#' for training and the current block is used for validation.
#'
#' Importantly, all p response coordinates corresponding to the same time
#' point are assigned to the same training or validation set.
#'
#' @param yy Vectorised response used by glmnet. For a VAR(d), this is ordered
#'   by response coordinate, with time running backwards within each coordinate:
#'   \eqn{x_{1,n}, \ldots, x_{1,d+1}, x_{2,n}, \ldots, x_{2,d+1},
#'   \ldots, x_{p,n}, \ldots, x_{p,d+1}}.
#'
#' @param xxs Sparse design matrix corresponding to yy.
#'
#' @param n Number of time observations in the original time series.
#'
#' @param p Number of component series.
#'
#' @param var.order VAR order d.
#'
#' @param n.blocks Number of consecutive time blocks used for rolling
#'   cross-validation. Since the first block is used only for initial training,
#'   n.blocks gives n.blocks - 1 validation steps.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{lambda.min}: lambda minimising the rolling CV error;
#'   \item \code{lambda}: common lambda path considered;
#'   \item \code{cvm}: mean validation error for each lambda;
#'   \item \code{glmnet.fit}: glmnet fit to the full dataset;
#'   \item \code{blocks}: time indices belonging to each CV block.
#' }
#'
#' @keywords internal
rolling.cv.glmnet <- function(yy,
                              xxs,
                              n,
                              p,
                              var.order,
                              n.blocks = 5) {

  n.time <- n - var.order

  if(n.blocks < 2)
    stop("n.blocks must be at least 2.")

  if(n.blocks > n.time)
    stop("n.blocks cannot exceed the number of usable time points.")


  ## Full-data fit gives the common lambda path
  full.fit <- glmnet::glmnet(
    y = yy,
    x = xxs,
    intercept = FALSE,
    standardize = FALSE
  )

  lambda.path <- full.fit$lambda


  ## Consecutive blocks defined according to time
  times <- (var.order + 1):n

  block.id <- cut(
    seq_along(times),
    breaks = n.blocks,
    labels = FALSE
  )

  blocks <- split(times, block.id)


  ## Map time indices to rows in the vectorised glmnet regression
  get.rows <- function(tt) {

    unlist(
      lapply(seq_len(p), function(j) {
        (j - 1) * n.time + (n - tt + 1)
      }),
      use.names = FALSE
    )
  }


  cv.sse <- rep(0, length(lambda.path))
  n.validation <- 0


  ## Expanding-window CV
  for(k in 2:n.blocks) {

    train.times <- unlist(blocks[1:(k - 1)])
    test.times <- blocks[[k]]

    train.rows <- get.rows(train.times)
    test.rows <- get.rows(test.times)


    fit.k <- glmnet::glmnet(
      y = yy[train.rows],
      x = xxs[train.rows, , drop = FALSE],
      lambda = lambda.path,
      intercept = FALSE,
      standardize = FALSE
    )


    pred.k <- stats::predict(
      fit.k,
      newx = xxs[test.rows, , drop = FALSE],
      s = lambda.path,
      type = "response"
    )

    pred.k <- as.matrix(pred.k)


    err.k <- sweep(
      pred.k,
      MARGIN = 1,
      STATS = yy[test.rows],
      FUN = "-"
    )

    cv.sse <- cv.sse + colSums(err.k^2)

    n.validation <- n.validation + length(test.rows)
  }


  cv.error <- cv.sse / n.validation

  lambda.min <- lambda.path[which.min(cv.error)]


  list(
    lambda.min = lambda.min,
    lambda = lambda.path,
    cvm = cv.error,
    glmnet.fit = full.fit,
    blocks = blocks
  )
}









#' @title Internal function for checking robust arguments
#' @description Checks the list of arguments used for selecting
#' the truncation parameter and fills in missing entries
#' with their default values.
#' @param robust.args A list specifying the arguments used
#' for truncation parameter selection.
#' @return A list containing the validated arguments.
#' @keywords internal
check.robust.args <- function(robust.args) {

  if (!is.list(robust.args)) {
    stop("robust.args must be a list")
  }

  # Check argument names

  if (length(robust.args) > 0L) {

    if (is.null(names(robust.args)) ||
        anyNA(names(robust.args)) ||
        any(!names(robust.args) %in% c("tuning", "LepskiC")) ||
        anyDuplicated(names(robust.args))) {

      stop("Invalid argument names in robust.args")

    }

  }

  # Fill in missing arguments

  if (is.null(robust.args$tuning)) {
    robust.args$tuning <- "cv"
  }

  if (is.null(robust.args$LepskiC)) {
    robust.args$LepskiC <- 0.072
  }

  # Check the tuning method

  robust.args$tuning <- match.arg(
    robust.args$tuning,
    c("cv", "lepski")
  )

  # Check the Lepski constant

  if (!is.numeric(robust.args$LepskiC) ||
      length(robust.args$LepskiC) != 1L ||
      !is.finite(robust.args$LepskiC) ||
      robust.args$LepskiC <= 0) {

    stop("LepskiC must be a positive finite number")

  }

  return(robust.args)
}



