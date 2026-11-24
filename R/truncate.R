#' Sample Covariance Function Without Centering
#'
#' This function computes the sample covariance matrix with an optional lag.
#'
#' @param data A numeric matrix.
#' @param lag An integer specifying the lag.
#' @return A numeric matrix representing the covariance matrix.
#' @useDynLib fnets
#' @importFrom Rcpp sourceCpp
acf_no_center <- function(data, lag = 0) {
  .Call("_fnets_acf_no_center", data, lag)
}


#' Sample Covariance Function Without Centering, after truncating data.
#'
#' This function truncated the data and computes the sample covariance matrix with an optional lag.
#'
#' @param data A numeric matrix.
#' @param tau A numeric vector.
#' @param lag An integer specifying the lag.
#' @return A numeric matrix representing the covariance matrix.
#' @useDynLib fnets
#' @importFrom Rcpp sourceCpp
truncateAndComputeCovariance_lag <- function(data, tau, lag = 0) {
  .Call("_fnets_truncateAndComputeCovariance_lag", data, tau, lag)
}


#' @title Truncate data, with truncation parameter chosen by cross-validation.
#' @param data Input time series each column representing a time series variable.
#' @param n_tau An integer that determines the number of taus to use in grid used for cross-validation
#' @param lag This is an integer argument that is used when the \code{cv_trunc} function is used for (auto)covariance estimation, of particular lag.
#'  The lag determines which (auto)covariance matrix is used in tuning.
#' @param cv_lag An integer argument, that is used when the \code{cv_trunc} function is used with data modeled as a factor-adjusted VAR.
#' The integer determines up to what lag auto-covariance matrix is used in the cv-measure.
#'  In implementation, this will be set as default to be the VAR order. When \code{cv_lag = 0}, only the (auto)covariance matrix, of lag determined
#'  by the \code{lag} argument, is used.
#' @param standardise boolean; whether to scale up the truncation parameter for each series by the MAD of the corresponding series.
#' @export
cv_trunc = function(data, n_tau = 60, lag = 0, cv_lag = 0, standardise = TRUE){

  data = as.matrix(data)

  if(cv_lag==0){
    tau_scores = cross_val(data, n_tau, lag, standardise = standardise)
  }else{
    tau_scores = cross_val_lag(data, n_tau, lag = cv_lag, standardise = standardise)
  }
  min_tau = as.numeric(names(which.min(tau_scores[[1]])))

  if(standardise){
    mad_data = mad_variables(data)
    tau_s = mad_data * min_tau
  } else {
    tau_s = rep(1, ncol(data)) * min_tau
  }

  for(i in 1:nrow(data)){
    x = data[i,]
    data[i,] = ifelse(abs(x) > abs(tau_s), sign(x)*tau_s, x)
  }

  return(list(data = data, tau = min_tau, tau_standardised = tau_s))

}

#' @title internal function that carries out cross validation to choose tuning parameter for data truncation
#' @keywords internal
cross_val = function(data, n_tau = 60, lag, standardise = TRUE){
  # get tau grid, standardised or not
  tau_l = tau_grid_stand_fun(data = data, n_steps = n_tau, standardise = standardise)

  n = nrow(data)
  p = ncol(data)

  # split data
  half_1 = data[1:(n/2), , drop = FALSE]
  half_2 = data[((n/2) + 1):n, , drop = FALSE]

  half_1_trunc = plyr::alply(tau_l$tau_grid, 1, truncateAndComputeCovariance_lag, data = half_1, lag = lag)
  half_2_trunc = plyr::alply(tau_l$tau_grid, 1, truncateAndComputeCovariance_lag, data = half_2, lag = lag)

  half_1_sample = acf_no_center(data = half_1, lag = lag)
  half_2_sample = acf_no_center(data = half_2, lag = lag)

  half_1_trunc_err = lapply(half_1_trunc, function(x){norm((x - half_2_sample), "M")})
  half_2_trunc_err = lapply(half_2_trunc, function(x){norm((x - half_1_sample), "M")})

  scores_per_tau = 1/2 * (unlist(half_1_trunc_err) + unlist(half_2_trunc_err))
  names(scores_per_tau) = tau_l$tau_values[[1]]

  return(list(scores_per_tau, tau_l$tau_values[[2]]))

}

#' @title internal function that carries out cross validation to choose tuning parameter for data truncation
#' @keywords internal
cross_val_lag = function(data, n_tau, lag, standardise = standardise){

  tau_l = tau_grid_stand_fun(data = data, n_steps = n_tau, standardise = standardise)

  n = nrow(data)
  half_1 = data[1:(n/2),]
  half_2 = data[((n/2) + 1):n,]

  scores_per_tau = vector(mode = "list", length = lag+1)

  for(i in 0:lag){
    half_1_trunc = plyr::alply(tau_l$tau_grid, 1, truncateAndComputeCovariance_lag, data = half_1, lag = i)
    half_2_trunc = plyr::alply(tau_l$tau_grid, 1, truncateAndComputeCovariance_lag, data = half_2, lag = i)

    half_1_sample = acf_no_center(data = half_1, lag = i)
    half_2_sample = acf_no_center(data = half_2, lag = i)

    half_1_trunc_err = lapply(half_1_trunc, function(x){norm((x - half_2_sample), "M")})
    half_2_trunc_err = lapply(half_2_trunc, function(x){norm((x - half_1_sample), "M")})

    scores_per_tau[[i+1]] = 1/2 * (unlist(half_1_trunc_err) + unlist(half_2_trunc_err))
  }

  scores_per_tau = do.call(pmax, scores_per_tau)

  names(scores_per_tau) = tau_l$tau_values[[1]]

  return(list(scores_per_tau, tau_l$tau_values[[2]]))
}





#' @title internal function for \code{cv_trunc}
#' @keywords internal
mad_variables = function(data){
  mad_data = vector(length = ncol(data))
  for(i in 1:ncol(data)){
    mad_data[i] = stats::mad(data[,i])
  }
  if(any(mad_data == 0)){
    for(i in 1:ncol(data)){
      mad_data[i] = stats::sd(data[,i])
    }
  }
  return(mad_data)
}

#' @title internal function for \code{cv_trunc}
#' @keywords internal
tau_grid_fun = function(data, n_steps){
  max_dat = max(abs(data))
  median_dat = stats::median(abs(data))
  tau_grid = seq(median_dat, max_dat, length.out = n_steps)
  tau_grid_quant = sapply(tau_grid, function(x){stats::ecdf(abs(data))(x)})

  return(list(tau_grid, tau_grid_quant))
}

#' @title internal function for \code{cv_trunc}
#' @keywords internal
tau_grid_stand_fun = function(data, n_steps, standardise){
  if(standardise){
    mad_data = mad_variables(data)
    tau_values = tau_grid_fun( t(t(data)/mad_data), n_steps = n_steps)
    mad_data_mat = matrix(mad_data, nrow = length(tau_values[[1]]), ncol = length(mad_data), byrow = TRUE)
    tau_grid = mad_data_mat * tau_values[[1]]
  } else {
    # just a matrix of 1s no scaling
    mad_data_mat = matrix(1, nrow = n_steps, ncol = ncol(data), byrow = TRUE)
    tau_values = tau_grid_fun(data, n_steps = n_steps)
    tau_grid = mad_data_mat * tau_values[[1]]
  }

  return(list(tau_values = tau_values, tau_grid = tau_grid))
}




#' @title Internal function for constructing the Lepski grid
#' @description Constructs a geometric grid of candidate truncation
#' parameters, with endpoints determined by the median and maximum
#' of the absolute observations.
#' @param data A numeric matrix, with each column representing a
#' time series variable.
#' @return A numeric vector of candidate truncation parameters
#' in increasing order.
#' @keywords internal
lepski_grid_fun <- function(data) {

  abs_data <- abs(data)

  tau_min <- stats::median(abs_data)
  tau_max <- max(abs_data)

  if (tau_max == 0) {
    stop("Cannot construct Lepski grid for all-zero data")
  }

  # Ensure a strictly positive lower endpoint
  if (tau_min == 0) {
    tau_min <- stats::median(abs_data[abs_data > 0])
  }

  # Number of candidate truncation parameters
  J <- ceiling(log2(tau_max / tau_min)) + 1L

  # Construct the geometric grid
  tau_grid <- tau_min * 2^(seq_len(J) - 1L)

  return(tau_grid)
}



#' @title Internal function for Lepski-based truncation parameter selection
#' @description Selects the truncation parameter using Lepski's method,
#' by comparing truncated sample (auto)covariance matrices obtained
#' over a geometric grid of candidate truncation parameters.
#' @param data A numeric matrix, with each column representing a
#' time series variable.
#' @param lag A non-negative integer specifying the maximum lag
#' of the autocovariance matrices used in the selection criterion.
#' @param C A positive constant used in the Lepski selection criterion.
#' @return A list containing the selected truncation parameter,
#' its grid index, the truncated data, the candidate grid,
#' the Lepski tolerances, the pairwise differences between
#' (auto)covariance estimators and the selected (auto)covariance matrices.
#' @keywords internal
lepski_tau <- function(data, lag = 1, C = 0.072) {

  data <- as.matrix(data)
  storage.mode(data) <- "double"

  n <- nrow(data)
  p <- ncol(data)

  stopifnot(
    all(is.finite(data)),
    p > 1,
    length(lag) == 1L,
    is.finite(lag),
    lag >= 0,
    lag == floor(lag),
    n > lag,
    length(C) == 1L,
    is.finite(C),
    C > 0
  )

  # STEP 1: Generate candidate truncation parameters

  tau_grid <- lepski_grid_fun(data)

  J <- length(tau_grid)

  # STEP 2: Calculate truncated covariance and autocovariance
  # matrices for each candidate tau and each lag 0,...,lag

  Gamma_hat <- lapply(tau_grid, function(tau) {

    lapply(0:lag, function(h) {

      truncateAndComputeCovariance_lag(
        data = data,
        tau = rep(tau, p),
        lag = h
      )

    })

  })

  # STEP 3: Calculate the Lepski tolerance R(tau)

  R_tau <- tau_grid^2 * log(n)^2 * log(p) / n

  # STEP 4: Calculate the maximum element-wise difference
  # between each pair of estimators, across all lags

  differences <- matrix(0, nrow = J, ncol = J)

  for (j in seq_len(J - 1L)) {

    for (l in (j + 1L):J) {

      differences[j, l] <- max(
        vapply(seq_len(lag + 1L), function(h) {

          max(abs(
            Gamma_hat[[j]][[h]] -
              Gamma_hat[[l]][[h]]
          ))

        }, numeric(1))
      )

    }

  }

  # STEP 5: Select the smallest j satisfying the Lepski
  # criterion for every larger candidate l

  j_hat <- J

  for (j in seq_len(J - 1L)) {

    larger_indices <- (j + 1L):J

    if (all(
      differences[j, larger_indices] <=
      C * R_tau[larger_indices]
    )) {

      j_hat <- j
      break

    }

  }

  tau_hat <- tau_grid[j_hat]

  # STEP 6: Truncate the data using the selected tau

  data_truncated <- pmax(pmin(data, tau_hat), -tau_hat)

  return(list(
    tau = tau_hat,
    j_hat = j_hat,
    data = data_truncated,
    tau_grid = tau_grid,
    R_tau = R_tau,
    differences = differences,
    Gamma_hat = Gamma_hat[[j_hat]]
  ))
}






#' @title Truncate data, with truncation parameter chosen by Lepski's method
#' @description Performs element-wise data truncation, with the truncation
#' parameter selected using Lepski's method based on comparisons of
#' truncated sample (auto)covariance matrices.
#' @param data Input time series, with each column representing a
#' time series variable.
#' @param lag A non-negative integer specifying the maximum lag
#' of the autocovariance matrices used in the Lepski criterion.
#' @param C A positive constant used in the Lepski selection criterion.
#' The default value is 0.072.
#' @param standardise Boolean; whether to scale up the truncation
#' parameter for each series by the MAD of the corresponding series.
#' @details The candidate truncation parameters are constructed on a
#' geometric grid, with successive candidates differing by a factor of two.
#' The lower and upper endpoints are determined by the median and maximum
#' of the absolute observations, respectively.
#'
#' When \code{standardise = TRUE}, the Lepski selection procedure is
#' applied to the MAD-standardised data, and the selected truncation
#' parameter is rescaled for each time series.
#'
#' The smallest candidate truncation parameter for which the corresponding
#' (auto)covariance estimators are sufficiently close to those obtained
#' using all larger candidate values is selected.
#' @return A list containing:
#' \item{data}{The truncated data on the original scale.}
#' \item{tau}{The selected truncation parameter.}
#' \item{tau_standardised}{A vector of series-specific truncation
#' parameters on the original scale.}
#' \item{j_hat}{The index of the selected truncation parameter.}
#' \item{tau_grid}{The grid of candidate truncation parameters.}
#' \item{R_tau}{The Lepski tolerances associated with the candidate
#' truncation parameters, before multiplication by \code{C}.}
#' \item{differences}{A matrix containing the maximum element-wise
#' differences between pairs of (auto)covariance estimators,
#' across the considered lags. Only the upper triangular entries
#' are computed.}
#' \item{Gamma_hat}{A list of the selected covariance and autocovariance
#' estimators, from lag zero to \code{lag}, on the scale used
#' for the Lepski selection procedure.}
#' @references Loh, P. L. (2021). Scale calibration for high-dimensional
#' robust regression. Electronic Journal of Statistics, 15(2), 5933--5994.
#' @seealso \link{cv_trunc}
#' @export
lepski_trunc <- function(data, lag = 1, C = 0.072,
                         standardise = TRUE) {

  data <- as.matrix(data)

  stopifnot(
    is.numeric(data),
    all(is.finite(data)),
    is.logical(standardise),
    length(standardise) == 1L,
    !is.na(standardise)
  )

  # Obtain the scale of each time series

  if (standardise) {

    scales <- mad_variables(data)

  } else {

    scales <- rep(1, ncol(data))

  }

  if (any(!is.finite(scales) | scales <= 0)) {
    stop("Invalid scale values for Lepski standardisation")
  }

  # Standardise the data

  data_scaled <- sweep(data, 2, scales, "/")

  # Apply Lepski's method

  fit <- lepski_tau(
    data = data_scaled,
    lag = lag,
    C = C
  )

  # Obtain series-specific truncation parameters

  tau_s <- fit$tau * scales

  # Truncate the original data using the selected thresholds

  data_truncated <- data

  for (i in seq_len(ncol(data))) {

    data_truncated[, i] <- pmax(
      pmin(data[, i], tau_s[i]),
      -tau_s[i]
    )

  }

  return(list(
    data = data_truncated,
    tau = fit$tau,
    tau_standardised = tau_s,
    j_hat = fit$j_hat,
    tau_grid = fit$tau_grid,
    R_tau = fit$R_tau,
    differences = fit$differences,
    Gamma_hat = fit$Gamma_hat
  ))
}
