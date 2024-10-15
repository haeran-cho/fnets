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





