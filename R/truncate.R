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


#' @title internal function that carries out data truncation
#' @param n_tau An integer that determines the number of taus to use in grid used for cross-validation
#' @param lag This is an integer argument that is used when the \code{cv_trunc} function is used for (auto)covariance estimation, of particular lag.
#'  The lag determines which (auto)covariance matrix is used in tuning.
#' @param cv_lag an integer vector denoting any group structure of the vertices
#' @param standardise an integer vector denoting any group structure of the vertices
#' @keywords internal
cv_trunc = function(data, n_tau = 60, lag = 0, cv_lag= F, standardise = T){

  if(!cv_lag){
    tau_scores = cross_val(data, n_tau, lag, standardise = standardise)
  }else{
    # default set to getting cv error using up to lag 1, if cv_lag = T
    tau_scores = cross_val_lag(data, n_tau, lag = 1, standardise = standardise)
  }
  # tau_scores = cross_val(data, n_tau, lag, trim, trim_d, trim_all, max = max)
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

  return(data)

}

#' @title internal function that carries out cross validation to choose tuning parameter for data truncation
#' @description
#' @keywords internal
cross_val = function(data, n_tau = 60, lag, standardise = T){

  # get tau grid, standardised or not
  tau_l = tau_grid_stand_fun(data = data, n_steps = n_tau, standardise = standardise)
  # split data
  n = nrow(data)
  half_1 = data[1:(n/2),]
  half_2 = data[((n/2) + 1):n,]

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
#' @description
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


