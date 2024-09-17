#include <Rcpp.h>
using namespace Rcpp;

// Define sample covariance function without centering, with lag
// [[Rcpp::export]]
NumericMatrix acf_no_center(NumericMatrix data, int lag = 0) {
  int n = data.nrow(); // Number of rows in the data matrix
  int p = data.ncol(); // Number of columns in the data matrix
  NumericMatrix cov_mat(p, p); // Covariance matrix

  // Loop through each pair of columns
  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      double sum = 0.0; // Sum of product of elements
      // Compute sum of product of elements for each pair of rows with lag
      for (int k = 0; k < n - lag; k++) {
        sum += data(k+lag, i) * data(k, j);
      }
      // Compute sample covariance without centering with lag
      cov_mat(i, j) = sum / (n - lag);
    }
  }

  return cov_mat;
}



// [[Rcpp::export]]
NumericMatrix truncateAndComputeCovariance_lag(NumericMatrix data, NumericVector tau, int lag = 0) {
  int n = data.nrow(); // Number of rows in the data matrix
  int p = data.ncol(); // Number of columns in the data matrix
  NumericMatrix data_copy = Rcpp::clone(data); // Create a copy of the data matrix
  NumericMatrix covariance(p, p); // Matrix to store the covariance

  // Loop through each row of the data matrix
  for (int i = 0; i < n; i++) {
    // Loop through each column of the data matrix
    for (int j = 0; j < p; j++) {
      // Truncate the data element using the corresponding tau value
      if (std::abs(data_copy(i, j)) > tau[j]) {
        data_copy(i, j) = tau[j] * std::copysign(1.0, data_copy(i, j)); // Truncate to tau with sign
      }
    }
  }

  // Compute the covariance without estimating the sample mean
  if (lag > 0) {
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < p; j++) {
        double sum = 0.0;
        for (int k = lag; k < n; k++) {
          sum += data_copy(k, i) * data_copy(k-lag, j);
        }
        covariance(i, j) = sum / (n - lag);
      }
    }
  } else {
    for (int i = 0; i < p; i++) {
      for (int j = i; j < p; j++) { // Only need to compute the upper triangular part
        double sum = 0.0;
        for (int k = 0; k < n; k++) {
          sum += data_copy(k, i) * data_copy(k, j);
        }
        covariance(i, j) = sum / (n);
        covariance(j, i) = covariance(i, j); // Covariance matrix is symmetric
      }
    }
  }

  return covariance;
}

