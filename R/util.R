#' @title Convert to matrix
#' @description Convert input to matrix if it is a vector.
#' @param F_hat Input data.
#' @return Matrix of input data.
#' @rdname to_matrix
#' @keywords internal
to_matrix <- \(F_hat) {
  ret <- F_hat
  if (!is.matrix(F_hat) && is.vector(F_hat)) {
    ret <- as.matrix(F_hat)
  }
  ret
}

#' @title Laplace transformation
#' @description Transform data to Laplace margins.
#' @param F_hat Matrix of CDF values.
#' @param tol Tolerance to avoid issues at 0 and 1.
#' @return Matrix of Laplace-transformed values.
#' @rdname plaplace
#' @export
dlaplace <- \(F_hat, tol = .Machine$double.eps) {
  apply(to_matrix(F_hat), 2, \(x) {
    y <- pmin(pmax(x, tol), 1 - tol)
    ifelse(y < 0.5, log(2 * y), -log(2 * (1 - y)))
  })
}

#' @title Inverse Laplace transformation
#' @description Transform data from Laplace margins back to original scale.
#' @param F_hat Matrix of Laplace-transformed values.
#' @return Matrix of CDF values.
#' @rdname plaplace
#' @export
plaplace <- \(F_hat) {
  apply(to_matrix(F_hat), 2, \(x) {
    ifelse(x < 0, exp(x) / 2, 1 - exp(-x) / 2)
  })
}

#' @title Laplace quantile function
#' @description Quantile function for Laplace distribution.
#' @param p Vector of probabilities.
#' @return Vector of quantiles.
#' @rdname qlaplace
#' @export
qlaplace <- \(p) {
  ifelse(p < 0.5, log(2 * p), -log(2 * (1 - p)))
}

#' @title Laplace random generation
#' @description Generate random samples from Laplace distribution.
#' @param n Number of samples to generate.
#' @return Vector of random samples.
#' @rdname rlaplace
#' @export
rlaplace <- \(n) {
  u <- stats::runif(n)
  qlaplace(u)
}
