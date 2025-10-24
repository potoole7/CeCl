#' @title Generate Random Variates from the Generalized Pareto Distribution
#' @description
#' Generates random variates from the Generalized Pareto Distribution (GPD)
#' using the inverse transform sampling method.
#' @param n Integer. The number of random variates to generate.
#' @param u Numeric. The threshold (location parameter) of the GPD.
#' @param sigma Numeric. The scale parameter of the GPD (must be positive).
#' @param xi Numeric. The shape parameter of the GPD.
#' @return A numeric vector of length `n` containing random variates from the
#' GPD.
#' @rdname rgpd
#' @keywords internal
rgpd <- \(n, u, sigma, xi) {
  p <- stats::runif(n)
  if (abs(xi) > 1e-6) {
    u + sigma / xi * ((1 - p)^(-xi) - 1)
  } else {
    u - sigma * log(1 - p)
  }
}
