#' @title Calculate within cluster sum of distances
#' @description Calculate within cluster sum of distances
#' @param k Number of clusters
#' @param distance_matrix Distance matrix
#' @param fun Clustering function
#' @return Total within-cluster sum of distances
#' @keywords internal
# compute the total within-cluster sum of distances
# TODO: Create methods for the below functions to differ for PAM vs k-means
within_cluster_sum <- function(k, distance_matrix, fun = cluster::pam, ...) {
  clust_res <- fun(distance_matrix, k, ...)
  if (inherits(clust_res, "kmeans")) {
    return(clust_res$tot.withinss)
  } else if (inherits(clust_res, "pam")) {
    return(clust_res$objective[1])
  } else {
    stop("Clustering class not currently supported")
  }
}


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

#' @title Generalized Pareto Distribution density function
#' @description PDF for Generalized Pareto Distribution (GPD).
#' @param x Vector of quantiles.
#' @param u Threshold parameter.
#' @param sigma Scale parameter.
#' @param xi Shape parameter.
#' @return Vector of density values.
#' @rdname dgpd
#' @export
dgpd <- \(x, u, sigma, xi) {
  if (abs(xi) > 1e-6) {
    dens <- (1 / sigma) * (1 + xi * (x - u) / sigma)^(-1 / xi - 1)
    dens[x < u] <- 0
    dens[(1 + xi * (x - u) / sigma) <= 0] <- 0
    dens
  } else {
    dens <- (1 / sigma) * exp(-(x - u) / sigma)
    dens[x < u] <- 0
    dens
  }
}

#' @title Generalized Pareto Distribution cumulative distribution function
#' @description CDF for Generalized Pareto Distribution (GPD).
#' @param q Vector of quantiles.
#' @param u Threshold parameter.
#' @param sigma Scale parameter.
#' @param xi Shape parameter.
#' @return Vector of CDF values.
#' @rdname pgpd
#' @export
pgpd <- \(q, u, sigma, xi) {
  if (abs(xi) > 1e-6) {
    cdf <- 1 - (1 + xi * (q - u) / sigma)^(-1 / xi)
    cdf[q < u] <- 0
    cdf[(1 + xi * (q - u) / sigma) <= 0] <- 1
    cdf
  } else {
    cdf <- 1 - exp(-(q - u) / sigma)
    cdf[q < u] <- 0
    cdf
  }
}

#' @title Generalized Pareto Distribution quantile function
#' @description Quantile function for Generalized Pareto Distribution (GPD).
#' @param p Vector of probabilities.
#' @param u Threshold parameter.
#' @param sigma Scale parameter.
#' @param xi Shape parameter.
#' @return Vector of quantiles.
#' @rdname qgpd
#' @export
qgpd <- \(p, u, sigma, xi) {
  if (abs(xi) > 1e-6) {
    u + (sigma / xi) * ((1 - p)^(-xi) - 1)
  } else {
    u - sigma * log(1 - p)
  }
}

#' @title Generalized Pareto Distribution random generation
#' @description Generate random samples from Generalized Pareto Distribution
#' (GPD).
#' @param n Number of samples to generate.
#' @param u Threshold parameter.
#' @param sigma Scale parameter.
#' @param xi Shape parameter.
#' @return Vector of random samples from GPD.
#' @rdname rgpd
#' @export
rgpd <- \(n, u, sigma, xi) {
  p <- stats::runif(n)
  # use qgpd to get samples
  qgpd(p, u, sigma, xi)
}
