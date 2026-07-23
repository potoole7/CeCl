#' @title Generic clustering function
#' @description Generic clustering function for different object classes.
#' @param x Object to cluster.
#' @param ... Additional arguments passed to methods.
#' @return Clustering results.
#' @rdname cecl_clust
#' @export
cecl_clust <- \(x, ...) {
  UseMethod("cecl_clust")
}

# TODO How to only cluster for one variable? Replace vars with var argument?
# TODO Create another method which just clusters, no distance calculation
# TODO Document
#' @title Cluster based on skew-geometric Jensen-Shannon divergence for the
#' conditional extremes model.
#' @description Function to calculate and cluster on the skew-geometric
#' Jensen-Shannon Divergence for the conditional extremes model.
#' When comparing conditional extremes fits for a single variable using the
#' JS-G between any two locations, we need to use the same "data" for each
#' location. Therefore, we look at values from the maximum of the thresholds
#' at each location, to some multiple `dat_max_mult` of this, with
#' `n_dat` equally spaced points.
#' \code{cluster:pam} is used to cluster the Jensen-Shannon divergence distance
#' matrix between all locations, summed across all variables.
#' @param x Object of class `cecl_dep`.
#' @param marg_obj Object of class `cecl_marg`.
#' @param k Number of clusters.
#' @param cond_var Optional conditioning variable name to cluster
#' on, if not all variables, Default: NULL.
#' @param cluster_mem Optional vector of true cluster memberships
#' to evaluate clustering solution, Default: NULL.
#' @inheritParams cecl_dist
#' @return List containing the clustering results and, if `cluster_mem` is
#' provided, the adjusted Rand index.
#' @rdname cecl_clust
#' @export
#' @method cecl_clust cecl_dep
cecl_clust.cecl_dep <- \(
  x,
  marg_obj,
  k,
  cond_var = NULL,
  cluster_mem = NULL,
  laplace_cap = 0.99,
  n_mc = 500,
  ncores = 1,
  par_dist = FALSE,
  seed = NULL,
  ...
) {
  stopifnot(inherits(x, "cecl_dep"))
  stopifnot(inherits(marg_obj, "cecl_marg"))
  
  # Calculate distance matrix
  dist_obj <- cecl_dist(
    dep_obj     = x,
    marg_obj    = marg_obj,
    
    # Previously this was passed as `cond_var`, but cecl_dist() calls the
    # conditioning-variable selector `var`.
    # cond_var    = cond_var,
    var          = cond_var,
    
    laplace_cap = laplace_cap,
    n_mc        = n_mc,
    ncores      = ncores,
    par_dist    = par_dist,
    seed        = seed,
    ...
  )
  
  # Cluster based on distance matrix
  ret <- cecl_clust.cecl_dist(
    x = dist_obj,
    k = k,
    cluster_mem = cluster_mem,
    ...
  )
  
  # for cecl_dep method, also return dist obj as user may want to inspect
  ret <- c(
    ret,
    list("dist_obj" = dist_obj)
  )
  
  class(ret) <- c("cecl_clust", class(ret))
  return(ret)
}

#' @title Cluster based on skew-geometric Jensen-Shannon divergence
#' distance matrix for the conditional extremes model.
#' @description Function to cluster on the skew-geometric
#' Jensen-Shannon Divergence distance matrix for the conditional
#' extremes model.
#' @param x Object of class `cecl_dist`.
#' @inheritParams cecl_clust.cecl_dep
#' @return List containing the clustering results and, if `cluster_mem` is
#' provided, the adjusted Rand index.
#' @rdname cecl_clust
#' @export
#' @method cecl_clust cecl_dist
cecl_clust.cecl_dist <- \(
  x,
  k,
  cond_var = NULL,
  cluster_mem = NULL,
  ...
) {
  stopifnot(inherits(x, "cecl_dist"))
  
  dist_mat <- x$dist_mat
  if (!is.null(cond_var)) {
    if (!cond_var %in% names(x$dist_mats)) {
      stop(
        sprintf(
          "Conditioning variable '%s' is not available in `x$dist_mats`.",
          cond_var
        ),
        call. = FALSE
      )
    }
    dist_mat <- x$dist_mats[[cond_var]]
  }
  
  # fit clustering
  ret <- cluster::pam(dist_mat, k = k)
  
  # add specific distance matrix used to output
  ret <- list(
    "pam" = ret,
    "dist_mat" = dist_mat,
    "dep_df" = x$dep_df
  )
  
  # evaluate quality of clustering solution if true membership provided
  if (!is.null(cluster_mem)) {
    adj_rand <- mclust::adjustedRandIndex(
      ret$pam$clustering,
      cluster_mem
    )
    ret$adj_rand <- adj_rand
  }
  
  class(ret) <- c("cecl_clust", class(ret))
  return(ret)
}

#' @title Print summary of `cecl_clust` object
#' @description Print a summary of a fitted `cecl_clust` object.
#' @param x Object of class `cecl_clust`.
#' @param ... Additional arguments (not used).
#' @rdname print.cecl_clust
#' @export
#' @method print cecl_clust
print.cecl_clust <- \(x, ...) {
  stopifnot(inherits(x, "cecl_clust"))
  cat("Clustering results of class 'cecl_clust'\n")
  cat("Number of locations:", attr(x$dist_mat, "Size"), "\n")
  cat("Number of clusters:", length(unique(x$pam$clustering)), "\n")
  if (!is.null(x$adj_rand)) {
    cat("Adjusted Rand index:", round(x$adj_rand, 4), "\n")
  }
  invisible(x)
}

# TODO Add model/matrix used in clustering to cecl_clust object
# TODO Change `pam` name in cluster object to `cluster` (to be more abstract)
#' @title Extract clustering solution from `cecl_clust` object
#' @description Extract clustering solution from a fitted
#' `cecl_clust` object.
#' @param object Object of class `cecl_clust`.
#' @param ... Additional arguments (not used).
#' @return Data frame of cluster allocations for each location.
#' @rdname coef.cecl_clust
#' @export
#' @method coef cecl_clust
coef.cecl_clust <- \(object, ...) {
  stopifnot(inherits(object, "cecl_clust"))
  
  ret <- data.frame(object$pam$clustering)
  ret$name <- row.names(ret)
  names(ret)[1] <- "cluster"
  rownames(ret) <- NULL
  
  class(ret) <- c("coef.cecl_clust", class(ret))
  return(ret)
}

#' @title Summary of `cecl_clust` object
#' @description Summarise a fitted `cecl_clust` object.
#' @param object Object of class `cecl_clust`.
#' @param ... Additional arguments (not used).
#' @return Summary of clustering results.
#' @rdname summary.cecl_clust
#' @export
#' @method summary cecl_clust
summary.cecl_clust <- \(object, ...) {
  stopifnot(inherits(object, "cecl_clust"))
  summary(object$pam)
}

# TODO Do I want `...` to proxy::dist or jsg_div??
# TODO Need details section here with extra info on laplace_cap and MC sampling
#' @title Calculate skew-geometric Jensen-Shannon divergence distance matrix
#' for the conditional extremes model.
#' @description Function to calculate the skew-geometric
#' Jensen-Shannon Divergence distance matrix for the conditional
#' extremes model.
#' @param dep_obj Object of class `cecl_dep`.
#' @param marg_obj Object of class `cecl_marg`.
#' @param var Optional conditioning variable name to calculate distance matrix
#' for, if not all variables. Default: NULL.
#' @param laplace_cap Upper quantile used to determine the Laplace-distribution
#' truncation point. Default: 0.99.
#' @param laplace_cap_val Optional fixed upper truncation value on the Laplace
#' scale. If supplied, this is used instead of calculating the upper value from
#' `laplace_cap`.
#' @param laplace_sample Optional numeric vector of pre-generated conditioning
#' values.
#' @param n_mc Number of Monte Carlo samples to use in distance
#' calculation. Default: 500.
#' @param ncores Number of cores to use for parallel computation. Default: 1.
#' @param par_dist Logical, whether to parallelise the distance computation
#' rather than parallelising over conditioning variables. Default: FALSE.
#' @param seed Seed number for random number generation in Monte Carlo sampling.
#' @param ... Additional arguments passed to the divergence calculation.
#' @return Object of class `cecl_dist`.
#' @rdname cecl_dist
#' @export
# TODO Add lambda parameter for weights in distance calculation
cecl_dist <- \(
  dep_obj,
  marg_obj,
  var = NULL,
  laplace_cap = 0.99,
  laplace_cap_val = NULL,
  laplace_sample = NULL,
  n_mc = 500,
  ncores = 1,
  par_dist = FALSE,
  seed = NULL,
  ...
) {
  stopifnot(inherits(dep_obj, "cecl_dep"))
  stopifnot(inherits(marg_obj, "cecl_marg"))
  
  # laplace_cap determines a data quantile, while laplace_cap_val supplies
  # the corresponding Laplace-scale value directly. Because laplace_cap has
  # a non-NULL default, supplying laplace_cap_val overrides it below.
  if (!is.null(laplace_cap) &&
      (!is.numeric(laplace_cap) ||
       length(laplace_cap) != 1L ||
       is.na(laplace_cap) ||
       laplace_cap <= 0 ||
       laplace_cap >= 1)
  ) {
    stop(
      "`laplace_cap` must be a single probability strictly between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (!is.null(laplace_cap_val) &&
      (!is.numeric(laplace_cap_val) ||
       length(laplace_cap_val) != 1L ||
       is.na(laplace_cap_val))
  ) {
    stop(
      "`laplace_cap_val` must be a single numeric value.",
      call. = FALSE
    )
  }
  
  # Only want a single conditioning variable, if provided.
  stopifnot(
    "`var` must be NULL or a single conditioning-variable name." =
      is.null(var) || length(var) == 1L
  )
  
  # check that laplace_sample is correct, if provided
  if (!is.null(laplace_sample)) {
    stopifnot(
      "`laplace_sample` must be a numeric vector." =
        is.vector(laplace_sample) && is.numeric(laplace_sample)
    )
    stopifnot(
      "`laplace_sample` must be a numeric vector of length `n_mc`." =
        length(laplace_sample) == n_mc
    )
  }
  
  # pull transformed data
  trans <- marg_obj$transformed
  
  # <<<<<<< Updated upstream
  #   dependence <- dep_obj$dependence
  #
  #   # pull parameter values for each location
  #   params <- lapply(dep_obj$dependence, pull_params)
  #
  #   # pull Laplace-scale threshold values for each location
  #   thresh <- lapply(dep_obj$dependence, pull_thresh_trans)
  #
  #   # list of locs containing vars -> list of vars, each containing all locs
  #   params <- purrr::transpose(params)
  #
  #   # If only want a single *conditioning* variable
  #   if (!is.null(var)) {
  #     params <- params[var]
  #     thresh <- thresh[var]
  #   }
  #
  #   # take maximum Laplace thresholds; want to generate points above this
  #   thresh_max <- lapply(dplyr::bind_rows(thresh), max)
  # =======
  #   # pull parameter values for each location
  #   params <- lapply(dep_obj$dependence, pull_params)
  #
  #   # list of locs containing vars -> list of vars, each containing all locs
  #   params <- purrr::transpose(params)
  #
  #   # pull Laplace-scale threshold values for each location
  #   if (is.null(dth)) {
  #     thresh <- lapply(dep_obj$dependence, pull_thresh_trans)
  #   } else {
  #     # TODO Maybe ensure this is done by variable?? Otherwise will lead to error
  #     thresh <- dth
  #   }
  #
  #   # If only want specific dependent variables
  #   if (!is.null(cond_var)) {
  #     params <- params[cond_var]
  #     if (is.null(dth)) {
  #       thresh <- lapply(thresh, \(x) x[[cond_var]])
  #     } else {
  #       thresh <- thresh[[cond_var]]
  #     }
  #   }
  #
  #   # take maximum Laplace thresholds; want to generate points above this
  #   if (is.null(dth)) {
  #     thresh_max <- lapply(dplyr::bind_rows(thresh), max)
  #   } else {
  #     thresh_max <- list(dth)
  #   }
  # >>>>>>> Stashed changes
  
  # Pull parameters and thresholds for every location. Each location contains
  # an outer list named by RHS conditioning variable. Each parameter matrix has
  # columns named by LHS conditioned variable.
  params_by_loc <- lapply(dep_obj$dependence, pull_params)
  thresh_by_loc <- lapply(dep_obj$dependence, pull_thresh_trans)
  
  # Transpose from:
  #   location -> conditioning variable
  # to:
  #   conditioning variable -> location
  params <- purrr::transpose(params_by_loc)
  thresh <- purrr::transpose(thresh_by_loc)
  
  if (length(params) == 0L) {
    stop(
      "`dep_obj` does not contain any conditioning-variable fits.",
      call. = FALSE
    )
  }
  
  if (is.null(names(params))) {
    stop(
      "The dependence fits must be named by conditioning variable.",
      call. = FALSE
    )
  }
  
  # If only one conditioning variable is requested, retain its outer element.
  if (!is.null(var)) {
    if (!var %in% names(params)) {
      stop(
        sprintf(
          "Conditioning variable '%s' is not available in `dep_obj`.",
          var
        ),
        call. = FALSE
      )
    }
    
    params <- params[var]
    thresh <- thresh[var]
  }
  
  # For each conditioning variable, use the maximum fitted threshold across
  # all locations and all corresponding conditioned-variable columns.
  thresh_max <- lapply(thresh, \(thresholds_i) {
    threshold_values <- unlist(thresholds_i, use.names = FALSE)
    
    if (length(threshold_values) == 0L ||
        all(is.na(threshold_values))
    ) {
      stop(
        "No non-missing thresholds were found for a conditioning variable.",
        call. = FALSE
      )
    }
    
    max(threshold_values, na.rm = TRUE)
  })
  
  # Sample values from a Laplace upper tail truncated above at y_max.
  rlaplace_trunc <- \(n, thresh_max, trans_x, upper_quant = 0.99, y_max) {
    if (is.null(y_max)) {
      y_max <- stats::quantile(
        trans_x,
        upper_quant,
        na.rm = TRUE,
        names = FALSE
      )
    }
    
    stopifnot(
      "y_max must be greater than thresh_max" = y_max > thresh_max
    )
    
    # Above a positive Laplace threshold, the excess follows an exponential
    # distribution. Truncate the excess at y_max - thresh_max.
    p_max <- 1 - exp(-(y_max - thresh_max))
    U <- stats::runif(n, min = 0, max = p_max)
    W <- -log(1 - U)
    
    list(
      "y" = thresh_max + W,
      "y_max" = y_max
    )
  }
  
  # optionally set seed to ensure reproducibility
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      has_seed <- TRUE
    } else {
      has_seed <- FALSE
    }
    
    set.seed(seed)
    
    on.exit(
      {
        if (has_seed) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else if (exists(
          ".Random.seed",
          envir = .GlobalEnv,
          inherits = FALSE
        )) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
  }
  
  # Generate common RHS conditioning values for every location fitted for each
  # conditioning variable.
  y <- lapply(seq_along(thresh_max), \(i) {
    conditioning_var_i <- names(thresh_max)[[i]]
    
    if (!all(vapply(
      trans,
      \(trans_i) conditioning_var_i %in% colnames(trans_i),
      logical(1)
    ))) {
      stop(
        sprintf(
          "Conditioning variable '%s' is absent from transformed data.",
          conditioning_var_i
        ),
        call. = FALSE
      )
    }
    
    trans_x <- unlist(
      lapply(
        trans,
        \(trans_i) trans_i[, conditioning_var_i, drop = TRUE]
      ),
      use.names = FALSE
    )
    
    if (is.null(laplace_sample)) {
      res <- rlaplace_trunc(
        n = n_mc,
        thresh_max = thresh_max[[i]],
        trans_x = trans_x,
        upper_quant = laplace_cap,
        y_max = laplace_cap_val
      )
      
      # Check whether the shared upper value exceeds any site-wise maxima.
      site_max <- vapply(
        trans,
        \(trans_i) {
          max(trans_i[, conditioning_var_i, drop = TRUE], na.rm = TRUE)
        },
        numeric(1)
      )
      n_above <- sum(res$y_max > site_max)
      
      if (n_above > 0L) {
        message(paste0(
          "Extrapolation performed for ",
          n_above,
          " groups for conditioning variable ",
          conditioning_var_i,
          ". Consider using a lower `laplace_cap`, ",
          "or a minima-of-site-wise-maxima approach."
        ))
      }
      
      return(res$y)
    }
    
    laplace_sample
  })
  names(y) <- names(params)
  
  # also want to return maximum value in y for each variable
  y_max <- vapply(y, max, numeric(1))
  
  # calculate distance matrices for each conditioning variable
  dist_chunk <- \(chunk_idx, lst, thresh_max, n, y_spec) {
    proxy::dist(
      lst[chunk_idx],
      lst,
      method     = jsg_div,
      thresh_max = thresh_max,
      n          = n,
      y_spec     = y_spec,
      upper      = TRUE,
      ...        = ...
    )
  }
  
  # parallel setup for embarrassingly parallel mclapply use
  apply_fun <- if (ncores == 1L) lapply else parallel::mclapply
  ext_args <- NULL
  if (ncores > 1L) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }
  
  # if not paralleling distance computation
  if (par_dist == FALSE) {
    dist_mats <- loop_fun(seq_along(params), \(i) {
      mat <- dist_chunk(
        chunk_idx = seq_len(length(params[[i]])),
        lst = params[[i]],
        thresh_max = thresh_max[[i]],
        n = n_mc,
        y_spec = y[[i]]
      )
      
      class(mat) <- "matrix"
      stats::as.dist(mat)
    })
  } else {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      varlist = c("dist_chunk", "jsg_div", "jsg_gauss"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, library(proxy))
    
    nlocs <- length(params[[1]])
    chunks <- split(
      seq_len(nlocs),
      sort(rep(seq_len(ncores), length.out = nlocs))
    )
    
    dist_mats <- lapply(seq_along(params), \(i) {
      mat <- do.call(
        rbind,
        parallel::parLapply(cl, chunks, \(chunk_idx) {
          dist_chunk(
            chunk_idx = chunk_idx,
            lst = params[[i]],
            thresh_max = thresh_max[[i]],
            n = n_mc,
            
            # Previously this passed the complete list `y`, which gave the
            # divergence function the wrong object.
            # y_spec = y,
            y_spec = y[[i]]
          )
        })
      )
      
      class(mat) <- "matrix"
      stats::as.dist(mat)
    })
    
    parallel::stopCluster(cl)
    cl <- NULL
  }
  
  # average distance matrices over conditioning variables
  dist_mat <- Reduce(`+`, dist_mats) / length(dist_mats)
  
  # name distance matrices by RHS conditioning variable
  names(dist_mats) <- names(params)
  
  dist_ret <- list(
    "dist_mat"  = dist_mat,
    "dist_mats" = dist_mats,
    "y"         = y,
    "y_max"     = y_max,
    "call"      = match.call(),
    "dep_df"    = stats::coef(dep_obj)
  )
  
  class(dist_ret) <- c("cecl_dist", class(dist_ret))
  dist_ret
}

#' @title Print summary of `cecl_dist` object
#' @description Print a summary of a fitted `cecl_dist` object.
#' @param x Object of class `cecl_dist`.
#' @param ... Additional arguments (not used).
#' @rdname print.cecl_dist
#' @export
#' @method print cecl_dist
print.cecl_dist <- \(x, ...) {
  stopifnot(inherits(x, "cecl_dist"))
  cat("Distance matrix of class 'cecl_dist'\n")
  cat("Number of locations:", attr(x$dist_mat, "Size"), "\n")
  cat(
    "Conditioning variables:",
    paste(names(x$dist_mats), collapse = ", "),
    "\n"
  )
  invisible(x)
}

#' @title Summary of `cecl_dist` object
#' @description Summarise a fitted `cecl_dist` object.
#' @param object Object of class `cecl_dist`.
#' @param var Optional conditioning-variable name to summarise, if not using
#' the overall distance matrix. Default: NULL.
#' @param ... Additional arguments (not used).
#' @return Summary of distance matrix.
#' @rdname summary.cecl_dist
#' @export
#' @method summary cecl_dist
summary.cecl_dist <- \(object, var = NULL, ...) {
  stopifnot(inherits(object, "cecl_dist"))
  
  dist_mat <- object$dist_mat
  if (!is.null(var)) {
    if (!var %in% names(object$dist_mats)) {
      stop(
        sprintf(
          "Conditioning variable '%s' is not available.",
          var
        ),
        call. = FALSE
      )
    }
    dist_mat <- object$dist_mats[[var]]
  }
  
  summary(dist_mat)
}