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
#' @param var Optional conditioning variable name to cluster
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
  var = NULL,
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
    var         = var,
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
  var = NULL,
  cluster_mem = NULL,
  ...
) {
  stopifnot(inherits(x, "cecl_dist"))

  dist_mat <- x$dist_mat
  if (!is.null(var)) {
    dist_mat <- x$dist_mats[[var]]
  }

  # fit clustering
  ret <- cluster::pam(dist_mat, k = k)

  # add specific distance matrix used to output
  # TODO Need for three calls to list here?
  ret <- c(
    list("pam" = ret),
    list("dist_mat" = dist_mat),
    list("dep_df" = x$dep_df) # output for plot methods
  )

  # evaluate quality of clustering solution if true membership provided
  if (!is.null(cluster_mem)) {
    adj_rand <- mclust::adjustedRandIndex(
      ret$pam$clustering,
      cluster_mem
    )
    ret <- c(ret, "adj_rand" = adj_rand)
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
# TODO Need details section here with extra info on laplace_cap and MC
# sampling
#' @title Calculate skew-geometric Jensen-Shannon divergence distance matrix
#' for the conditional extremes model.
#' @description Function to calculate the skew-geometric
#' Jensen-Shannon Divergence distance matrix for the conditional
#' extremes model.
#' @param dep_obj Object of class `cecl_dep`.
#' @param marg_obj Object of class `cecl_marg`.
#' @param var Optional conditioning variable name to calculate distance matrix
#' for, if not all variables, Default: NULL.
#' @param laplace_cap Upper quantile to sample Laplace distribution
#' truncation point from, Default: 0.99.
#' @param n_mc Number of Monte Carlo samples to use in distance
#' calculation, Default: 500.
#' @param ncores Number of cores to use for parallel computation, Default: 1.
#' @param par_dist Logical, whether to parallelise the distance computation,
#' rather than embarrassingly parallel over each variable, which is useful for a
#' high number of locations, but has significant overheads, Default: FALSE.
#' @param seed Seed number for random number generation in Monte Carlo sampling,
#' for reproducibility.
#' @param ... Additional arguments passed to `proxy::dist()`.
#' @return List of distance matrices for each variable.
#' @rdname cecl_dist
#' @export
# TODO Add lambda parameter for weights in distance calculation
cecl_dist <- \(
  dep_obj,
  marg_obj,
  var = NULL,
  laplace_cap = 0.99,
  n_mc = 500,
  ncores = 1,
  par_dist = FALSE,
  seed = NULL,
  ...
) {
  stopifnot(inherits(dep_obj, "cecl_dep"))
  stopifnot(inherits(marg_obj, "cecl_marg"))

  n <- NULL

  # Only want a single variable, if provided
  stopifnot(is.null(var) || length(var == 1))

  # pull transformed data
  trans <- marg_obj$transformed

  dependence <- dep_obj$dependence

  # pull parameter values for each location
  params <- lapply(dep_obj$dependence, pull_params)

  # pull Laplace-scale threshold values for each location
  thresh <- lapply(dep_obj$dependence, pull_thresh_trans)

  # list of locs containing vars -> list of vars, each containing all locs
  params <- purrr::transpose(params)

  # If only want a single *conditioning* variable
  if (!is.null(var)) {
    params <- params[var]
    thresh <- thresh[var]
  }

  # take maximum Laplace thresholds; want to generate points above this
  thresh_max <- lapply(dplyr::bind_rows(thresh), max)

  # TODO Move calculating y values to separate function? As above
  # TODO Move this to separate function anyway!
  rlaplace_trunc <- \(n, thresh_max, trans_x, upper_quant = 0.99) {
    # get maximum point
    y_max <- stats::quantile(trans_x, upper_quant, na.rm = TRUE)
    stopifnot(
      "y_max must be greater than thresh_max" = y_max > thresh_max
    )
    # get probability of being below this point from exponential CDF
    p_max <- 1 - exp(-(y_max - thresh_max))
    # sample from uniform distribution below this point
    U <- stats::runif(n, min = 0, max = p_max) # min=0 as we push up by thresh
    # inversion sampling from exponential distribution
    W <- -log(1 - U)
    # shift to the right by the threshold to get samples from truncated Laplace
    return(list(
      "y"     = thresh_max + W,
      "y_max" = y_max # to check if y_max exceeds any site-wise maxima
    ))
  }

  # optionally set seed to ensure reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # loop through variables
  y <- lapply(seq_along(thresh_max), \(i) {
    # get transformed data for this variable
    trans_x <- unlist(lapply(trans, \(x) x[, i, drop = TRUE]))
    # sample from truncated Laplace distribution
    res <- rlaplace_trunc(
      n_mc, thresh_max[[i]], trans_x,
      upper_quant = laplace_cap
    )

    # check if y_max exceeds any site-wise maxima
    # If so, will be performing extrapolation, may want to warn user
    site_max <- vapply(trans, \(x) max(x[, i], na.rm = TRUE), numeric(1))
    n_above <- sum(res$y_max > site_max)
    if (n_above > 0) {
      message(paste0(
        "Extrapolation performed for ",
        n_above,
        " groups for variable ",
        names(thresh_max)[i],
        ". Consider using a lower `laplace_cap`, ",
        " or 'minima of site-wise maxima' approach."
      ))
    }

    res$y
  })
  names(y) <- names(params)

  # also want to return maximum value in y for each variable
  y_max <- vapply(y, max, numeric(1))

  # calculate distance matrices for each variable
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
  # TODO Functionalise, used elsewhere
  apply_fun <- ifelse(ncores == 1, lapply, parallel::mclapply)
  ext_args <- NULL
  if (ncores > 1) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }

  # if not paralleling distance computation (reserved for high # locations)
  if (par_dist == FALSE) {
    dist_mats <- loop_fun(seq_along(params), \(i) {
      mat <- dist_chunk(
        seq_len(length(params[[i]])), # don't chunk
        params[[i]],
        thresh_max[[i]],
        n = n,
        y_spec = y[[i]],
        ... = ...
      )
      # convert from crossidst to dist object
      class(mat) <- "matrix"
      return(stats::as.dist(mat))
    })
    # if paralleling distance computation, need to split into chunks
  } else {
    cl <- parallel::makeCluster(ncores)
    # Export cluster objects
    parallel::clusterExport(
      cl,
      varlist = c("dist_chunk", "jsg_div"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, library(proxy))
    # split distance calcuation for locations into chunks
    nlocs <- length(params[[1]])
    chunks <- split(
      seq_len(nlocs),
      sort(rep(seq_len(ncores), length.out = nlocs))
    )

    # Compute distances in parallel
    dist_mats <- lapply(seq_along(params), \(i) {
      mat <- do.call(rbind, parallel::parLapply(cl, chunks, \(chunk_idx) {
        dist_chunk(
          chunk_idx,
          lst = params[[i]],
          thresh_max = thresh_max[[i]],
          n,
          y_spec = y
        )
      }))
      # convert from cross.dist to dist
      class(mat) <- "matrix"
      return(stats::as.dist(mat))
    })
    parallel::stopCluster(cl) # stop cluster
  }

  # average distance matrices over different variables together
  dist_mat <- Reduce(`+`, dist_mats) / length(dist_mats)

  # name distance matrices
  names(dist_mats) <- names(params)

  # list to output
  dist_ret <- list(
    "dist_mat"  = dist_mat,
    "dist_mats" = dist_mats,
    "y"         = y,
    "y_max"     = y_max,
    "call"      = match.call(),
    "dep_df"    = stats::coef(dep_obj) # include dep params for plot methods
  )

  class(dist_ret) <- c("cecl_dist", class(dist_ret))
  dist_ret
}

#' @ title Print summary of `cecl_dist` object
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
  cat("Variables:", paste(names(x$dist_mats), collapse = ", "), "\n")
  invisible(x)
}

#' @title Summary of `cecl_dist` object
#' @description Summarise a fitted `cecl_dist` object.
#' @param object Object of class `cecl_dist`.
#' @param var Optional variable name to summarise distance matrix for,
#' if not overall distance matrix, Default: NULL.
#' @param ... Additional arguments (not used).
#' @return Summary of distance matrix.
#' @rdname summary.cecl_dist
#' @export
#' @method summary cecl_dist
summary.cecl_dist <- \(object, var = NULL, ...) {
  stopifnot(inherits(object, "cecl_dist"))
  dist_mat <- object$dist_mat
  if (!is.null(var)) {
    dist_mat <- object$dist_mats[[var]]
  }
  summary(dist_mat)
}

#' @title Plotting function for `cecl_clust` object
#' @description Create plots from a fitted `cecl_clust` object.
#' @param x Object of class `cecl_clust`.
#' @param which Character string specifying which plot to produce.
#' Either `"image"` for distance matrix image/heatmap, or `"scatter"` for
#' dependence parameters scatter plot, both with labels coloured by
#' cluster membership.
#' @param ... Additional arguments passed to plotting functions.
#' @rdname plot.cecl_clust
#' @export
#' @method plot cecl_clust
plot.cecl_clust <- \(
  x,
  which = c("image", "scatter"),
  ...
) {
  stopifnot(inherits(x, "cecl_clust"))
  which <- match.arg(which)

  if (which == "image") {
    plot_image(x, type = "plot", ...)
  } else if (which == "scatter") {
    plot_scatter(x, type = "plot", ...)
  }
}

#' @title ggplot function for `cecl_clust` object
#' @description Create ggplot objects from a fitted `cecl_clust` object.
#' @param data Object of class `cecl_clust`.
#' @param mapping ggplot2 aesthetic mappings.
#' @inheritParams plot.cecl_clust
#' @method ggplot cecl_clust
#' @export
ggplot.cecl_clust <- \(
  data = NULL,
  mapping = ggplot2::aes(),
  which = c("image", "scatter"),
  ...
) {
  stopifnot(inherits(data, "cecl_clust"))
  which <- match.arg(which)

  if (which == "image") {
    plot_image(data, type = "ggplot", ...)
  } else if (which == "scatter") {
    plot_scatter(data, type = "ggplot", ...)
  }
}

# TODO Keep dep_obj as arg here, or add to cecl_dist and cecl_clust objects?
#' @title Plot scatter plot from `cecl_dep` object
#' @description Plot scatter plot of dependence parameters from a fitted
#' `cecl_dep` object.
#' @param x Object of class `cecl_clust`.
#' @param var Conditioned variable name to plot.
#' @param cond_var Conditioning variable name to plot against.
#' @param labels List mapping variable names to plot labels, e.g.,
#' `list("rain" = "Precipitation", "wind" = "Wind Speed")`, for use in axis
#' labels. Default is `NULL`, which uses variable names as is.
#' @param type Type of plot to return. Either `"ggplot"` (default) or `"plot"`.
#' @param ... Additional arguments to pass to plotting functions.
#' @return ggplot object of scatter plot.
#' @rdname plot_scatter
#' @export
#' @method plot_scatter cecl_clust
plot_scatter.cecl_clust <- \(
  x, var, cond_var, labels = NULL, type = c("ggplot", "plot"), ...
) {
  stopifnot(inherits(x, "cecl_clust"))
  type <- match.arg(type)
  stopifnot("dep_df" %in% names(x))

  a <- b <- clust <- name <- NULL

  # pull dependence parameters for all locs for pecific var/cond_var
  dep_params_spec <- x$dep_df |>
    dplyr::filter(var == !!var, cond_var == !!cond_var) |>
    # add colour variable based on clustering
    dplyr::left_join(
      data.frame(
        "clust" = factor(x$pam$clustering),
        "name" = names(x$pam$clustering)
      ),
      by = "name"
    )

  # For plotting, tidy up variable names
  var_lab <- var
  cond_var_lab <- cond_var
  if (!is.null(labels)) {
    var_lab <- labels[[var]]
    cond_var_lab <- labels[[cond_var]]
  }
  dep_params_spec$facet_lab <- paste0(var_lab, " | ", cond_var_lab)

  if (type == "ggplot") {
    plot <- dep_params_spec |>
      ggplot2::ggplot(ggplot2::aes(x = a, y = b, colour = clust)) +
      ggplot2::geom_point(...) +
      ggplot2::facet_wrap(~facet_lab) +
      cecl_theme() +
      ggplot2::labs(
        x = expression(a),
        y = expression(b),
      ) +
      ggplot2::guides(colour = "none")

    # add labels if ggrepel is installed
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot <- plot +
        ggrepel::geom_text_repel(ggplot2::aes(label = name))
    } else {
      plot <- plot +
        ggplot2::geom_text(ggplot2::aes(label = name), vjust = -0.5)
    }

    return(plot)
    # base plot
  } else {
    cl_levels <- levels(dep_params_spec$clust)
    n_cl <- length(cl_levels)
    if (n_cl < 9) {
      cols <- ggsci::pal_nejm()(n_cl)
      # fallback if there are too many colours
    } else {
      cols <- grDevices::rainbow(n_cl)
    }

    # add alpha/transparency to match ggplot semi-transparent points
    cols_alpha <- grDevices::adjustcolor(cols, alpha.f = 0.5)

    # map colours to each row by cluster
    col_map <- stats::setNames(cols_alpha, cl_levels)
    point_cols <- col_map[as.character(dep_params_spec$clust)]

    # plot points coloured by cluster
    plot(
      dep_params_spec$a,
      dep_params_spec$b,
      xlab = expression(alpha),
      ylab = expression(beta),
      main = paste0(var_lab, " | ", cond_var_lab),
      pch = 16,
      col = point_cols,
      xlim = c(-1, 1),
      ylim = c(min(dep_params_spec$b) - 0.1, max(dep_params_spec$b) + 0.1),
      ...
    )

    # add labels, using same colour as points (but slightly darker for contrast)
    label_cols <- grDevices::adjustcolor(cols, alpha.f = 1)[
      as.integer(dep_params_spec$clust)
    ]
    graphics::text(
      dep_params_spec$a,
      dep_params_spec$b,
      labels = dep_params_spec$name,
      pos = 3,
      col = label_cols
    )
  }
}

# TODO Change to be S3 method like others?
#' @title Scree plot for `cecl_dist` object
#' @description Create scree plot of total within-group sum of squares
#' for different k values from a fitted `cecl_dist` object.
#' @param dist_mat Distance matrix of class `dist`.
#' @param scree_k Vector of k values to use in scree plot, Default is `1:5`.
#' @param type Character string specifying which plot to produce.
#' Either `"ggplot"` for ggplot object, or `"plot"` for base R plot.
#' Default is `"ggplot"`.
#' @param ... Additional arguments passed to plotting functions.
#' @return ggplot object or base R plot of scree plot, and total within-group
#' sum of squares values.
#' @rdname plot_scree
#' @keywords internal
# plot_scree <- \(dist_obj, scree_k = 1:5, type = c("ggplot", "plot"), ...) {
plot_scree <- \(dist_mat, scree_k = 1:5, type = c("ggplot", "plot"), ...) {
  # stopifnot(inherits(dist_obj, "cecl_dist"))
  type <- match.arg(type)

  k <- twgss <- NULL

  # dist_mat <- dist_obj$dist_mat

  total_within_ss <- vapply(
    scree_k, within_cluster_sum, dist_mat,
    fun = cluster::pam, FUN.VALUE = numeric(1)
  )

  if (type == "plot") {
    plot(
      scree_k,
      total_within_ss,
      type = "b",
      xlab = "k",
      ylab = "TWGSS",
      ...
    )
    return(total_within_ss)
  } else {
    df <- data.frame(
      "k"     = scree_k,
      "twgss" = total_within_ss
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = twgss)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      cecl_theme() +
      ggplot2::labs(
        x = "k",
        y = "TWGSS"
      )
    return(list(
      "plot" = p,
      "twgss" = total_within_ss
    ))
  }
}

#' @title Generic image plot function
#' @description Generic image plot function for different object classes.
#' @param x Object to plot.
#' @param ... Additional arguments passed to methods.
#' @return Plot of object.
#' @rdname plot_image
#' @keywords internal
plot_image <- \(x, ...) {
  UseMethod("plot_image")
}

#' @title Image plot for `cecl_dist` object
#' @description Create image/heatmap plot of distance matrix
#' from a fitted `cecl_dist` object.
#' @param x Object of class `cecl_dist`.
#' @param type Character string specifying which plot to produce.
#' Either `"ggplot"` for ggplot object, or `"plot"` for base R plot.
#' Default is `"ggplot"`.
#' @param ... Additional arguments passed to plotting functions.
#' @return ggplot object or base R plot of distance matrix image/heatmap.
# #' @rdname plot_image
# #' @method plot_image cecl_dist
#' @export
plot_image.cecl_dist <- \(x, type = c("ggplot", "plot"), ...) {
  stopifnot(inherits(x, "cecl_dist"))
  type <- match.arg(type)

  Var1 <- Var2 <- Freq <- NULL

  dist_matrix <- as.matrix(x$dist_mat) # convert from dist to matrix
  # dist_matrix <- as.matrix(dist_mat)

  if (type == "plot") {
    # extract row and column tick labels
    x_names <- rownames(dist_matrix)
    y_names <- rev(colnames(dist_matrix))
    n <- length(x_names)

    dist_plt <- t(dist_matrix[x_names, y_names])

    # plot
    graphics::image(
      1:n, 1:n, dist_plt,
      axes = FALSE,
      zlim = c(0, max(dist_matrix)), # ensure 0 included
      xlab = "", ylab = "",
      ...
    )

    graphics::axis(1, at = 1:n, labels = x_names, las = 2)
    graphics::axis(2, at = 1:n, labels = y_names, las = 2)
    graphics::box()
  } else {
    # convert matrix to data frame
    df <- as.data.frame(as.table(dist_matrix)) |>
      # ensure correct ordering of factors
      dplyr::mutate(
        Var1 = factor(Var1, levels = rownames(dist_matrix)),
        Var2 = factor(Var2, levels = rev(colnames(dist_matrix)))
      )
    p <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = Freq)) +
      ggplot2::geom_tile() +
      cecl_theme(nejm_pal = FALSE) +
      ggplot2::labs(
        x = "",
        y = "",
        fill = "Distance"
      ) +
      # TODO Change colour scheme to one in paper?
      ggplot2::scale_fill_viridis_c() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    return(p)
  }
}

# TODO Add some kind of highlighting for locations within the same cluster?
# TODO Allow custom upper bound to colour scale? Large distances may dominate
#' @title Image plot for `cecl_clust` object
#' @description Create image/heatmap plot of distance matrix
#' from a fitted `cecl_clust` object.
#' @param x Object of class `cecl_clust`.
#' @param type Character string specifying which plot to produce.
#' Either `"ggplot"` for ggplot object, or `"plot"` for base R plot.
#' Default is `"ggplot"`.
#' @param ... Additional arguments passed to plotting functions.
#' @return ggplot object or base R plot of distance matrix image/heatmap.
#' @rdname plot_image
#' @method plot_image cecl_clust
#' @export
plot_image.cecl_clust <- \(
  x,
  type = c("ggplot", "plot"),
  order_by_cluster = TRUE,
  label_colours = NULL,
  ...
) {
  type <- match.arg(type)
  stopifnot(inherits(x, "cecl_clust"))

  Var1 <- Var2 <- value <- Freq <- NULL

  #  extract distance matrix and PAM clustering
  dist_matrix <- as.matrix(x$dist_mat)
  clustering <- x$pam$clustering

  # ensure clustering is named and matches matrix dimnames if possible
  if (is.null(names(clustering))) {
    # try to attach names from rownames of dist matrix if lengths match
    if (length(clustering) == nrow(dist_matrix)) {
      names(clustering) <- rownames(dist_matrix)
    }
  }

  # if clustering names exist but don't match the matrix names, try to align
  if (!is.null(names(clustering))) {
    common <- intersect(names(clustering), rownames(dist_matrix))
    if (length(common) == 0) {
      message(paste(
        "Cluster membership names do not match distance-matrix rownames.",
        "Will assume same order."
      ))
    } else if (!identical(common, rownames(dist_matrix))) {
      # reorder matrix to the ordering of clustering (for safety)
      keep <- names(clustering)[names(clustering) %in% rownames(dist_matrix)]
      if (length(keep) == nrow(dist_matrix)) {
        dist_matrix <- dist_matrix[keep, keep, drop = FALSE]
        clustering <- clustering[keep]
      } else {
        # partial overlap: keep intersection order of matrix
        ord <- rownames(dist_matrix)
        clustering <- clustering[ord]
      }
    }
  } else {
    # no names -- assume same order
    if (length(clustering) != nrow(dist_matrix)) {
      stop(paste(
        "Length of clustering vector does not match distance-matrix",
        "dimensions and clustering has no names."
      ))
    }
  }

  # optional ordering by cluster
  if (order_by_cluster) {
    # preserve original order within clusters
    ord <- order(clustering, seq_along(clustering))
    dist_matrix <- dist_matrix[ord, ord, drop = FALSE]
    clustering <- clustering[ord]
  } else {
    # ensure clustering is in same order as dist_matrix rows
    if (!is.null(names(clustering))) {
      clustering <- clustering[rownames(dist_matrix)]
    }
  }

  # prepare label colours
  clusters <- sort(unique(as.vector(clustering)))
  k <- length(clusters)
  if (is.null(label_colours)) {
    # default palette
    if (k < 9) {
      palette_cols <- ggsci::pal_nejm()(k)
    } else {
      palette_cols <- grDevices::rainbow(k)
    }
  } else {
    if (length(label_colours) < k) {
      stop("label_colours must provide at least one colour per cluster.")
    }
    palette_cols <- label_colours[seq_len(k)]
  }
  # map cluster id to colour
  cluster_to_col <- stats::setNames(palette_cols, clusters)

  # plotting
  if (type == "plot") {
    # base R image plot and coloured axis text (drawn with text())
    x_names <- rownames(dist_matrix)
    y_names <- rev(colnames(dist_matrix))
    n <- length(x_names)

    # plotting matrix (match ggplot geom_tile orientation)
    dist_plt <- t(dist_matrix[x_names, y_names, drop = FALSE])

    # compute plotting zlim
    zlim <- c(min(dist_matrix, na.rm = TRUE), max(dist_matrix, na.rm = TRUE))

    # draw image without axes
    graphics::image(1:n, 1:n, dist_plt,
      axes = FALSE,
      zlim = zlim,
      xlab = "", ylab = "",
      ...
    )
    # TODO Convert these comments from Chat GPT
    # Draw axes first (so we keep tick marks and baseline label placement)
    graphics::axis(1, at = 1:n, labels = FALSE, las = 2) # ticks only
    graphics::axis(2, at = 1:n, labels = FALSE, las = 2)

    usr <- graphics::par("usr")
    # offsets for label placement
    x_off <- 0.8
    y_off <- 0.3

    # label colors from cluster assignment
    x_cols <- cluster_to_col[as.character(clustering[x_names])]
    # y_sample_names <- rev(y_names)
    # y_cols <- cluster_to_col[as.character(clustering[y_sample_names])]
    y_cols <- cluster_to_col[as.character(clustering[y_names])]

    # Draw x-axis labels (colored) — slightly below the ticks
    graphics::text(
      x = 1:n,
      y = usr[3] - x_off,
      labels = x_names,
      srt = 90,
      adj = 1,
      xpd = TRUE,
      col = x_cols,
      cex = graphics::par("cex.axis")
    )

    # Draw y-axis labels (colored) — slightly left of ticks
    graphics::text(
      x = usr[1] - y_off,
      y = 1:n,
      labels = y_names,
      adj = 1,
      xpd = TRUE,
      col = y_cols,
      cex = graphics::par("cex.axis")
    )

    graphics::box()
  } else {
    # preserve the ordering in the matrix: x_names in matrix row order
    x_names <- rownames(dist_matrix)
    y_names <- colnames(dist_matrix)

    # build heatmap dataframe
    df <- as.data.frame(as.table(dist_matrix)) |>
      # ensure factor order is the current matrix order
      dplyr::mutate(
        Var1 = factor(Var1, levels = x_names),
        Var2 = factor(Var2, levels = rev(y_names))
      )

    # add clustering
    clust_df <- dplyr::as_tibble(clustering) |>
      dplyr::rename(cluster = value) |>
      dplyr::mutate(name = names(clustering))

    # match colours to clustering for x and y axis labels
    plot_cols_x <- palette_cols[
      clust_df$cluster[match(levels(df$Var1), clust_df$name)]
    ]
    plot_cols_y <- palette_cols[
      clust_df$cluster[match(levels(df$Var2), clust_df$name)]
    ]

    # Split into diagonal and off-diagonal dataframes
    # (want to have white for NAs without affecting fill legend)
    diag_df <- dplyr::filter(df, Var1 == Var2)
    off_diag_df <- dplyr::filter(df, Var1 != Var2)

    p <- off_diag_df |>
      # ggplot(aes(x = Var1, y = Var2, fill = Distance)) +
      # TODO optionally bin here
      ggplot2::ggplot(ggplot2::aes(x = Var1, y = Var2, fill = Freq)) +
      ggplot2::geom_tile() +
      ggplot2::geom_tile(
        data = diag_df,
        # aes(x = Column, y = Row),
        ggplot2::aes(x = Var1, y = Var2),
        fill = "white",
        show.legend = FALSE
      ) +
      # ggplot2::scale_fill_viridis_d(
      #   option = "A",
      #   direction = -1
      # ) +
      ggplot2::coord_fixed() + # keep squares square
      ggplot2::labs(
        x = "", y = "",
        fill = "Dissimilarity"
      ) +
      ggplot2::theme(
        # remove x-axis labels, as they are repeats
        # axis.text.x      = ggplot2::element_blank(),
        # axis.ticks.x     = ggplot2::element_blank(),
        # Colour by cluster
        axis.text.x = ggplot2::element_text(
          colour = plot_cols_x,
          size = 11.5,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = ggplot2::element_text(
          colour = plot_cols_y,
          size = 11.5
        ),
        panel.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        legend.title = ggplot2::element_text(size = 15),
        legend.text = ggplot2::element_text(size = 14)
      ) +
      NULL

    # colour based on if colours are binned or not
    p <- p +
      ggplot2::scale_fill_viridis_c(
        option = "A",
        direction = -1
      )

    return(p)
  }

  invisible(NULL)
}

#' @title Plot from `cecl_dist` object
#' @description Plot distance matrix or scree plot from a fitted `cecl_dist`
#' object.
#' @param x Object of class `cecl_dist`.
#' @param which Character string specifying which plot to produce.
#' Either `"image"` for distance matrix image/heatmap, or `"scree"` for
#' scree plot.
#' @param var Optional variable name to plot distance matrix for,
#' if not overall distance matrix, Default: NULL.
#' @param scree_k Vector of k values to use in scree plot, if `which` is
#' `scree`. Default is `1:5`.
#' @param ... Additional arguments passed to plotting functions.
#' @return ggplot object of specified plot.
#' @rdname plot.cecl_dist
#' @export
#' @method plot cecl_dist
plot.cecl_dist <- \(
  x,
  which = c("image", "scree"),
  var = NULL,
  scree_k = 1:5,
  ...
) {
  stopifnot(inherits(x, "cecl_dist"))
  which <- match.arg(which)

  if (which == "image") {
    # plot_image(dist_mat, type = "plot", ...)
    plot_image(x, type = "plot", ...)
  } else if (which == "scree") {
    # pull appropriate distance matrix
    dist_mat <- x$dist_mat
    if (!is.null(var)) {
      dist_mat <- x$dist_mats[[var]]
    }
    plot_scree(dist_mat, type = "plot", scree_k = scree_k, ...)
  }
}

#' @title ggplot from `cecl_dist` object
#' @description Create ggplot of distance matrix or scree plot from a fitted
#' `cecl_dist` object.
#' @param data Object of class `cecl_dist`.
#' @param mapping Not used.
#' @param which Character string specifying which plot to produce.
#' Either `"image"` for distance matrix image/heatmap, or `"scree"` for
#' scree plot.
#' @inheritParams plot.cecl_dist
#' @param ... Additional arguments passed to plotting functions.
#' @param environment Not used.
#' @return ggplot object of specified plot.
#' @rdname ggplot.cecl_dist
#' @export
#' @method ggplot cecl_dist
ggplot.cecl_dist <- \(
  data = NULL,
  mapping = ggplot2::aes(),
  which = c("image", "scree"),
  var = NULL,
  scree_k = 1:5,
  ...,
  environment = parent.frame()
) {
  stopifnot(inherits(data, "cecl_dist"))
  which <- match.arg(which)

  if (which == "image") {
    p <- plot_image(data, type = "ggplot", ...)
  } else if (which == "scree") {
    # pull appropriate distance matrix
    dist_mat <- data$dist_mat
    if (!is.null(var)) {
      dist_mat <- data$dist_mats[[var]]
    }
    p <- plot_scree(dist_mat, type = "ggplot", scree_k = scree_k, ...)
  }
  return(p)
}

# TODO Update docs
#' @title Pull parameters for conditional extremes model
#' @description Function to pull the parameters (a, b, m and s) for the
#' conditional extremes model.
#' @param dep List of `mexDependence` objects for each location.
#' @return List of named vectors containing the parameters for each location.
#' @keywords internal
pull_params <- \(dep) {
  # parameters to extract
  pars <- c("a", "b", "m", "s")
  # fail if not list of `mex` objects for single location
  stopifnot(is.list(dep))
  # loop through conditioning variables for single location
  return(lapply(dep, \(x) {
    return(x[rownames(x) %in% pars, , drop = FALSE])
  }))
}

#' @title Pull thresholds for conditional extremes model
#' @description Function to pull the thresholds for the conditional extremes
#' model.
#' @param dep List of `mexDependence` objects for each location.
#' @return List of named vectors containing the thresholds for each location.
#' @keywords internal
pull_thresh_trans <- \(dep) {
  # fail if not list for single location
  stopifnot(is.list(dep))
  # return quantile of transformed data (already calculated in texmex)
  # return(lapply(dep, \(x) x["dth", ]))
  # use setNames so even if x has 1 col, we still output a named vector
  return(lapply(dep, \(x) stats::setNames(x["dth", ], colnames(x))))
}

# TODO Document arguments (and function itself more fully)
# TODO Allow different function for calculating mu and sigma
# TODO Enable for non-diagonal covariance matrices
#' @title Calculate skew-geometric Jensen-Shannon divergence for each data point
#' @description Function to calculate skew-geometric Jensen-Shannon divergence
#' for each data point.
#' @param params_x Parameters for the first Gaussian distribution.
#' @param params_y Parameters for the second Gaussian distribution.
#' @param thresh_max Maximum threshold value across all locations.
#' @param lambda Skew parameter for skew-geometric JS divergence (default 0.5).
#' @param y_spec Values of the conditioning variable at which to evaluate the
#' divergence.
#' @param agg_fun Function to aggregate divergence values across variables
#' (default `mean`).
#' @param ... Additional arguments to pass to `jsg_gauss()`.
#' @return Skew-geometric Jensen-Shannon divergence value.
#' @keywords internal
jsg_div <- \(
  params_x,
  params_y,
  thresh_max,
  lambda = 0.5,
  y_spec,
  agg_fun = mean,
  ...
) {
  # test that input vectors have correct conditional extremes parameters
  stopifnot(is.matrix(params_x) && is.matrix(params_y))
  stopifnot(all(
    c(rownames(params_x), rownames(params_y)) == rep(c("a", "b", "m", "s"), 2)
  ))

  # check that parameters are not equal (as in diagonal of distance matrix)
  if (all(params_x == params_y)) {
    return(0)
  }

  # calculate vector mu and Sigma for normal dist as in 5.2 of Heff & Tawn '04
  # TODO Have as arguments, or include Sig_fun for non-diagonal covariance
  mu_fun <- \(ce, y) {
    return((ce["a", ] * y) + ((y^ce["b", ]) * ce["m", ]))
  }
  Sig_fun <- \(ce, y) {
    return((y^ce["b", ]) * ce["s", ]) # assumed to be diagonal covariance
  }

  # loop across variables, corresponding to columns in params_x, sum JS
  mu_vec <- lapply(list(params_x, params_y), mu_fun, y = y_spec)
  Sig_vec <- lapply(list(params_x, params_y), Sig_fun, y = y_spec)

  # calculate JS divergence for each variable
  return(jsg_gauss(
    mu1 = mu_vec[[1]],
    mu2 = mu_vec[[2]],
    Sig1 = Sig_vec[[1]],
    Sig2 = Sig_vec[[2]],
    lambda = lambda,
    agg_fun = agg_fun,
    ... = ...
  ))
}

# TODO Replace scalar variance with covariance matrix!
#' @title Calculate skew-geometric Jensen-Shannon divergence between two
#' multivariate Gaussian distributions with diagonal covariance.
#' @description Function to calculate the skew-geometric Jensen-Shannon
#' divergence between two multivariate Gaussian distributions with diagonal
#' covariance.
#' This metric is symmetric, as required for clustering.
#' Formula from slide 18 of
#' https://franknielsen.github.io/M-JS/Slide-JSSymmetrization.pdf
#' @param mu1 Mean of the first Gaussian distribution.
#' @param mu2 Mean of the second Gaussian distribution.
#' @param Sig1 Variance of the first Gaussian distribution.
#' @param Sig2 Variance of the second Gaussian distribution.
#' @param lambda Skew parameter for skew-geometric JS divergence (default 0.5).
#' @param agg_fun Function to aggregate divergence values across variables
#' (default `mean`).
#' @return The skew-geometric Jensen-Shannon divergence between two
#' multivariate Gaussian distributions.
#' @keywords internal
jsg_gauss <- \(mu1, mu2, Sig1, Sig2, lambda = 0.5, agg_fun = mean, ...) {
  # Compute the effective variance
  Sig_a <- 1 / (((1 - lambda) * (1 / Sig1)) + (lambda * (1 / Sig2)))

  # Compute the weighted mean
  mu_a <- Sig_a *
    (((1 - lambda) * (1 / Sig1) * mu1) + (lambda * (1 / Sig2) * mu2))

  # Compute the divergence contribution for each dimension, according to
  # https://franknielsen.github.io/M-JS/Slide-JSSymmetrization.pdf (slide 18)
  # TODO Change to use squared terms duhh
  # TODO Could we speed this up?
  div_per_dim <- 0.5 * (
    ((1 - lambda) * mu1 * (1 / Sig1) * mu1) +
      (lambda * mu2 * (1 / Sig2) * mu2) -
      (mu_a * (1 / Sig_a) * mu_a) +
      log(((Sig1^(1 - lambda)) * (Sig2^lambda)) / Sig_a)
    # (1 - lambda) * log(Sig1) + lambda * log(Sig2) - log(Sig_a)
  )

  # # Sum per-dimension divs to get scalar
  # Take the MC average of the divergence contributions
  # TODO (optionally) return MC uncertainty?
  return(agg_fun(div_per_dim, ...))
}
