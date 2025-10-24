#' @title Fit marginal models for multivariate conditional extremes
#' @description Fit marginal models for multivariate conditional extremes.
#' @param data List of dataframes at each location containing data for each
#' variable. Optionally, can be an object of class `evc_marg` from a previous
#' call to `fit_ce`, with elements `marginal`, `data_thresh` and `original`.
#' @param mult_col Name of column specifying multivariate structure, default
#' "name".
#' @param vars Names of variable columns.
#' @param thresh_method Method to select marginal thresholds, one of
#' "value" (fixed value), "quantile" (quantile thresholding),
#' "regression" (covariate dependent thresholding via `qgam`) or "none", for no
#' thresholding (or where you have already thresholded the data prior to input).
#' @param thresh_args Arguments to be passed to thresholding function. For
#' "value", a numeric value (scalar for shared value or vector for each
#' `vars`) for simple value thresholding. For "quantile", a numeric value
#' (scalar for shared value or vector for each `vars`) for quantile
#' thresholding. For "regression", a list of arguments to be passed to
#' `qgam_thresh` (see documentation for details).
# not implemented from here
#' @param thresh_only If TRUE, only threshold data and do not fit marginal
#' models, Default: FALSE.
#' @param marg_method Method to fit marginal models, one of "ecdf" (empirical
#' CDF), "ismev" (using `ismev::gpd.fit`) or "evgam" (using
#' `evgam::evgam`).
#' @param marg_args Arguments to be passed to marginal fitting function.
#' @param ncores Number of cores to use for parallel computation, Default: 1.
#' @return Object of type `cecl_marg` for each location.
#' @rdname cecl_marg
#' @importFrom rlang .data :=
#' @export
# TODO Allow for any "name" column
cecl_marg <- \(
  data,
  mult_col = "name",
  vars = NULL,
  thresh_method = c("value", "quantile", "regression", "none"),
  thresh_args = NULL, # TODO Expand description
  thresh_only = FALSE,
  marg_method = c("ecdf", "ismev", "evgam"),
  marg_args = NULL,
  ncores = 1
) {
  ## Initial ##

  # initialise to remove `devtools::check()` note
  name <- NULL

  # Parallel setup
  apply_fun <- ifelse(ncores == 1, lapply, parallel::mclapply)
  ext_args <- NULL
  if (ncores > 1) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }

  ## Argument Checks ##

  thresh_method <- match.arg(thresh_method)
  marg_method <- match.arg(marg_method)

  # TODO Check arguments correctly specified for each thresh_method
  # Need to do regression!
  stopifnot(
    length(thresh_method) == 1 &&
      thresh_method %in% c("value", "quantile", "regression", "none")
  )

  if (thresh_method == "regression") {
    stopifnot(is.list(thresh_args))
    stopifnot("f" %in% names(thresh_args))
    stopifnot(all(names(thresh_args) %in% c("f", "qu", "jitter")))
  }

  # prepare data frame
  prep_out <- prep_df(data, mult_col, vars)
  data_df <- prep_out[[1]]
  vars <- prep_out[[2]] # in case vars = NULL in call

  # ensure that values or quantiles for thresholding are correct length
  if (thresh_method %in% c("value", "quantile")) {
    stopifnot(
      length(thresh_args) %in% c(1, length(vars))
    )
  }

  # threshold data using specified function
  thresh_out <- marg_thresh(
    data_df       = data_df,
    mult_col      = mult_col,
    vars          = vars,
    thresh_method = thresh_method,
    thresh_args   = thresh_args
  )
  data_thresh_out <- data_thresh <- thresh_out[[1]]
  locs_keep <- thresh_out[[2]] # locations with exceedances for all variables
  names(data_thresh_out) <- names(data_thresh) <- vars

  # Only keep locs with exceedances for all vars, otherwise can't do CE!
  data_df <- dplyr::filter(data_df, name %in% locs_keep)

  # original data split by location (for returning later)
  orig_dat <- data_df |>
    dplyr::group_by(.data[[mult_col]]) |>
    dplyr::group_split(.keep = TRUE)
  # add names to list elements
  names(orig_dat) <- purrr::map_chr(
    orig_dat,
    ~ as.character(.x[[mult_col]][1])
  )

  # return thresholded data only if specified
  # TODO Make into S3 class??
  if (thresh_only) {
    ret <- list(
      "data_thresh" = data_thresh,
      "original"    = orig_dat,
      "vars"        = vars,
      "call"        = match.call()
    )
    class(ret) <- c(
      "cecl_thresh",
      paste0("cecl_thresh_", thresh_method),
      class(ret)
    )

    return(ret)
  }

  # reverse splitting of data_thresh into locations
  data_thresh <- lapply(data_thresh, \(x) dplyr::bind_rows(unname(x)))

  # fit marginal models
  marginal <- fit_marg(
    data_df      = data_df, # TODO data_df will be weird if thresh = none
    data_thresh  = data_thresh,
    mult_col     = mult_col,
    vars         = vars,
    marg_method  = marg_method,
    marg_args    = marg_args,
    loop_fun     = loop_fun
  )
  if (marg_method == "evgam") {
    evgam_fit <- marginal$evgam_fit
    marginal <- marginal$marginal
  }
  names(marginal) <- locs_keep

  # Transform data to Laplace margins
  marginal_trans <- trans_marg(marginal, data_df, mult_col, vars)
  names(marginal_trans) <- locs_keep

  # return
  ret <- list(
    "marginal"    = marginal,
    "data_thresh" = data_thresh_out,
    "original"    = orig_dat,
    "transformed" = marginal_trans,
    "vars"        = vars,
    "call"        = match.call()
  )
  # remove marginal fits if ecdf method used
  if (marg_method == "ecdf") {
    ret$marginal <- NULL
  }
  # add evgam fit object if fitted
  if (exists("evgam_fit", envir = environment())) {
    names(evgam_fit) <- vars
    ret$evgam_fit <- evgam_fit
  }

  # make ret object of class `evc_marg`
  class(ret) <- c(
    "cecl_marg",
    paste0("cecl_marg_", marg_method),
    class(ret)
  )
  return(ret)
}

#' @title Prepare data frame for cecl_marg
#' @description Prepare data frame for cecl_marg: Convert list of data into
#' dataframe with `mult_col` column specifying multivariate structure.
#' @inheritParams cecl_marg
#' @return Data frame ready for `cecl_marg`.
#' @rdname prep_df
#' @keywords internal
prep_df <- \(data, mult_col = "name", vars) {
  # if already thresholded, just return data
  if (inherits(data, "cecl_thresh")) {
    return(list(data, data$vars))
  }

  nvars <- length(vars)
  # convert to data frame, if required
  if (!is.data.frame(data) && is.list(data)) {
    data_df <- dplyr::bind_rows(lapply(seq_along(data), \(i) {
      ret <- as.data.frame(data[[i]])
      # Add name column if not in list already
      if (!mult_col %in% names(ret)) {
        ret <- ret |>
          dplyr::mutate(!!mult_col := paste0("location_", i))
      }
    }))
    if (!is.null(vars)) {
      names(data_df)[seq_len(nvars)] <- vars
    } else {
      # Note: assuming all but last column are variables!
      vars <- paste("var", seq_len(ncol(data_df) - 1), sep = "")
      names(data_df)[seq_len(ncol(data_df) - 1)] <- vars
    }
  } else {
    data_df <- as.data.frame(data)
    # set vars if NULL
    if (is.null(vars)) {
      vars <- names(data_df)[!names(data_df) == mult_col]
    }
    # convert matrix names if required
    names(data_df)[names(data_df) %in% paste0("V", seq_len(nvars))] <- vars
  }

  # must have names column, and all of vars must be columns in data_df
  stopifnot("Must have a `mult_col` column" = "name" %in% names(data_df))
  stopifnot("All of `vars` must be in data" = all(vars %in% names(data_df)))

  # make name a factor in order of appearance
  data_df[[mult_col]] <- forcats::fct_inorder(data_df[[mult_col]])

  return(list(data_df, vars))
}

#' @title Threshold data for cecl_marg
#' @description Threshold data for cecl_marg using specified method.
#' @inheritParams cecl_marg
#' @return List of thresholded data frames for each variable.
#' @rdname marg_thresh
#' @keywords internal
marg_thresh <- \(
  data_df,
  mult_col = "name",
  vars,
  thresh_method,
  thresh_args
) {
  quantile <- excess <- NULL

  # locations (before thresholding)
  locs <- unique(data_df[[mult_col]])

  if (inherits(data_df, "cecl_thresh")) {
    return(list(data_df$data_thresh, locs))
  }

  if (thresh_method == "none") {
    data_thresh <- lapply(vars, \(x) {
      ret <- data_df |>
        # remove other responses, will be joined together after
        dplyr::select(-dplyr::all_of(vars[vars != x])) |>
        dplyr::mutate(
          thresh = NA_real_,
          excess = !!rlang::sym(x)
        ) |>
        # also split by location
        dplyr::group_split(.data[[mult_col]], .keep = TRUE)
      names(ret) <- purrr::map_chr(
        ret,
        ~ as.character(.x[[mult_col]][1])
      )
      ret
    })
    names(data_thresh) <- vars

    return(list(data_thresh, locs))
  }

  # If marg_val not specified, calculate thresh as quantile across all locs
  if (thresh_method %in% c("value", "quantile")) {
    if (thresh_method == "quantile") {
      marg_val <- mapply(
        quantile,
        data_df[, vars],
        thresh_args, # can be scalar or vector
        MoreArgs = list(na.rm = TRUE)
      )
    } else {
      marg_val <- thresh_args
      if (length(marg_val) == 1) {
        marg_val <- rep(marg_val, length(vars))
      }
    }
    names(marg_val) <- vars

    # for each variable, calculate excess over threshold
    data_thresh <- lapply(vars, \(x) {
      data_df |>
        # remove other responses, will be joined together after
        dplyr::select(-dplyr::all_of(vars[vars != x])) |>
        dplyr::mutate(
          thresh = marg_val[x],
          excess = !!rlang::sym(x) - marg_val[x]
        ) |>
        dplyr::filter(excess > 0) |>
        # also split by location
        identity()
    })
    # If thresh is a list, assume it is arguments to qgam_thresh
    # TODO: May be easier to just copy each vars column as response in data_df
    # Would allow for simpler formula specification
  }

  if (thresh_method == "regression") {
    data_thresh <- lapply(vars, \(x) {
      print(paste0("thresholding ", x))

      # Change formula to include response in question
      spec_params <- thresh_args
      spec_params$f <- lapply(thresh_args$f, \(f_spec) {
        stats::formula(stringr::str_replace_all(f_spec, "response", x))
      })
      # allow different thresholds for each variable
      if (length(spec_params$qu) > 1) {
        spec_params$qu <- spec_params$qu[vars == x]
      }

      # Run thresholding function for each response with specified args
      # TODO Also return evgam fits to ald
      quantile_fits <- do.call(
        qgam_thresh,
        args = c(
          list(data = data_df, response = x), # data args
          spec_params
        )
      )
      ret <- quantile_fits$data_thresh

      # return message where no exceedances are observed for any locations
      # TODO Implement for "normal" thresholding as well
      thresh_locs <- unique(ret[[mult_col]])
      if (length(thresh_locs) < length(locs)) {
        loc_missing <- setdiff(locs, thresh_locs)
        message(paste0(
          "No exceedances for variable ", x, " at: ",
          paste(loc_missing, collapse = ", "),
          ", removing for all variables"
        ))
      }
      ret
    })
  }

  # locations with exceedances for all variables
  locs_keep <- Reduce(
    intersect,
    lapply(data_thresh, \(x) unique(x[[mult_col]]))
  )

  # remove duplicate threshold rows kept through floating point errors
  data_thresh <- lapply(data_thresh, \(x) {
    ret <- dplyr::filter(x, .data[[mult_col]] %in% locs_keep) |>
      dplyr::group_by(
        .data[[mult_col]],
        dplyr::across(dplyr::any_of(c("date", !!vars)))
      ) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      # also split by location
      dplyr::group_split(.data[[mult_col]], .keep = TRUE)

    names(ret) <- purrr::map_chr(
      ret,
      ~ as.character(.x[[mult_col]][1])
    )
    ret
  })
  return(list(data_thresh, locs_keep))
}

#' @title `qgam` varying threshold
#' @description Fit varying threshold using quantile regression via `qgam`.
#' @param data Dataframe for one location which we wish to
#' threshold.
#' @param response Name of variable to threshold.
#' @param f Formula for `evgam` model, specified in `cecl_marg`.
#' @param qu Quantile to threshold at (see \link[qgam]{qgam} for details),
#' specified in `cecl_marg`, default 0.95.
#' @param jitter Add jitter to data to remove 0s, specified in `cecl_marg`,
#' default TRUE.
#' @param thresh return thresholded (i.e. filtered) data if TRUE, default TRUE.
#' @return Dataframe with `thresh` and `excess` columns, optionally thresholded.
#' @rdname qgam_thresh
#' @export
qgam_thresh <- \(
  data,
  response,
  f,
  qu = .95,
  jitter = TRUE,
  thresh = TRUE
) {
  excess <- NULL

  # jitter, if specified, to remove 0s when calculating quantiles
  # TODO Change jitter to match magnitude, may be too large for small data
  if (jitter == TRUE) {
    data <- data |>
      dplyr::mutate(dplyr::across(
        dplyr::all_of(response), ~ . + abs(
          stats::rnorm(dplyr::n(), 0, 1e-6)
        )
      ))
  }

  # fit the quantile regression model at qu'th percentile
  qgam_fit <- qgam::qgam(
    f,
    data,
    qu = qu
  )

  # add threshold to data
  predictors <- names(qgam_fit$var.summary)
  predictions <- data |>
    dplyr::mutate(thresh = qgam_fit$fitted.values) |>
    dplyr::distinct(dplyr::across(c(dplyr::all_of(predictors), thresh)))

  data_thresh <- data |>
    dplyr::left_join(predictions, by = predictors) |>
    dplyr::mutate(excess = !!rlang::sym(response) - thresh)
  # threshold if desired
  if (thresh == TRUE) {
    data_thresh <- dplyr::filter(data_thresh, excess > 0)
  }

  # return model fit and thresholded data
  return(list(
    "m"           = qgam_fit,
    "data_thresh" = data_thresh
  ))
}


#' @title Fit marginal models for cecl_marg
#' @description Fit marginal models for cecl_marg using specified method.
#' @inheritParams cecl_marg
#' @return List of marginal model fits for each location.
#' @rdname fit_marg
#' @keywords internal
fit_marg <- \(
  data_df,
  data_thresh,
  mult_col = "name",
  vars,
  marg_method,
  marg_args,
  loop_fun
) {
  thresh <- shape <- sigma <- xi <- NULL

  if (marg_method == "ecdf") {
    locs <- unique(data_df[[mult_col]])
    # dummy GPD fits
    marginal <- lapply(locs, \(loc) {
      gpd_y <- lapply(vars, \(var) {
        list(
          "sigma" = NA_real_,
          "xi" = NA_real_,
          # take as threshold the maximum value (to transform using only ECDF)
          "thresh" = data_thresh[[var]] |>
            dplyr::filter(.data[[mult_col]] == loc) |>
            dplyr::arrange(!!rlang::sym(var)) |>
            dplyr::slice(dplyr::n()) |>
            dplyr::pull(!!rlang::sym(var))
        )
      })
      names(gpd_y) <- vars
      gpd_y
    })
    names(marginal) <- locs
  }

  if (marg_method == "ismev") {
    # calculate for all locations
    marginal <- data_df |>
      dplyr::group_split(.data[[mult_col]], .keep = TRUE) |>
      loop_fun(\(x) {
        # pull marginal thresholds
        mth <- vapply(data_thresh, \(y) {
          y |>
            # need thresh for correct loc
            dplyr::filter(.data[[mult_col]] == x[[mult_col]][[1]]) |>
            dplyr::slice(1) |>
            dplyr::pull(thresh)
        }, numeric(1))

        # fit
        gpd_fits <- lapply(seq_along(vars), \(i) {
          fit <- ismev::gpd.fit(
            x[[vars[i]]],
            threshold = mth[i],
            show      = FALSE
          )
          ret <- list(
            "sigma"     = fit$mle[1],
            "xi"        = fit$mle[2],
            "thresh"    = fit$threshold[[1]],
            "name"      = x[[mult_col]][1]
          )
          names(ret)[4] <- mult_col
          ret
        })
        names(gpd_fits) <- vars
        gpd_fits
      })

    # add names (correctly!)
    names(marginal) <- purrr::map_chr(
      marginal,
      ~ as.character(.x[[1]][[mult_col]])
    )
  }

  # fit evgam model for each marginal
  if (marg_method == "evgam") {
    evgam_fit <- loop_fun(data_thresh, \(x) {
      fit_evgam(
        data      = x,
        pred_data = data_df,
        f         = marg_args$f
      )
    })

    # pull scale, shape and threshold for each variable
    marginal <- lapply(seq_along(evgam_fit), \(i) {
      # pull for each location (or predictor combo)
      # using signif to account for floating point errors
      params_df <- dplyr::select(
        evgam_fit[[i]]$predictions,
        !!mult_col,
        sigma = scale,
        xi = shape
      ) |>
        dplyr::distinct(signif(sigma, 6), signif(xi, 6), .keep_all = TRUE) |>
        # also pull in thresholds
        dplyr::left_join(
          dplyr::select(data_thresh[[i]], !!mult_col, thresh) |>
            dplyr::distinct(
              .data[[mult_col]],
              signif(thresh, 6),
              .keep_all = TRUE
            )
        ) |>
        dplyr::select(-dplyr::matches("signif"))

      # split into list by name as dependence pars will also be this way
      out <- params_df |>
        dplyr::relocate(!!mult_col, .after = dplyr::last_col()) |>
        dplyr::group_split(.data[[mult_col]], .keep = TRUE) |>
        lapply(as.vector, mode = "list")
      # pull names correctly from list objects
      names(out) <- purrr::map_chr(
        out,
        ~ as.character(.x[[mult_col]][1])
      )
      out
    })
    names(marginal) <- vars

    # transpose list from variables -> locations to locations -> variables
    marginal <- purrr::transpose(marginal)

    # also return evgam fits
    marginal <- list(
      "marginal" = marginal,
      "evgam_fit" = evgam_fit
    )
  }

  return(marginal)
}

#' @title Transform data to Laplace margins for cecl_marg
#' @description Transform data to Laplace margins for cecl_marg.
#' @param marginal List of fitted marginal models for each location.
#' @param data_df Data frame of original data.
#' @param mult_col Name of column specifying multivariate structure, default
#' "name".
#' @param vars Names of variable columns.
#' @return List of matrices transformed to Laplace margins for each location.
trans_marg <- \(
  marginal,
  data_df,
  mult_col = NULL,
  vars
) {
  # Calculate dependence from marginals (default output object)
  # first, transform margins to Laplace
  lapply(seq_along(marginal), \(i) {
    # semi-parametric CDF
    # TODO: More efficient to also split data_df by name and subset with i
    F_hat <- data_df |>
      dplyr::filter(.data[[mult_col]] == names(marginal)[i]) |>
      dplyr::select(dplyr::all_of(vars)) |>
      p_gpd_ecdf(marginal[[i]])
    # Laplace transform
    Y <- dlaplace(F_hat)
    colnames(Y) <- vars
    Y
  })
}

# TODO Move below to `utils.R` ??
#' @title Semi-parametric CDF for marginal models
#' @description Calculate semi-parametric CDF for marginal models:
#' empirical CDF below threshold, GPD above threshold.
#' @param dat Data matrix of observations.
#' @param gpd List of fitted GPD parameters for each variable.
#' @param n Number of observations.
#' @return Matrix of semi-parametric CDF values.
#' @rdname semi_par_cdf
#' @keywords internal
p_gpd_ecdf <- \(dat, gpd, n = nrow(dat)) {
  # As in Heff & Tawn '04, semiparametric mod uses ecdf below thresh, GPD above
  return(vapply(seq_along(gpd), \(i) {
    dat_spec <- dat[, i, drop = TRUE]
    spec_sigma <- gpd[[i]]$sigma
    spec_xi <- gpd[[i]]$xi
    spec_loc <- gpd[[i]]$thresh

    # TODO: Replace with ecdf fun from evc
    # order and sort data
    dat_spec_ord <- order(dat_spec)
    dat_spec_sort <- dat_spec[dat_spec_ord]

    # calculate ECDF
    m <- length(dat_spec)
    ecdf_vals <- (seq_len(m)) / (m + 1)
    # convert back to original order
    ecdf_dat_ord <- numeric(m)
    ecdf_dat_ord[dat_spec_ord] <- ecdf_vals

    # initialise
    cdf <- numeric(n)
    # ecdf (i.e. non-parametric) below threshold
    cdf[dat[, i] <= spec_loc] <- ecdf_dat_ord[dat_spec <= spec_loc]
    # GPD above threshold (parametric part of model)
    if (any(dat_spec > spec_loc)) {
      para <- pmax(
        0,
        1 + spec_xi * (dat_spec[dat_spec > spec_loc] - spec_loc) / spec_sigma
      )^(-1 / spec_xi)
      cdf[dat[, i] > spec_loc] <- 1 - (mean(dat_spec > spec_loc) * para)
    }
    cdf
  }, FUN.VALUE = numeric(n)))
}

#' @title Inverse semi-parametric CDF for marginal models
#' @description Calculate inverse semi-parametric CDF for marginal models:
#' empirical CDF below threshold, GPD above threshold.
#' @param F_hat Matrix of semi-parametric CDF values.
#' @param dat Data matrix of observations.
#' @param gpd List of fitted GPD parameters for each variable.
#' @return Matrix of reconstructed data values.
#' @rdname inv_semi_par_cdf
#' @keywords internal
d_gpd_ecdf <- \(F_hat, dat, gpd) {
  vapply(seq_along(gpd), \(i) {
    dat_spec <- dat[, i, drop = TRUE]
    stopifnot(names(gpd[[i]]) == c("sigma", "xi", "thresh"))

    spec_sigma <- gpd[[i]][[1]]
    spec_xi <- gpd[[i]][[2]]
    spec_loc <- gpd[[i]][[3]]

    n <- length(dat_spec)
    probs <- (1:n) / (n + 1) # Empirical CDF probabilities

    # Find closest probability match for each F_hat value
    px <- vapply(F_hat[, i], \(x, p) {
      p[[which.min(abs(x - p))]] # Nearest empirical CDF probability
    }, 0, p = probs)

    px <- as.integer(round(px * (1 + n)))
    res <- sort(dat_spec)[px] # Get corresponding data values

    # Adjust upper tail using GPD if above threshold
    i_F <- F_hat[, i] >= mean(dat_spec <= spec_loc) # Upper tail condition
    i_res <- res > spec_loc # Above threshold in reconstructed values
    i_adjust <- i_F & i_res # Both conditions met

    if (sum(i_adjust) > 0) {
      # Compute inverse GPD transformation
      p_above <- (1 - F_hat[i_adjust, i]) / mean(dat_spec > spec_loc)
      gpd_vals <- spec_loc + (spec_sigma / spec_xi) *
        ((pmax(0, p_above)^(-spec_xi)) - 1)

      # Order properly
      ordered_res <- res[i_adjust]
      order_idx <- order(ordered_res)
      ordered_res <- ordered_res[order_idx]
      ordered_res[
        length(ordered_res):(length(ordered_res) - length(gpd_vals) + 1)
      ] <- rev(sort(gpd_vals))
      ordered_res <- ordered_res[order(order_idx)]

      res[i_adjust] <- ordered_res
    }

    # Ensure final ordering matches input ordering
    res[order(F_hat[, i])] <- sort(res)

    res
  }, FUN.VALUE = numeric(nrow(F_hat)))
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
#' @keywords internal
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
#' @keywords internal
plaplace <- \(F_hat) {
  apply(to_matrix(F_hat), 2, \(x) {
    ifelse(x < 0, exp(x) / 2, 1 - exp(-x) / 2)
  })
}

#' @title Fit `evgam` model
#' @description Fit and generate predictions from `evgam` model.
#' @param data Dataframe for one location.
#' @param pred_data Dataframe for one location to predict on.
#' @param f Formula for `evgam` model.
#' @return List with model `m` and predictions `predictions`.
#' @rdname fit_evgam
#' @keywords internal
fit_evgam <- \(
  data,
  pred_data,
  mult_col = "name",
  # formula used in evgam, fitting to both scale and shape parameters
  f = list(
    excess ~ s(lon, lat), # increase smoothing on scale parameter
    ~ s(lon, lat) # shape parameter
  )
) {
  # ensure f is a formula
  f <- lapply(f, stats::formula)
  # fit evgam model
  m <- evgam::evgam(f, data = data, family = "gpd")

  # create predictions for unique rows in pred_data (ensures one pred per loc)
  predictors <- m$predictor.names
  pred_dat_distinct <- pred_data |>
    dplyr::distinct(
      .data[[mult_col]],
      dplyr::across(dplyr::all_of(predictors))
    )
  predictions <- cbind(
    pred_dat_distinct,
    stats::predict(m, pred_dat_distinct, type = "response")
  )

  # return model fit and predictions
  return(list(
    "m"           = m,
    "predictions" = predictions
  ))
}

#### Methods for cecl_marg class ####

#' @title `coef` method for `cecl_marg` class
#' @description Extract marginal model coefficients from `cecl_marg` object.
#' @param object Object of class `cecl_marg`.
#' @param ... Additional arguments (not used).
#' @return Data frame of marginal model coefficients.
#' @rdname coef.cecl_marg
#' @method coef cecl_marg
#' @export
# TODO What to do with varying threshold???
coef.cecl_marg <- \(object, ...) {
  stopifnot(inherits(object, "cecl_marg"))

  if (inherits(object, "cecl_marg_ecdf")) {
    stop(paste(
      "No marginal model coefficients to extract for 'ecdf' marg_method.",
      "Use 'thresh_only = TRUE' in 'cecl_marg' to only threshold data."
    ))
  }

  if (
    inherits(object, "cecl_marg_ismev") || inherits(object, "cecl_marg_evgam")
  ) {
    coefs <- lapply(object$marginal, \(loc) {
      do.call(rbind, lapply(loc, as.data.frame))
    })
    coefs_df <- do.call(rbind, lapply(names(coefs), \(loc_name) {
      loc_df <- coefs[[loc_name]]
      loc_df$name <- loc_name
      loc_df
    }))
    rownames(coefs_df) <- NULL
    coefs_df
  }

  stop("coef method not implemented for this marg_method")
}

#' @title `print` method for `cecl_marg` class
#' @description Print summary of `cecl_marg` object.
#' @param x Object of class `cecl_marg`.
#' @param ... Additional arguments (not used).
#' @return Printed summary of `cecl_marg` object.
#' @rdname print.cecl_marg
#' @method print cecl_marg
#' @export
print.cecl_marg <- \(x, ...) {
  stopifnot(inherits(x, "cecl_marg"))

  cat("Conditional Extremes Marginal Model Fit\n")
  cat("Number of sites:", length(x$original), "\n")
  cat("Variables:", paste(x$vars, collapse = ", "), "\n")

  if (inherits(x, "cecl_marg_ecdf")) {
    cat("Marginal method: Empirical CDF (no parametric fit)\n")
  } else if (inherits(x, "cecl_marg_ismev")) {
    cat("Marginal method: GPD fit via ismev::gpd.fit\n")
  } else if (inherits(x, "cecl_marg_evgam")) {
    cat("Marginal method: GPD fit via evgam::evgam\n")
  }
}

#' @title Summary method for `cecl_marg` objects
#' @description Print a summary of a `cecl_marg` object, including marginal
#' model coefficients.
#' @param object Object of class `cecl_marg`.
#' @param n Number of rows of coefficients to display (default all).
#' @param ... Additional arguments (not used).
#' @return Invisibly returns the input object after printing summary
#' information.
#' @method summary cecl_marg
#' @rdname summary.cecl_marg
#' @export
summary.cecl_marg <- \(object, n, ...) {
  stopifnot(inherits(object, "cecl_marg"))

  print(object)

  if (inherits(object, "cecl_marg_ecdf")) {
    stop(paste(
      "No marginal model coefficients to summarise for 'ecdf' marg_method.",
      "Use 'thresh_only = TRUE' in 'cecl_marg' to only threshold data."
    ))
  }

  coefs_df <- coef.cecl_marg(object)
  if (!missing(n)) {
    coefs_df <- utils::head(coefs_df, n)
    cat("\nMarginal Model Coefficients (first", n, "rows):\n")
  } else {
    cat("\nMarginal Model Coefficients:\n")
  }

  print(coefs_df)

  invisible(object)
}

#' @title Plot method for `cecl_marg` objects
#' @description Generate diagnostic plots for a `cecl_marg` object, including
#' QQ, PP, histogram, and return level plots.
#' @param x Object of class `cecl_marg`.
#' @param which Type of plot to generate. Choices are `"qq"` for QQ plot,
#' `"pp"` for PP plot, `"hist"` for histogram of residuals, or `"return"` for
#' return level plot.
#' @param loc Location name to plot.
#' @param mult_col Name of the column representing locations, default `"name"`.
#' @param var Variable name to plot.
#' @param ... Additional arguments (not used)
#' @return For QQ, PP, and histogram plots: invisible NULL after plotting.
#' For return level plots: a `ggplot` object.
#' @method plot cecl_marg
#' @rdname plot.cecl_marg
#' @export
plot.cecl_marg <- \(
  x, which = c("qq", "pp", "hist", "return"), loc, mult_col = "name", var, ...
) {
  stopifnot(inherits(x, "cecl_marg"))

  if (inherits(x, "cecl_marg_ecdf")) {
    stop(paste(
      "No residuals to plot for 'ecdf' marg_method.",
      "Use 'thresh_only = TRUE' in 'cecl_marg' to only threshold data."
    ))
  }

  which <- match.arg(which)

  if (missing(loc) || missing(var)) {
    stop("Please specify both 'loc' and 'var' to plot.")
  }

  if (!loc %in% names(x$original)) {
    stop(paste("Location", loc, "not found in the cecl_marg object."))
  }
  if (!var %in% x$vars) {
    stop(paste("Variable", var, "not found in the cecl_marg object."))
  }

  # Extract original and thresholded data for specified location and variable
  orig_data <- x$original[[loc]] |>
    dplyr::select(dplyr::all_of(var))
  thresh_data <- x$data_thresh[[var]][[loc]] |>
    dplyr::select(dplyr::all_of(var))

  # Calculate residuals based on marginal method
  if (inherits(x, "cecl_marg_ismev")) {
    gpd_params <- x$marginal[[loc]][[var]]
    residuals <- (thresh_data[[var]] - gpd_params$thresh) / gpd_params$sigma
  } else if (inherits(x, "cecl_marg_evgam")) {
    evgam_fit <- x$evgam_fit[[which(x$vars == var)]]
    pred_row <- evgam_fit$predictions |>
      dplyr::filter(.data[[x$mult_col]] == loc)
    sigma <- pred_row$scale
    residuals <- (thresh_data[[var]] - pred_row$thresh) / sigma
  } else {
    stop("Plot method not implemented for this marg_method")
  }

  # Generate specified plot
  if (which == "qq") {
    stats::qqplot(
      stats::qexp(stats::ppoints(length(residuals))),
      residuals,
      main = paste("QQ Plot for", var, "at", loc),
      xlab = "Theoretical Quantiles",
      ylab = "Sample Quantiles"
    )
    graphics::abline(0, 1, col = "red")
  } else if (which == "pp") {
    plot(
      stats::ppoints(length(residuals)),
      stats::pexp(sort(residuals)),
      main = paste("PP Plot for", var, "at", loc),
      xlab = "Theoretical Probabilities",
      ylab = "Sample Probabilities"
    )
    graphics::abline(0, 1, col = "red")
  } else if (which == "hist") {
    graphics::hist(
      residuals,
      breaks = 20,
      main = paste("Histogram of Residuals for", var, "at", loc),
      xlab = "Residuals"
    )
  } else if (which == "return") {
    #  Extract fitted parameters
    gpd_params <- x$marginal[[loc]][[var]]
    u <- gpd_params$thresh
    sigma <- gpd_params$sigma
    xi <- gpd_params$xi

    # Estimate exceedance rate
    n_total <- nrow(orig_data)
    n_exc <- nrow(thresh_data)
    lambda_u <- n_exc / n_total

    # Define return periods
    T_vals <- c(1.5, 2, 5, 10, 20, 50, 100, 200)
    z_T <- if (abs(xi) > 1e-6) {
      u + (sigma / xi) * ((T_vals * lambda_u)^xi - 1)
    } else {
      u + sigma * log(T_vals * lambda_u)
    }

    plot(
      T_vals, z_T,
      type = "b", pch = 19,
      main = paste("Return Level Plot for", var, "at", loc),
      xlab = "Return Period",
      ylab = "Return Level"
    )

    nboot <- 500
    zT_boot <- matrix(NA, nrow = nboot, ncol = length(T_vals))

    for (b in seq_len(nboot)) {
      sim_data <- rgpd(
        n     = nrow(thresh_data),
        u     = gpd_params$thresh,
        sigma = gpd_params$sigma,
        xi    = gpd_params$xi
      )

      fit_b <- ismev::gpd.fit(
        sim_data,
        threshold = gpd_params$thresh, show = FALSE
      )

      if (inherits(fit_b, "try-error")) next # skip failed fit

      sigma_b <- fit_b$mle[1]
      xi_b <- fit_b$mle[2]
      lambda_u_b <- lambda_u

      # skip if invalid MLEs
      if (any(is.na(fit_b$mle))) next
      xi_b <- ifelse(abs(xi_b) < 1e-6, 1e-6, xi_b)

      zT_boot[b, ] <- gpd_params$thresh +
        (sigma_b / xi_b) * ((T_vals * lambda_u_b)^xi_b - 1)
    }

    ci <- apply(zT_boot, 2, stats::quantile, probs = c(0.025, 0.975))
    z_T_lower <- ci[1, ]
    z_T_upper <- ci[2, ]
    graphics::lines(
      T_vals, z_T_upper,
      lty = 2, col = ggsci::pal_nejm()(n = 2)[2]
    )
    graphics::lines(
      T_vals, z_T_lower,
      lty = 2, col = ggsci::pal_nejm()(n = 2)[2]
    )
  }
}

#' @title CECL ggplot theme
#' @description Custom ggplot theme for CECL plots.
#' @param legend.position Position of legend in plot, default "bottom".
#' @param nejm_pal Logical indicating whether to use NEJM color palette,
#' default TRUE.
#' @return List of ggplot theme elements.
#' @rdname cecl_theme
#' @export
cecl_theme <- \(legend.position = "bottom", nejm_pal = TRUE) {
  ret <- ggplot2::theme_bw() + ggplot2::theme(
    legend.position = legend.position,
    plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
    axis.text = ggplot2::element_text(size = 12),
    axis.title = ggplot2::element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(size = 12),
    strip.text = ggplot2::element_text(size = 13, face = "bold"),
    strip.background = ggplot2::element_rect(fill = NA, colour = "black"),
    plot.tag = ggplot2::element_text(size = 16, face = "bold"),
    panel.background = ggplot2::element_rect(fill = NA, colour = "black")
  )
  ret <- list(ret)
  if (nejm_pal == TRUE) {
    ret <- c(ret, list(ggsci::scale_colour_nejm(), ggsci::scale_fill_nejm()))
  }
  return(ret)
}

#' @title `ggplot` method for `cecl_marg` class
#' @description Generate ggplot diagnostic plots for `cecl_marg` object.
#' @param data object of class `cecl_marg`.
#' @param mapping Not used.
#' @param ... Additional arguments (not used).
#' @param environment Parent frame environment.
#' @inheritParams plot.cecl_marg
#' @return ggplot diagnostic plot for `cecl_marg` object.
#' @rdname ggplot.cecl_marg
#' @method ggplot cecl_marg
#' @export
#' @importFrom ggplot2 ggplot
ggplot.cecl_marg <- \(
  data = NULL,
  mapping = ggplot2::aes(),
  which = c("qq", "pp", "hist", "return"),
  loc,
  mult_col = "name",
  var,
  ...,
  environment = parent.frame()
) {
  stopifnot(inherits(data, "cecl_marg"))

  if (inherits(data, "cecl_marg_ecdf")) {
    stop(paste(
      "No residuals to plot for 'ecdf' marg_method.",
      "Use 'thresh_only = TRUE' in 'cecl_marg' to only threshold data."
    ))
  }

  which <- match.arg(which)

  if (missing(loc) || missing(var)) {
    stop("Please specify both 'loc' and 'var' to plot.")
  }

  if (!loc %in% names(data$original)) {
    stop(paste("Location", loc, "not found in the cecl_marg object."))
  }
  if (!var %in% data$vars) {
    stop(paste("Variable", var, "not found in the cecl_marg object."))
  }

  quantile <- lower <- upper <- NULL

  # Extract original and thresholded data for specified location and variable
  orig_data <- data$original[[loc]] |>
    dplyr::select(dplyr::all_of(var))
  thresh_data <- data$data_thresh[[var]][[loc]] |>
    dplyr::select(dplyr::all_of(var))

  # Calculate residuals based on marginal method
  if (inherits(data, "cecl_marg_ismev")) {
    gpd_params <- data$marginal[[loc]][[var]]
    residuals <- (thresh_data[[var]] - gpd_params$thresh) / gpd_params$sigma
  } else if (inherits(data, "cecl_marg_evgam")) {
    evgam_fit <- data$evgam_fit[[which(data$vars == var)]]
    pred_row <- evgam_fit$predictions |>
      dplyr::filter(.data[[data$mult_col]] == loc)
    sigma <- pred_row$scale
    residuals <- (thresh_data[[var]] - pred_row$thresh) / sigma
  } else {
    stop("ggplot method not implemented for this marg_method")
  }

  res_df <- data.frame(residuals = residuals)

  # Generate specified ggplot
  if (which == "qq") {
    ggplot2::ggplot(res_df, ggplot2::aes(
      sample = residuals,
    )) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line(col = "red") +
      ggplot2::labs(
        title = paste("QQ Plot for", var, "at", loc),
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
      ) +
      cecl_theme()
  } else if (which == "pp") {
    ggplot2::ggplot(res_df, ggplot2::aes(
      x = stats::ppoints(length(residuals)),
      y = stats::pexp(sort(residuals))
    )) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(slope = 1, intercept = 0, col = "red") +
      ggplot2::labs(
        title = paste("PP Plot for", var, "at", loc),
        x = "Theoretical Probabilities",
        y = "Sample Probabilities"
      ) +
      cecl_theme()
  } else if (which == "hist") {
    ggplot2::ggplot(res_df, ggplot2::aes(x = residuals)) +
      ggplot2::geom_histogram(
        bins = 20,
        fill = ggsci::pal_nejm()(1)[1],
        color = "black"
      ) +
      ggplot2::labs(
        title = paste("Histogram of Residuals for", var, "at", loc),
        x = "Residuals"
      ) +
      cecl_theme()
  } else if (which == "return") {
    # Extract parameters
    gpd_params <- data$marginal[[loc]][[var]]
    u <- gpd_params$thresh
    sigma <- gpd_params$sigma
    xi <- gpd_params$xi

    # Exceedance rate
    n_total <- nrow(orig_data)
    n_exc <- nrow(thresh_data)
    lambda_u <- n_exc / n_total

    # Return periods
    T_vals <- c(1.5, 2, 5, 10, 20, 50, 100, 200)
    z_T <- if (abs(xi) > 1e-6) {
      u + (sigma / xi) * ((T_vals * lambda_u)^xi - 1)
    } else {
      u + sigma * log(T_vals * lambda_u)
    }

    df <- data.frame(
      T_val = T_vals,
      z_T = z_T
    )

    # bootstrap CI
    nboot <- 200
    zT_list <- list()
    for (b in seq_len(nboot)) {
      sim_data <- rgpd(
        n     = nrow(thresh_data),
        u     = u,
        sigma = sigma,
        xi    = xi
      )
      fit_b <- try(
        ismev::gpd.fit(sim_data, threshold = u, show = FALSE),
        silent = TRUE
      )
      if (inherits(fit_b, "try-error") || any(is.na(fit_b$mle))) next
      sigma_b <- fit_b$mle[1]
      xi_b <- ifelse(abs(fit_b$mle[2]) < 1e-6, 1e-6, fit_b$mle[2])
      zT_list[[length(zT_list) + 1]] <- u +
        (sigma_b / xi_b) * ((T_vals * lambda_u)^xi_b - 1)
    }
    if (length(zT_list) > 0) {
      zT_boot <- do.call(rbind, zT_list)
      df$lower <- apply(zT_boot, 2, stats::quantile, probs = 0.025)
      df$upper <- apply(zT_boot, 2, stats::quantile, probs = 0.975)
    }

    # ggplot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = T_vals, y = z_T)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::labs(
        title = paste("Return Level Plot for", var, "at", loc),
        x = "Return Period",
        y = "Return Level"
      ) +
      cecl_theme()

    # Add CI ribbon if available
    if ("lower" %in% names(df) && "upper" %in% names(df)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        alpha = 0.2,
        fill = ggsci::pal_nejm()(1)
      )
    }

    return(p)
  }
}
