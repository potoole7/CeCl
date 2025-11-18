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
#' "qgam" (covariate dependent thresholding via `qgam`) or "none", for no
#' thresholding (or where you have already thresholded the data prior to input).
#' @param thresh_args Arguments to be passed to thresholding function. For
#' "value", a numeric value (scalar for shared value or vector for each
#' `vars`) for simple value thresholding. For "quantile", a numeric value
#' (scalar for shared value or vector for each `vars`) for quantile
#' thresholding. For "qgam", a list of arguments to be passed to
#' `qgam_thresh` (see documentation for details).
# not implemented from here
#' @param thresh_only If TRUE, only threshold data and do not fit marginal
#' models, Default: FALSE.
#' @param marg_method Method to fit marginal models, one of "ecdf" (empirical
#' CDF), "ismev" (using `ismev::gpd.fit`) or "evgam" (using
#' `evgam::evgam`).
#' @param marg_args Arguments to be passed to marginal fitting function.
#' @param ret_obj If TRUE and `thresh_method` is "qgam" or
#' `marg_method` is "ismev" or "evgam", return fitted
#' marginal models. Default: TRUE.
#' @param ncores Number of cores to use for parallel computation, Default: 1.
#' @return Object of type `cecl_marg` for each location.
#' @examples
#' # simulate some data
#' set.seed(123)
#' n_locs <- 5
#' n <- 1000
#' vars <- c("X1", "X2")
#' df <- do.call(
#'   rbind,
#'   lapply(1:n_locs, function(i) {
#'     data.frame(
#'       X1 = rgpd(n, u = 0, sigma = 1, xi = 0.2),
#'       X2 = rgpd(n, u = 0, sigma = 2, xi = -0.1),
#'       name = paste0("loc_", i)
#'     )
#'   })
#' )
#' # fit marginal models with quantile thresholding and ISMEV GPD fits
#' marg_fit <- cecl_marg(
#'   df,
#'   thresh_method = "quantile",
#'   thresh_args = 0.9,
#'   marg_method = "ismev",
#'   ncores = 1
#' )
#' # inspect output
#' marg_fit
#' summary(marg_fit)
#'
#' # fit marginal models using evgam
#' marg_fit_evgam <- cecl_marg(
#'   df,
#'   thresh_method = "quantile",
#'   thresh_args = 0.9,
#'   marg_method = "evgam",
#'   marg_args = list(f = list("excess ~ name", "~ name")),
#'   ncores = 1
#' )
#'
#' marg_fit_evgam
#' @rdname cecl_marg
#' @importFrom rlang .data :=
#' @export
cecl_marg <- \(
  data,
  mult_col = "name",
  vars = NULL,
  thresh_method = c("value", "quantile", "qgam", "none"),
  thresh_args = NULL, # TODO Expand description
  thresh_only = FALSE,
  marg_method = c("ecdf", "ismev", "evgam"),
  marg_args = NULL,
  ret_obj = TRUE,
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

  # data type
  stopifnot(is.list(data) || inherits(data, "evc_marg"))

  # mult_col and vars
  stopifnot(is.character(mult_col), length(mult_col) == 1)
  if (!is.null(vars)) stopifnot(is.character(vars))

  # thresholding method
  if (thresh_method %in% c("value", "quantile")) {
    stopifnot(!is.null(thresh_args))
    stopifnot(is.numeric(thresh_args))
    if (!is.null(vars) && !(length(thresh_args) %in% c(1, length(vars)))) {
      stop("Length of 'thresh_args' must be 1 or equal to number of variables.")
    }
  } else if (thresh_method == "qgam") {
    stopifnot(is.list(thresh_args))
    stopifnot("f" %in% names(thresh_args))
    stopifnot(
      "Currently only support the same `qgam` formula for all variables." =
        length(thresh_args$f) >= 1
    )
  }

  # marginal method
  stopifnot(marg_method %in% c("ecdf", "ismev", "evgam"))

  # thresh_only
  stopifnot(is.logical(thresh_only), length(thresh_only) == 1)

  # ncores
  stopifnot(is.numeric(ncores), length(ncores) == 1, ncores >= 1)

  # TODO Check arguments correctly specified for each thresh_method
  # Need to do qgam!
  stopifnot(
    length(thresh_method) == 1 &&
      thresh_method %in% c("value", "quantile", "qgam", "none")
  )

  if (thresh_method == "qgam") {
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
    thresh_args   = thresh_args,
    ret_obj       = ret_obj
  )
  data_thresh_out <- data_thresh <- thresh_out[[1]]
  locs_keep <- thresh_out[[2]] # locations with exceedances for all variables
  if (thresh_method == "qgam" && ret_obj) {
    qgam_fit <- thresh_out[[3]]
  }
  names(data_thresh_out) <- names(data_thresh) <- vars

  # Only keep locs with exceedances for all vars, otherwise can't do CE!
  data_df <- dplyr::filter(data_df, name %in% locs_keep)
  if (length(locs_keep) == 0) {
    stop("No locations have exceedances for all variables after thresholding.")
  }

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
    if (thresh_method == "qgam" && ret_obj) {
      ret$qgam_fit <- qgam_fit
    }
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
    loop_fun     = loop_fun,
    ret_obj      = ret_obj
  )
  if (marg_method == "evgam" && ret_obj) {
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
  # add qgam_fit & evgam fit object if fitted
  if (exists("qgam_fit", envir = environment())) {
    ret$qgam_fit <- qgam_fit
  }
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
  thresh_args,
  ret_obj = TRUE
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
        dplyr::filter(excess > 0)
    })
    # If thresh is a list, assume it is arguments to qGam_thresh
    # TODO: May be easier to just copy each vars column as response in data_df
    # Would allow for simpler formula specification
  }

  if (thresh_method == "qgam") {
    data_thresh <- lapply(vars, \(x) {
      print(paste0("thresholding ", x))

      # add thresh_args and ret_obj
      spec_params <- c(thresh_args)
      # Change formula to include response in question
      spec_params$f <- stats::formula(
        stringr::str_replace_all(spec_params$f, "response", x)
      )
      # allow different thresholds for each variable
      if (length(spec_params$qu) > 1) {
        spec_params$qu <- spec_params$qu[vars == x]
      }

      # Run thresholding function for each response with specified args
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
      if (ret_obj) {
        ret <- list(
          "data_thresh" = ret,
          "qgam_fit"    = quantile_fits$m
        )
      }
      ret
    })

    if (ret_obj) {
      qgam_fit <- lapply(data_thresh, \(x) x$qgam_fit)
      names(qgam_fit) <- vars
      data_thresh <- lapply(data_thresh, \(x) x$data_thresh)
      names(data_thresh) <- vars
    }
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

  ret <- list("data_thresh" = data_thresh, "locs_keep" = locs_keep)
  if (thresh_method == "qgam" && ret_obj) {
    ret$qgam_fit <- qgam_fit
  }
  class(ret) <- c(
    "cecl_thresh",
    paste0("cecl_thresh_", thresh_method),
    class(ret)
  )
  return(ret)
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
#' @keywords internal
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

  # return thresholded data (and model fit if desired)
  ret <- list(
    "m" = qgam_fit,
    "data_thresh" = data_thresh
  )

  return(ret)
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
  loop_fun,
  ret_obj = TRUE
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
            "fit"       = fit,
            "sigma"     = fit$mle[1],
            "xi"        = fit$mle[2],
            "thresh"    = fit$threshold[[1]],
            "name"      = x[[mult_col]][1]
          )
          names(ret)[5] <- mult_col
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

    if (ret_obj == FALSE) {
      marginal <- lapply(marginal, \(x) {
        lapply(x, `[`, c("sigma", "xi", "thresh", "name"))
      })
    }
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
            ),
          by = mult_col # TODO Should this also include predictor columns?
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
    if (ret_obj == FALSE) {
      marginal <- marginal$marginal
    }
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
#' @export
trans_marg <- \(
  marginal,
  data_df,
  mult_col = "name",
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
      ret <- do.call(rbind, lapply(loc, \(var) {
        as.data.frame(var[names(var) != "fit"])
      }))
      ret$var <- names(loc)
      ret
    })
    coefs_df <- do.call(rbind, lapply(names(coefs), \(loc_name) {
      loc_df <- coefs[[loc_name]]
      loc_df$name <- loc_name
      loc_df
    }))
    rownames(coefs_df) <- NULL
    coefs_df <- coefs_df[, c("name", "var", "thresh", "sigma", "xi")]
    return(coefs_df)
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

  cat("Call:\n")
  print(x$call)
  cat("\nNumber of locations:", length(x$transformed), "\n")
  cat("Variables:", paste(x$vars, collapse = ", "), "\n")
  cat(
    "Marginal method:",
    sub("cecl_marg_", "", class(x)[2]),
    "\n"
  )

  invisible(x)
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

#' @title Bootstrap QQ/PP envelope for GPD exceedances
#' @description Generate bootstrap QQ/PP envelopes for GPD exceedances.
#' @param exceedances Numeric vector of exceedances (data - threshold).
#' @param gpd_params List with elements `thresh`, `sigma`, `xi`
#' specifying fitted GPD parameters.
#' @param type Type of envelope to generate: "qq" for QQ plot,
#' "pp" for PP plot, or both.
#' @param nboot Number of bootstrap simulations, set to NULL for no
#' envelopes, Default 200.
#' @param refit Logical; if TRUE, refit GPD model for each simulated dataset
#' (slower, includes parameter uncertainty). If FALSE, use fixed fitted
#' parameters.
#' @param ci_quantiles Quantiles for confidence intervals,
#' Default c(0.025, 0.975).
#'
bootstrap_pp_qq <- \(
  exceedances,
  gpd_params,
  type = c("qq", "pp"),
  nboot = 200,
  refit = FALSE,
  ci_quantiles = c(0.025, 0.975)
) {
  type <- match.arg(type, several.ok = TRUE)
  if (!is.numeric(exceedances) || length(exceedances) < 2) {
    stop("exceedances must be numeric length >= 2")
  }

  # fitted params
  u <- gpd_params$thresh
  sigma <- gpd_params$sigma
  xi <- gpd_params$xi

  n_exc <- length(exceedances)
  p <- stats::ppoints(n_exc)

  # observed quantities
  observed_sample_q <- sort(exceedances)
  # PITs under fitted model
  observed_pit <- pgpd(exceedances, u = 0, sigma = sigma, xi = xi)
  observed_pp_sorted <- sort(observed_pit)
  # theoretical quantiles for QQ x-axis
  theor_q <- qgpd(p, u = 0, sigma = sigma, xi = xi)

  # preallocate matrix to hold simulated quantiles/PITs
  sim_mat_qq <- NULL
  sim_mat_pp <- NULL
  ok_qq <- logical(nboot)
  ok_pp <- logical(nboot)
  if ("qq" %in% type) {
    sim_mat_qq <- matrix(NA_real_, nrow = nboot, ncol = n_exc)
  }
  if ("pp" %in% type) {
    sim_mat_pp <- matrix(NA_real_, nrow = nboot, ncol = n_exc)
  }

  # if nboot = 0, return NA envelopes
  if (nboot == 0) {
    res <- list(nboot = nboot, refit = refit, ci_quantiles = ci_quantiles)
    if ("qq" %in% type) {
      qq_df <- data.frame(
        p = p, theor = theor_q, sample = observed_sample_q,
        lower = NA_real_, upper = NA_real_
      )
      res$n_succ_qq <- 0L
      res$qq <- qq_df
    }
    if ("pp" %in% type) {
      pp_df <- data.frame(
        p = p, model = observed_pp_sorted,
        lower = NA_real_, upper = NA_real_
      )
      res$n_succ_pp <- 0L
      res$pp <- pp_df
    }
    return(res)
  }
  for (b in seq_len(nboot)) {
    # simulate n_exc raw values from rgpd with threshold = u
    sim_raw <- try(
      rgpd(n = n_exc, u = u, sigma = sigma, xi = xi),
      silent = TRUE
    )
    if (inherits(sim_raw, "try-error") || any(is.na(sim_raw))) {
      next
    }

    # convert to exceedances (Y = X - u)
    sim_exc <- sim_raw - u

    if (!refit) {
      # fixed-parameter envelope: evaluate simulated sample quantiles or PITs
      # under fitted params
      if ("qq" %in% type) {
        sim_mat_qq[b, ] <- sort(sim_exc) # sample quantiles of simulated data
        ok_qq[b] <- TRUE
      }
      if ("pp" %in% type) {
        sim_pit <- pgpd(sim_exc, u = 0, sigma = sigma, xi = xi)
        sim_mat_pp[b, ] <- sort(sim_pit)
        ok_pp[b] <- TRUE
      }
    } else {
      # refit to simulated raw data (threshold = u)
      fit_b <- try(
        ismev::gpd.fit(sim_raw, threshold = u, show = FALSE),
        silent = TRUE
      )
      if (inherits(fit_b, "try-error") || any(is.na(fit_b$mle))) {
        next
      }

      sigma_b <- fit_b$mle[1]
      xi_b <- fit_b$mle[2]

      # compute same quantities but using the fitted parameters where needed
      if ("qq" %in% type) {
        # Use sample quantiles of simulated exceedances
        sim_mat_qq[b, ] <- sort(sim_exc)
        ok_qq[b] <- TRUE
      }
      if ("pp" %in% type) {
        # PIT under the refit params
        sim_pit_b <- pgpd(sim_exc, u = 0, sigma = sigma_b, xi = xi_b)
        sim_mat_pp[b, ] <- sort(sim_pit_b)
        ok_pp[b] <- TRUE
      }
    }
  }

  # prepare result object
  res <- list(nboot = nboot, refit = refit, ci_quantiles = ci_quantiles)

  # QQ results (if requested)
  ci <- ci_quantiles
  if ("qq" %in% type) {
    if (!any(ok_qq)) {
      message("bootstrap_pp_qq: no successful QQ simulations")
      qq_df <- data.frame(
        p = p, theor = theor_q, sample = observed_sample_q,
        lower = NA_real_, upper = NA_real_
      )
      res$n_succ_qq <- 0L
    } else {
      sim_mat_qq_ok <- sim_mat_qq[ok_qq, , drop = FALSE]
      n_succ_qq <- nrow(sim_mat_qq_ok)
      lower_qq <- apply(
        sim_mat_qq_ok, 2, stats::quantile,
        probs = ci[1], na.rm = TRUE
      )
      upper_qq <- apply(
        sim_mat_qq_ok, 2, stats::quantile,
        probs = ci[2], na.rm = TRUE
      )
      qq_df <- data.frame(
        p = p, theor = theor_q, sample = observed_sample_q,
        lower = lower_qq, upper = upper_qq
      )
      res$n_succ_qq <- n_succ_qq
    }
    res$qq <- qq_df
  }

  # PP results (if requested)
  if ("pp" %in% type) {
    if (!any(ok_pp)) {
      message("bootstrap_pp_qq: no successful PP simulations")
      pp_df <- data.frame(
        p = p, model = observed_pp_sorted,
        lower = NA_real_, upper = NA_real_
      )
      res$n_succ_pp <- 0L
    } else {
      sim_mat_pp_ok <- sim_mat_pp[ok_pp, , drop = FALSE]
      n_succ_pp <- nrow(sim_mat_pp_ok)
      lower_pp <- apply(
        sim_mat_pp_ok, 2, stats::quantile,
        probs = ci[1], na.rm = TRUE
      )
      upper_pp <- apply(
        sim_mat_pp_ok, 2, stats::quantile,
        probs = ci[2], na.rm = TRUE
      )
      pp_df <- data.frame(
        p = p, model = observed_pp_sorted,
        lower = lower_pp, upper = upper_pp
      )
      res$n_succ_pp <- n_succ_pp
    }
    res$pp <- pp_df
  }

  return(res)
}

#' @title Plot method for `cecl_marg` objects
#' @description Generate diagnostic plots for a `cecl_marg` object, including
#' QQ, PP, histogram, and return level plots.
#' @param x Object of class `cecl_marg`.
#' @param which Type of plot to generate. Choices are `"qq"` for QQ plot,
#' `"pp"` for PP plot, `"hist"` for histogram of residuals, `"return"` for
#' return level plot, and `"transformed"` for a plot of the Laplace-scale
#' transformed data. Several plots can also be produced.
#' @param loc Location name to plot.
#' @param mult_col Name of the column representing locations, default `"name"`.
#' @param var Variable name to plot.
#' @param cond_var Optional conditioning variable for plotting Laplace
#' transformed data, Default NULL.
#' @param plot_dens Logical indicating whether to plot histogram on top of
#' density estimate, Default TRUE.
#' @param return_periods Return periods for return level plot,
#' Default c(1.5, 2.5, 5, 10, 20, 50, 100, 200).
#' @param nboot Number of bootstrap samples for confidence
#' intervals, set to NULL for none confidence interval, Default: 200.
#' @param refit Logical indicating whether to refit GPD for each bootstrap
#' sample when calculating uncertainty envelopes, Default FALSE.
#' @param ci_quantiles Quantiles for confidence intervals,
#' Default c(0.025, 0.975).
#' @param log_scale Logical indicating whether to use log scale for return
#' level plot x-axis, Default TRUE.
#' @param ... Additional arguments to pass to plotting.
#' @method plot cecl_marg
#' @rdname plot.cecl_marg
#' @export
#' @examples
#' # simulate some data
#' set.seed(123)
#' n_locs <- 5
#' n <- 1000
#' vars <- c("X1", "X2")
#' df <- do.call(
#'   rbind,
#'   lapply(1:n_locs, function(i) {
#'     data.frame(
#'       X1 = rgpd(n, u = 0, sigma = 1, xi = 0.2),
#'       X2 = rgpd(n, u = 0, sigma = 2, xi = -0.1),
#'       name = paste0("loc_", i)
#'     )
#'   })
#' )
#' # fit marginal models with quantile thresholding and ISMEV GPD fits
#' marg_fit <- cecl_marg(
#'   df,
#'   thresh_method = "quantile",
#'   thresh_args = 0.9,
#'   marg_method = "ismev",
#'   ncores = 1
#' )
#' # plot for a specific location and variable
#' plot(marg_fit, which = "qq", loc = "loc_1", var = "X1")
plot.cecl_marg <- \(
  x,
  which = c("qq", "pp", "hist", "return", "transformed"),
  loc,
  mult_col = "name",
  var,
  cond_var = NULL,
  plot_dens = TRUE,
  return_periods = c(1.5, 2.5, 5, 10, 20, 50, 100, 200),
  nboot = 200,
  refit = FALSE,
  ci_quantiles = c(0.025, 0.975),
  log_scale = TRUE,
  ...
) {
  stopifnot(inherits(x, "cecl_marg"))

  if (inherits(x, "cecl_marg_ecdf") && which != "transformed") {
    stop(paste(
      "No residuals to plot for 'ecdf' `marg_method`,",
      "only 'transformed' plot available."
    ))
  }

  which <- match.arg(which, several.ok = TRUE)

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
  resid_fun <- \(q, gpd_params) do.call(
    pgpd,
    c(
      list(q = q),
      c(stats::setNames(
        gpd_params[c("sigma", "xi")],
        c("sigma", "xi")
      ), "u" = 0)
    )
  )

  # Calculate residuals based on marginal method
  if (any(which != "transformed")) {
    gpd_params <- x$marginal[[loc]][[var]]
    if (inherits(x, "cecl_marg_ismev")) {
      sigma <- gpd_params$sigma
      xi <- gpd_params$xi
      exceedances <- thresh_data[[var]] - gpd_params$thresh
      residuals <- resid_fun(exceedances, gpd_params)
    } else if (inherits(x, "cecl_marg_evgam")) {
      # stop("Plot method not implemented for evgam marg_method yet.")
      evgam_fit <- x$evgam_fit[[which(x$vars == var)]]
      pred_row <- evgam_fit$predictions |>
        dplyr::filter(.data[[mult_col]] == loc)
      sigma <- pred_row$scale
      xi <- pred_row$shape
      # exceedances <- thresh_data[[var]] - pred_row$thresh
      exceedances <- thresh_data[[var]] - gpd_params$thresh
      residuals <- resid_fun(exceedances, list(sigma = sigma, xi = xi))
    } else {
      stop("Plot method not implemented for this marg_method")
    }
    if (!is.null(cond_var)) {
      message("Ignoring `cond_var` for which != 'transformed'")
    }
  }
  # For transformed plot
  if (any(which == "transformed")) {
    stopifnot("must specify `cond_var`" = !is.null(cond_var))
    stopifnot("`cond_var` must differ from `var`." = cond_var != var)
    stopifnot(
      "`cond_var` not found in cecl_marg object." = cond_var %in% x$vars
    )
    transformed <- x$transformed[[loc]][, c(var, cond_var), drop = FALSE]
  }

  # calculate uncertainty for PP and/or QQ plots
  if (any(which %in% c("pp", "qq"))) {
    env <- bootstrap_pp_qq(
      exceedances = exceedances,
      gpd_params = gpd_params,
      type = intersect(which, c("pp", "qq")),
      nboot = nboot,
      refit = refit,
      ci_quantiles = ci_quantiles
    )
  }

  # Ask if multiple plots needed
  mf <- graphics::par("mfrow")
  capacity <- mf[1] * mf[2]
  # save original par settings & ensure reset
  op <- graphics::par(ask = length(which) > capacity)
  on.exit(graphics::par(op))

  # Generate specified plots
  for (w in which) {
    if (w == "transformed") {
      plot(
        transformed[, cond_var],
        transformed[, var],
        xlab = cond_var,
        ylab = var,
        ...
      )
    } else if (w == "qq") {
      qq_df <- env$qq
      stats::qqplot(
        qq_df$theor,
        qq_df$sample,
        xlab = "Theoretical",
        ylab = "Sample",
        ...
      )
      # add uncertainty
      if (nboot > 0) {
        graphics::polygon(
          c(qq_df$theor, rev(qq_df$theor)),
          c(qq_df$lower, rev(qq_df$upper)),
          col = grDevices::rgb(0.7, 0.7, 0.7, 0.4),
          border = NA
        )
      }
      graphics::points(qq_df$theor, qq_df$sample)
      graphics::abline(0, 1, col = "red")
    } else if (w == "pp") {
      pp_df <- env$pp
      plot(
        pp_df$p,
        pp_df$model,
        xlab = "Theoretical",
        ylab = "Model"
      )
      if (nboot > 0) {
        graphics::polygon(
          c(pp_df$p, rev(pp_df$p)),
          c(pp_df$lower, rev(pp_df$upper)),
          col = grDevices::rgb(0.7, 0.7, 0.7, 0.4),
          border = NA
        )
      }
      graphics::points(pp_df$p, pp_df$model)
      graphics::abline(0, 1, col = "red")
    } else if (w == "hist") {
      # if main not specified, set to NULL (hist automatically adds title)
      plot_args <- list(
        x      = residuals,
        breaks = 20,
        xlab   = "Residuals",
        prob   = ifelse(plot_dens, TRUE, FALSE),
        ...
      )
      if (!"main" %in% names(plot_args)) {
        plot_args[["main"]] <- list(NULL)
      }
      do.call(graphics::hist, plot_args)
      # add density plot if specified
      if (plot_dens) {
        graphics::lines(stats::density(residuals), col = "red", lwd = 2)
      }
    } else if (w == "return") {
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
      T_vals <- return_periods
      z_T <- if (abs(xi) > 1e-6) {
        u + (sigma / xi) * ((T_vals * lambda_u)^xi - 1)
      } else {
        u + sigma * log(T_vals * lambda_u)
      }

      if (log_scale) {
        T_vals_plot <- log(T_vals)
        xlab <- "Return Period (log)"
      } else {
        T_vals_plot <- T_vals
        xlab <- "Return Period"
      }

      plot(
        T_vals_plot, z_T,
        type = "b", pch = 19,
        xlab = xlab,
        ylab = "Return Level",
        ...
      )

      zT_boot <- matrix(NA, nrow = nboot, ncol = length(T_vals))

      # return NULL and exit if uncertainty not requested
      if (nboot == 0) {
        next
      }

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

      ci <- apply(zT_boot, 2, stats::quantile, probs = ci_quantiles)
      z_T_lower <- ci[1, ]
      z_T_upper <- ci[2, ]

      graphics::lines(
        T_vals_plot, z_T_upper,
        lty = 2, col = ggsci::pal_nejm()(1)[1]
      )
      graphics::lines(
        T_vals_plot, z_T_lower,
        lty = 2, col = ggsci::pal_nejm()(1)[1]
      )
    }
  }
  return(invisible(NULL))
}

#' @title CECL ggplot theme
#' @description Custom ggplot theme for CECL plots.
#' @param legend.position Position of legend in plot, default "bottom".
#' @param nejm_pal Logical indicating whether to use NEJM colour palette,
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
#' @param ... Additional arguments for main `ggplot2` plotting function.
#' @param environment Parent frame environment.
#' @inheritParams plot.cecl_marg
#' @return ggplot diagnostic plot for `cecl_marg` object.
#' @rdname ggplot.cecl_marg
#' @method ggplot cecl_marg
#' @export
#' @importFrom ggplot2 ggplot
#' @examples
#' library(ggplot2)
#'
#' # simulate some data
#' set.seed(123)
#' n_locs <- 5
#' n <- 1000
#' vars <- c("X1", "X2")
#' df <- do.call(
#'   rbind,
#'   lapply(1:n_locs, function(i) {
#'     data.frame(
#'       X1 = rgpd(n, u = 0, sigma = 1, xi = 0.2),
#'       X2 = rgpd(n, u = 0, sigma = 2, xi = -0.1),
#'       name = paste0("loc_", i)
#'     )
#'   })
#' )
#' # fit marginal models with quantile thresholding and ISMEV GPD fits
#' marg_fit <- cecl_marg(
#'   df,
#'   thresh_method = "quantile",
#'   thresh_args = 0.9,
#'   marg_method = "ismev",
#'   ncores = 1
#' )
#' # plot for a specific location and variable
#' ggplot(marg_fit, which = "qq", loc = "loc_1", var = "X1")
ggplot.cecl_marg <- \(
  data = NULL,
  mapping = ggplot2::aes(),
  which = c("qq", "pp", "hist", "return", "transformed"),
  loc,
  var,
  cond_var = NULL,
  mult_col = "name",
  plot_dens = TRUE,
  return_periods = c(1.5, 2.5, 5, 10, 20, 50, 100, 200),
  nboot = 200,
  refit = FALSE,
  ci_quantiles = c(0.025, 0.975),
  log_scale = TRUE,
  ...,
  environment = parent.frame()
) {
  stopifnot(inherits(data, "cecl_marg"))

  if (inherits(data, "cecl_marg_ecdf") && which != "transformed") {
    stop(paste(
      "No residuals to plot for 'ecdf' `marg_method`,",
      "only 'transformed' plot available."
    ))
  }

  which <- match.arg(which, several.ok = TRUE)

  if (missing(loc) || missing(var)) {
    stop("Please specify both 'loc' and 'var' to plot.")
  }

  if (!loc %in% names(data$original)) {
    stop(paste("Location", loc, "not found in the cecl_marg object."))
  }
  if (!var %in% data$vars) {
    stop(paste("Variable", var, "not found in the cecl_marg object."))
  }

  quantile <- lower <- upper <- density <- x <- theor <- model <- NULL

  # Extract original and thresholded data for specified location and variable
  orig_data <- data$original[[loc]] |>
    dplyr::select(dplyr::all_of(var))
  thresh_data <- data$data_thresh[[var]][[loc]] |>
    dplyr::select(dplyr::all_of(var))

  # Calculate residuals based on marginal method
  resid_fun <- \(q, gpd_params) do.call(
    pgpd,
    c(
      list(q = q),
      c(stats::setNames(
        gpd_params[c("sigma", "xi")],
        c("sigma", "xi")
      ), "u" = 0)
    )
  )

  if (any(which != "transformed")) {
    gpd_params <- data$marginal[[loc]][[var]]
    if (inherits(data, "cecl_marg_ismev")) {
      exceedances <- thresh_data[[var]] - gpd_params$thresh
      residuals <- resid_fun(exceedances, gpd_params)
    } else if (inherits(data, "cecl_marg_evgam")) {
      evgam_fit <- data$evgam_fit[[which(data$vars == var)]]
      pred_row <- evgam_fit$predictions |>
        dplyr::filter(.data[[mult_col]] == loc)
      sigma <- pred_row$scale
      exceedances <- thresh_data[[var]] - gpd_params$thresh
      residuals <- resid_fun(
        exceedances, list(sigma = sigma, xi = pred_row$shape)
      )
    } else {
      stop("ggplot method not implemented for this marg_method")
    }
    # check residuals
    stopifnot(
      "NA residuals - check pgpd arguments/parameterization" =
        !any(is.na(residuals))
    )
    stopifnot(
      "Some residuals are outside [0,1] - check parameterization" =
        all(residuals >= 0 & residuals <= 1)
    )
    res_df <- data.frame(residuals = residuals)

    if (!is.null(cond_var)) {
      message("Ignoring `cond_var` for which != 'transformed'")
    }
  }

  if (any(which == "transformed")) {
    stopifnot("must specify `cond_var`." = !is.null(cond_var))
    stopifnot("`cond_var` must differ from `var`." = cond_var != var)
    stopifnot(
      "`cond_var` not found in cecl_marg object." = cond_var %in% x$vars
    )
    transformed <- data$transformed[[loc]][, c(var, cond_var), drop = FALSE]
  }

  if (any(which %in% c("pp", "qq"))) {
    env <- bootstrap_pp_qq(
      exceedances = exceedances,
      gpd_params = gpd_params,
      type = intersect(which, c("pp", "qq")),
      nboot = nboot,
      refit = refit,
      ci_quantiles = ci_quantiles
    )
  }

  ret <- vector(mode = "list", length = length(which))

  # Generate specified ggplot
  # for (w in which) {
  for (i in seq_along(which)) {
    w <- which[[i]]
    if (w == "transformed") {
      p <- ggplot2::ggplot(
        as.data.frame(transformed),
        ggplot2::aes_string(x = cond_var, y = var)
      ) +
        ggplot2::geom_point(...) +
        ggplot2::labs(
          x = cond_var,
          y = var,
        ) +
        cecl_theme()
    } else if (w == "qq") {
      qq_df <- env$qq
      x_min <- min(qq_df$theor, na.rm = TRUE)
      x_max <- max(qq_df$theor, na.rm = TRUE)

      p <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theor, y = sample))
      if (nboot > 0) {
        p <- p +
          ggplot2::geom_ribbon(
            ggplot2::aes(ymin = lower, ymax = upper),
            fill = "grey80", alpha = 0.5
          )
      }
      p <- p +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red") +
        ggplot2::geom_point() +
        ggplot2::labs(x = "Theoretical", y = "Sample") +
        cecl_theme() +
        ggplot2::scale_x_continuous(
          limits = c(
            min(qq_df$theor, na.rm = TRUE),
            max(qq_df$theor, na.rm = TRUE)
          ),
          expand = c(0.01, 0.01)
        )
    } else if (w == "pp") {
      pp_df <- env$pp

      p <- ggplot2::ggplot(pp_df, ggplot2::aes(x = p, y = model))
      if (nboot > 0) {
        p <- p +
          ggplot2::geom_ribbon(
            ggplot2::aes(ymin = lower, ymax = upper),
            fill = "grey80", alpha = 0.5
          )
      }
      p <- p +
        ggplot2::geom_abline(intercept = 0, slope = 1, colour = "red") +
        ggplot2::geom_point() +
        ggplot2::labs(x = "Theoretical", y = "Model") +
        cecl_theme()
    } else if (w == "hist") {
      p <- ggplot2::ggplot(res_df, ggplot2::aes(x = residuals)) +
        ggplot2::labs(x = "Residuals") +
        cecl_theme()

      if (plot_dens) {
        p <- p +
          ggplot2::geom_histogram(
            ggplot2::aes(y = ggplot2::after_stat(density)),
            fill = "grey",
            colour = "black",
            ...
          ) +
          ggplot2::geom_density(
            colour = "red",
            linewidth = 1
          )
      } else {
        p <- p +
          ggplot2::geom_histogram(
            fill = "grey",
            colour = "black",
            ...
          )
      }
    } else if (w == "return") {
      # TODO Add predict method for marginal fits! would make this easier
      # TODO Add checking for return level arguments if using

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
      T_vals <- return_periods

      return_level <- \(T, u, sigma, xi, tol = 1e-8) {
        TT <- T * lambda_u
        if (abs(xi) > tol) {
          u + (sigma / xi) * (TT^xi - 1)
        } else {
          u + sigma * log(TT)
        }
      }

      # Compute nominal (fitted) return levels
      z_T <- vapply(
        T_vals, return_level, numeric(1),
        u = u, sigma = sigma, xi = xi
      )

      df <- data.frame(T_vals = T_vals, z_T = z_T)

      # bootstrap CI
      zT_list <- vector(mode = "list", length = nboot)
      b_ok <- 0 # count of successful bootstraps
      for (b in seq_len(nboot)) {
        sim_exc <- try(
          rgpd(n = n_exc, u = 0, sigma = sigma, xi = xi),
          silent = TRUE
        )
        if (inherits(sim_exc, "try-error")) {
          # fallback to inverse cdf
          sim_exc <- qgpd(stats::runif(n_exc), u = 0, sigma = sigma, xi = xi)
        }
        # fit GPD to simulated exceedances
        fit_b <- try(
          ismev::gpd.fit(sim_exc, threshold = 0, show = FALSE),
          silent = TRUE
        )
        if (inherits(fit_b, "try-error") || any(is.na(fit_b$mle))) next
        if (any(is.na(fit_b$mle))) next

        # extract bootstrap MLEs
        sigma_b <- fit_b$mle[1]
        xi_b <- fit_b$mle[2]

        # compute return levels using robust formula (handle xi_b ~ 0)
        zTb <- vapply(T_vals, function(T) {
          return_level(T, u = u, sigma = sigma_b, xi = xi_b)
        }, numeric(1))

        b_ok <- b_ok + 1
        zT_list[[b_ok]] <- zTb
      }

      # remove unused trailing NULLs if any
      zT_list <- zT_list[seq_len(b_ok)]

      if (b_ok > 0) {
        zT_boot <- do.call(rbind, zT_list)
        df$lower <- apply(zT_boot, 2, stats::quantile, probs = ci_quantiles[1])
        df$upper <- apply(zT_boot, 2, stats::quantile, probs = ci_quantiles[2])
      } else {
        df$lower <- NA_real_
        df$upper <- NA_real_
        message("No successful bootstrap fits for CI estimation.")
      }

      # ggplot
      p <- ggplot2::ggplot(df, ggplot2::aes(x = T_vals, y = z_T)) +
        ggplot2::geom_line(...) +
        ggplot2::geom_point() +
        ggplot2::labs(x = "Return Period", y = "Return Level") +
        cecl_theme()

      # Add CI ribbon if available
      if (!all(is.na(df$lower)) && !all(is.na(df$upper))) {
        p <- p + ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower, ymax = upper),
          alpha = 0.2,
          fill = ggsci::pal_nejm()(1)
        )
      }

      # convert to natural log scale
      if (log_scale == TRUE) {
        p <- p +
          ggplot2::scale_x_continuous(breaks = T_vals, transform = "log") +
          ggplot2::annotation_logticks(sides = "b") +
          ggplot2::labs(x = "Return Period (log)", y = "Return Level")
      }
    }
    ret[[i]] <- p
  }
  names(ret) <- which
  if (length(ret) == 1) {
    return(ret[[1]])
  } else {
    return(ret)
  }
}

#' @title `as_cecl_marg` method
#' @description Convert an object to class `cecl_marg`. The user must ensure
#' that the data provided has already been transformed to Laplace margins,
#' using whatever method is appropriate.
#' @param x Object to be transformed. Can either be a list of matrices (where
#' each matrix contains the transformed data for a group/location),
#'
#' by `cecl_marg`) of groups/locations, each containing a matrix of transformed
#' data.
#' @param ... Additional arguments passed to methods.
#' @return Object of class `cecl_marg`.
#' @rdname as_cecl_marg
#' @export
as_cecl_marg <- \(x, ...) {
  UseMethod("as_cecl_marg")
}

#' @title `as_cecl_marg` method for data.frames and tibbles
#' @description Convert a data.frame or tibble to class `cecl_marg`.
#' @param x Data.frame to be transformed.
#' @param name_col Name of column representing groups/locations, default
#' "name".
#' @param ... Additional arguments (not used).
#' @return Object of class `cecl_marg` and `cecl_marg_misc`.
#' @rdname as_cecl_marg
#' @method as_cecl_marg data.frame
#' @export
as_cecl_marg.data.frame <- \(x, name_col = "name", ...) {
  stopifnot(inherits(x, "data.frame"))
  stopifnot(
    "name_col must be in column names of x" = name_col %in% colnames(x)
  )

  x_fact <- x |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(name_col), \(y) factor(y, levels = unique(x[[name_col]]))
    ))

  ret <- x_fact |>
    dplyr::group_split(.data[[name_col]], .keep = FALSE) |>
    lapply(as.matrix)

  names(ret) <- levels(x_fact[[name_col]])
  ret <- list(
    "transformed" = ret,
    "vars"        = names(ret)[names(ret) != name_col]
  )
  class(ret) <- c("cecl_marg", "cecl_marg_user")
  ret
}

# TODO Could this just call data.frame method?
#' @title `as_cecl_marg` method for lists
#' @description Convert a list of matrices or dataframes/tibbles to class
#' `cecl_marg`.
#' @param x List of matrices or dataframes/tibbles to be transformed.
#' @param ... Additional arguments (not used).
#' @return Object of class `cecl_marg` and `cecl_marg_misc`.
#' @rdname as_cecl_marg
#' @method as_cecl_marg list
#' @export
as_cecl_marg.list <- \(x, ...) {
  stopifnot(inherits(x, "list"))

  # for a list of matrices
  if (all(vapply(x, is.matrix, logical(1)))) {
    cecl_marg_obj <- list(
      "transformed" = x,
      "vars"        = colnames(x[[1]]) # NOTE: Assumes no name_col
    )
    # for dataframes or tibbles, convert to matrices
  } else if (all(vapply(x, \(y) {
    inherits(y, c("data.frame", "tbl_df"))
  }, logical(1)))) {
    # check that all columns for each are numeric
    stopifnot("All columns must be numeric" = all(vapply(x, \(y) {
      all(vapply(y, is.numeric, logical(1)))
    }, logical(1))))

    cecl_marg_obj <- list(
      "transformed" = lapply(x, as.matrix), # TODO: Add optional var names
      "vars"        = names(x)
    )
  } else {
    stop("Input list must contain only matrices or dataframes/tibbles.")
  }

  class(cecl_marg_obj) <- c("cecl_marg", "cecl_marg_user")
  cecl_marg_obj
}
