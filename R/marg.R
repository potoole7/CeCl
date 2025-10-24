#' @title Fit marginal models for multivariate conditional extremes
#' @description Fit marginal models for multivariate conditional extremes.
#' @param data List of dataframes at each location containing data for each
#' variable. Optionally, can be an object of class `evc_marg` from a previous
#' call to `fit_ce`, with elements `marginal`, `data_thresh` and `original`.
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
  vars = NULL,
  thresh_method = c("value", "quantile", "regression", "none"),
  thresh_args = NULL, # TODO Expand description
  thresh_only = FALSE,
  marg_method = c("ecdf", "ismev", "evgam"),
  marg_args = NULL,
  # marg_prob = list(
  #   f          = list("response ~ name", "~ name"), # must be as character
  #   tau        = .95,
  #   jitter     = TRUE
  # ),
  # marg_val = NULL,
  # f = list(excess ~ name, ~1), # keep shape constant for now
  ncores = 1
) {
  ## Initial ##

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

  # if (is.list(marg_prob)) {
  #   stopifnot(is.list(marg_prob$f))
  #   if (!all(vapply(marg_prob$f, is.character, logical(1)))) {
  #     stop(paste(
  #       "f should be a list of characters where 'response' is replaced by",
  #       "each specified 'vars'"
  #     ))
  #   }
  # }

  # # check marginal thresholds specified correctly
  # if (is.null(marg_val) && is.null(marg_prob)) { # must provide one
  #   stop("you must provide one of marg_val or marg_prob")
  # }
  # if (!is.null(marg_val) && !is.null(marg_prob)) { # must provide only one
  #   stop("you must provide precisely one of marg_val or marg_prob")
  # }

  # # must have marginal values for each variable
  # # TODO Expand so marg_val can be a list of locations like `start`
  # # TODO Expand so marg_val can be a list of dataframes with varying thresh
  # if (!is.null(marg_val) && is.numeric(marg_val)) {
  #   stopifnot(
  #     length(marg_val) == length(vars) && all(names(marg_val) == vars)
  #   )
  # }

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
  prep_out <- prep_df(data, vars)
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
    vars          = vars,
    thresh_method = thresh_method,
    thresh_args   = thresh_args
  )
  data_thresh <- thresh_out[[1]]
  locs_keep <- thresh_out[[2]] # locations with exceedances for all variables
  names(data_thresh) <- vars

  # Only keep locs with exceedances for all vars, otherwise can't do CE!
  data_df <- dplyr::filter(data_df, name %in% locs_keep)

  # original data split by location (for returning later)
  orig_dat <- data_df |>
    dplyr::group_by(name) |>
    dplyr::group_split(.keep = TRUE)
  # add names to list elements
  names(orig_dat) <- purrr::map_chr(orig_dat, ~ as.character(.x$name[1]))

  # return thresholded data only if specified
  # TODO Make into S3 class??
  if (thresh_only) {
    ret <- list(
      "data_thresh" = data_thresh,
      "original"    = orig_dat,
      "vars"        = vars
    )
    class(ret) <- c(
      "cecl_thresh",
      paste0("cecl_thresh_", thresh_method),
      class(ret)
    )

    return(ret)
  }

  # reverse splitting of data_thresh into locations
  data_thresh <- lapply(data_thresh, \(x) bind_rows(unname(x)))

  # fit marginal models
  marginal <- fit_marg(
    data_df      = data_df, # TODO data_df will be weird if thresh = none
    data_thresh  = data_thresh,
    vars         = vars,
    marg_method  = marg_method,
    marg_args    = marg_args,
    # f            = f,
    loop_fun     = loop_fun
  )
  names(marginal) <- locs_keep

  # Transform data to Laplace margins
  marginal_trans <- trans_marg(marginal, data_df, vars)
  names(marginal_trans) <- locs_keep

  # return
  ret <- list(
    "marginal"    = marginal,
    "data_thresh" = data_thresh,
    "original"    = orig_dat,
    "transformed" = marginal_trans,
    "vars"        = vars
  )
  # add evgam fit object if fitted
  # TODO Add to output of fit_marg
  # TODO Fix
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

## Data Prep ##

#' @title Prepare data frame for cecl_marg
#' @description Prepare data frame for cecl_marg: Convert list of data into
#' dataframe with `name` columne specifying multivariate structure.
#' @inheritParams cecl_marg
#' @return Data frame ready for `cecl_marg`.
#' @rdname prep_df
#' @keywords internal
prep_df <- function(data, vars) {
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
      if (!"name" %in% names(ret)) {
        ret <- ret |>
          dplyr::mutate(name = paste0("location_", i))
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
      vars <- names(data_df)[!names(data_df) %in% "name"]
    }
    # convert matrix names if required
    names(data_df)[names(data_df) %in% paste0("V", seq_len(nvars))] <- vars
  }

  # must have names column, and all of vars must be columns in data_df
  stopifnot("Must have a `name` column" = "name" %in% names(data_df))
  stopifnot("All of `vars` must be in data" = all(vars %in% names(data_df)))

  # make name a factor in order of appearance
  data_df$name <- forcats::fct_inorder(data_df$name)

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
  vars,
  thresh_method,
  thresh_args
) {
  # locations (before thresholding)
  locs <- unique(data_df$name)

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
        group_split(name, .keep = TRUE)
      names(ret) <- purrr::map_chr(ret, ~ as.character(.x$name[1]))
      return(ret)
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
        # group_split(name, .keep = TRUE) |>
        identity()
    })
    # If thresh is a list, assume it is arguments to thresh_fun (now qgam_thresh)
    # TODO: May be easier to just copy each vars column as response in data_df
    # Would allow for simpler formula specification
  }

  if (thresh_method == "regression") {
    data_thresh <- lapply(vars, \(x) {
      print(paste0("thresholding ", x))

      # Change formula to include response in question
      # spec_params <- marg_prob
      spec_params <- thresh_args
      # spec_params$f <- lapply(marg_prob$f, \(f_spec) {
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
      thresh_locs <- unique(ret$name)
      if (length(thresh_locs) < length(locs)) {
        loc_missing <- setdiff(locs, thresh_locs)
        message(paste0(
          "No exceedances for variable ", x, " at: ",
          paste(loc_missing, collapse = ", "),
          ", removing for all variables"
        ))
      }
      return(ret)
    })
  }

  # locations with exceedances for all variables
  locs_keep <- Reduce(intersect, lapply(data_thresh, \(x) unique(x$name)))

  # remove duplicate threshold rows kept through floating point errors
  data_thresh <- lapply(data_thresh, \(x) {
    ret <- dplyr::filter(x, name %in% locs_keep) |>
      dplyr::group_by(
        name,
        dplyr::across(dplyr::any_of(c("date", !!vars)))
      ) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      # also split by location
      group_split(name, .keep = TRUE)

    names(ret) <- purrr::map_chr(ret, ~ as.character(.x$name[1]))
    return(ret)
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
    dplyr::distinct(across(c(all_of(predictors), thresh)))

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
  vars,
  marg_method,
  marg_args,
  loop_fun
) {
  if (marg_method == "ecdf") {
    locs <- unique(data_df$name)
    # dummy GPD fits
    marginal <- lapply(locs, \(loc) {
      gpd_y <- lapply(vars, \(var) {
        list(
          "sigma" = NA_real_,
          "xi" = NA_real_,
          # take as threshold the maximum value (to transform using only ECDF)
          "thresh" = data_thresh[[var]] |>
            dplyr::filter(name == loc) |>
            arrange(!!rlang::sym(var)) |>
            dplyr::slice(n()) |>
            dplyr::pull(!!rlang::sym(var))
        )
      })
      names(gpd_y) <- vars
      return(gpd_y)
    })
    names(marginal) <- locs
  }

  if (marg_method == "ismev") {
    # calculate for all locations
    marginal <- data_df |>
      dplyr::group_split(name, .keep = TRUE) |>
      loop_fun(\(x) {
        # pull marginal thresholds
        mth <- vapply(data_thresh, \(y) {
          y |>
            # need thresh for correct loc
            dplyr::filter(name == x$name[[1]]) |>
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
          return(list(
            "sigma"     = fit$mle[1],
            "xi"        = fit$mle[2],
            "thresh"    = fit$threshold[[1]],
            # "name"      = fit$name[1]
            "name"      = x$name[1]
          ))
        })
        names(gpd_fits) <- vars
        return(gpd_fits)
      })

    # add names (correctly!)
    names(marginal) <- purrr::map_chr(marginal, ~ as.character(.x[[1]]$name))
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
        name,
        sigma = scale,
        xi = shape
      ) |>
        dplyr::distinct(signif(sigma, 6), signif(xi, 6), .keep_all = TRUE) |>
        # also pull in thresholds
        dplyr::left_join(
          dplyr::select(data_thresh[[i]], name, thresh) |>
            dplyr::distinct(name, signif(thresh, 6), .keep_all = TRUE)
        ) |>
        dplyr::select(-dplyr::matches("signif"))

      # split into list by name as dependence pars will also be this way
      # loc_names_spec <- unique(params_df$name)
      out <- params_df |>
        # dplyr::group_split(dplyr::row_number(), .keep = FALSE) |>
        # dplyr::group_split(name, .keep = FALSE) |>
        dplyr::relocate(name, .after = dplyr::last_col()) |>
        dplyr::group_split(name, .keep = TRUE) |>
        # setNames(loc_names_spec) |> # wrong, causes bug!!
        lapply(as.vector, mode = "list")
      # pull names correctly from list objects
      names(out) <- purrr::map_chr(out, ~ as.character(.x$name))
      return(out)
    })
    names(marginal) <- vars

    # transpose list from variables -> locations to locations -> variables
    marginal <- purrr::transpose(marginal)
  }

  return(marginal)
}

#' @title Transform data to Laplace margins for cecl_marg
#' @description Transform data to Laplace margins for cecl_marg.
#' @param marginal List of fitted marginal models for each location.
#' @param data_df Data frame of original data.
#' @param vars Names of variable columns.
#' @return List of data matrices transformed to Laplace margins for each location.
trans_marg <- \(
  marginal,
  data_df,
  vars
) {
  # Calculate dependence from marginals (default output object)
  # first, transform margins to Laplace
  lapply(seq_along(marginal), \(i) {
    # semi-parametric CDF
    # TODO: More efficient to also split data_df by name and subset with i
    F_hat <- data_df |>
      # dplyr::filter(name == locs_keep[i]) |>
      dplyr::filter(name == names(marginal)[i]) |>
      dplyr::select(dplyr::all_of(vars)) |>
      p_gpd_ecdf(marginal[[i]])
    # Laplace transform
    Y <- dlaplace(F_hat)
    colnames(Y) <- vars
    return(Y)
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
    return(cdf)
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
d_gpd_ecdf <- function(F_hat, dat, gpd) {
  return(vapply(seq_along(gpd), function(i) {
    dat_spec <- dat[, i, drop = TRUE]
    stopifnot(names(gpd[[i]]) == c("sigma", "xi", "thresh"))

    spec_sigma <- gpd[[i]][[1]]
    spec_xi <- gpd[[i]][[2]]
    spec_loc <- gpd[[i]][[3]]

    n <- length(dat_spec)
    probs <- (1:n) / (n + 1) # Empirical CDF probabilities

    # Find closest probability match for each F_hat value
    px <- vapply(F_hat[, i], function(x, p) {
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

    return(res)
  }, FUN.VALUE = numeric(nrow(F_hat))))
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
  return(ret)
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
    return(ifelse(y < 0.5, log(2 * y), -log(2 * (1 - y))))
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
#' @description Fit and generate predictions from `evgam` model
#' @param data Dataframe for one location.
#' @param pred_data Dataframe for one location to predict on.
#' @param f Formula for `evgam` model.
#' @return List with model `m` and predictions `predictions`.
#' @rdname fit_evgam
#' @keywords internal
fit_evgam <- \(
  data,
  pred_data,
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
    dplyr::distinct(name, dplyr::across(dplyr::all_of(predictors)))
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
