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
#' @param marg_prob Can be a list of arguments to `evgam_ald_thresh` if a
#' variable quantile is desired, or a numeric value (scalar for a shared value
#' or a vector for each `vars`) for simple quantile thresholding.
#' @param marg_val Explicit value for marginal thresholds for each variable,
#' only specified if marg_prob is NULL.
#' @param f Formula for `evgam` model.
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
  thresh_args, # TODO Expand description
  thresh_only = FALSE,
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
    stopifnot(!names(thresh_args) %in% c("f", "qu", "jitter"))
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
  data_thresh <- lapply(data_thresh, bind_rows)

  # fit marginal models
  marginal <- fit_marg(
    data_df      = data_df,
    data_thresh  = data_thresh,
    vars         = vars,
    f            = f,
    loop_fun     = loop_fun
  )

  # return
  ret <- list(
    "marginal"    = marginal,
    "data_thresh" = data_thresh,
    "original"    = orig_dat,
    "vars"        = vars
  )
  # add evgam fit object if fitted
  # TODO Add to output of fit_marg
  if (exists("evgam_fit", envir = environment())) {
    names(evgam_fit) <- vars
    ret$evgam_fit <- evgam_fit
  }

  # make ret object of class `evc_marg`
  class(ret) <- c(
    "cecl_marg",
    paste0("cecl_marg_", thresh_method),
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
        group_split(name, .keep = TRUE)
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


# put below into function
fit_marg <- \(
  data_df,
  data_thresh,
  vars,
  f,
  loop_fun
) {
  # If f NULL, fit ordinary marginal models with `ismev::gpd.fit` for each loc
  if (is.null(f)) {
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
        # TODO Could replace with mapply!
        gpd_fits <- lapply(seq_along(vars), \(i) {
          fit <- ismev::gpd.fit(
            x[[vars[i]]],
            threshold = mth[i],
            # threshold = quantile(x[[vars[i]]], 0.9),
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
    # TODO Unsure if this works, test later
    names(marginal) <- purrr::map_chr(marginal, ~ as.character(.x[[1]]$name))

    # fit evgam model for each marginal
  } else {
    evgam_fit <- loop_fun(data_thresh, \(x) {
      fit_evgam(
        data      = x,
        pred_data = data_df,
        f         = f
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
