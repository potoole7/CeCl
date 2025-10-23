#' @title ???
#' @param data List of dataframes at each location containing data for each
#' variable. Optionally, can be an object of class `evc_marg` from a previous
#' call to `fit_ce`, with elements `marginal`, `data_thresh` and `original`.
#' @param vars Variable names for each location.
#' @param marg_prob Can be a list of arguments to `evgam_ald_thresh` if a
#' variable quantile is desired, or a numeric value (scalar for a shared value
#' or a vector for each `vars`) for simple quantile thresholding.
#' @param thresh_fun Function to threshold data, default is `evgam_ald_thresh`,
#' alternative is `qgam_thresh`.
#' @param marg_val Explicit value for marginal thresholds for each variable,
#' only specified if marg_prob is NULL.
#' @param f Formula for `evgam` model.
#' @param ncores Number of cores to use for parallel computation, Default: 1.
#' @return Object of type `cecl_marg` for each location.
#' @rdname cecl_marg
#' @importFrom rlang .data :=
#' @export

cecl_marg <- \(
  data,
  vars = NULL,
  marg_prob = list(
    f          = list("response ~ name", "~ name"), # must be as character
    tau        = .95,
    jitter     = TRUE
  ),
  thresh_fun = qgam_thresh,
  marg_val = NULL,
  f = list(excess ~ name, ~1), # keep shape constant for now
  ncores = 1
) {
  # number of variables
  nvars <- length(vars)

  # Parallel setup
  apply_fun <- ifelse(ncores == 1, lapply, parallel::mclapply)
  ext_args <- NULL
  if (ncores > 1) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }

  # if marg_prob used as args to thresh_fun, check args correct
  if (is.list(marg_prob)) {
    stopifnot(all(names(marg_prob) %in% names(formals(thresh_fun))))
    stopifnot(is.list(marg_prob$f))
    if (!all(vapply(marg_prob$f, is.character, logical(1)))) {
      stop(paste(
        "f should be a list of characters where 'response' is replaced by",
        "each specified 'vars'"
      ))
    }
  }

  # check marginal thresholds specified correctly
  if (is.null(marg_val) && is.null(marg_prob)) { # must provide one
    stop("you must provide one of marg_val or marg_prob")
  }
  if (!is.null(marg_val) && !is.null(marg_prob)) { # must provide only one
    stop("you must provide precisely one of marg_val or marg_prob")
  }
  # must have marginal values for each variable
  # TODO Expand so marg_val can be a list of locations like `start`
  # TODO Expand so marg_val can be a list of dataframes with varying thresh
  if (!is.null(marg_val) && is.numeric(marg_val)) {
    stopifnot(
      length(marg_val) == length(vars) && all(names(marg_val) == vars)
    )
  }

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
    # convert matrix names if required
    names(data_df)[names(data_df) %in% paste0("V", seq_len(nvars))] <- vars
  }

  # must have names column, and all of vars must be columns in data_df
  stopifnot("Must have a `name` column" = "name" %in% names(data_df))
  stopifnot("All of `vars` must be in data" = all(vars %in% names(data_df)))

  # locations (before thresholding)
  locs <- unique(data_df$name)
  data_df$name <- forcats::fct_inorder(data_df$name)

  ## Threshold Selection ##

  # If marg_val not specified, calculate thresh as quantile across all locs
  # TODO Allow marg_val to change by location
  if (is.null(marg_prob) || (!is.list(marg_prob) && is.numeric(marg_prob))) {
    # TODO Could also calculate for each location/name??
    if (is.null(marg_val)) {
      # marg_val <- apply(
      #   data_df[, c(vars)], 2, stats::quantile, marg_prob,
      #   na.rm = TRUE
      # )
      marg_val <- mapply(
        quantile,
        data_df[, vars],
        marg_prob, # can be scalar or vector
        MoreArgs = list(na.rm = TRUE)
      )
      names(marg_val) <- vars
    }
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
    # If thresh is a list, assume it is arguments to thresh_fun
    # TODO: May be easier to just copy each vars column as response in data_df
    # Would allow for simpler formula specification
  } else if (is.list(marg_prob)) {
    data_thresh <- lapply(vars, \(x) {
      print(paste0("thresholding ", x))

      # Change formula to include response in question
      spec_params <- marg_prob
      spec_params$f <- lapply(marg_prob$f, \(f_spec) {
        stats::formula(stringr::str_replace_all(f_spec, "response", x))
      })
      # allow different thresholds for each variable
      if (length(spec_params$tau) > 1) {
        spec_params$tau <- spec_params$tau[vars == x]
      }

      # Run thresholding function for each response with specified args
      # TODO Also return evgam fits to ald
      quantile_fits <- do.call(
        thresh_fun,
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
  } else {
    stop("marg_prob must be numeric or arguments to specified marg_fun")
  }
  # Only keep locs with exceedances for all vars, otherwise can't do CE!
  # TODO evc only works for locations, add error if data not in this form!
  locs_keep <- Reduce(intersect, lapply(data_thresh, \(x) unique(x$name)))
  data_df <- dplyr::filter(data_df, name %in% locs_keep)
  # remove duplicate threshold rows kept through floating point errors
  data_thresh <- lapply(data_thresh, \(x) {
    dplyr::filter(x, name %in% locs_keep) |>
      dplyr::group_by(
        name,
        dplyr::across(dplyr::any_of(c("date", !!vars)))
      ) |>
      dplyr::slice(1) |>
      dplyr::ungroup()
  })
  # TODO Add option to only return thresholded data!
  if (thresh_only) {
    return(list("data_thresh" = data_thresh, "original" = data_df))
  }

  ## Marginal Model ##

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

  # original data
  orig_dat <- data_df |>
    dplyr::group_by(name) |>
    dplyr::group_split(.keep = TRUE)

  names(orig_dat) <- purrr::map_chr(orig_dat, ~ as.character(.x$name[1]))

  # names(marginal) <- locs_keep # TODO Check that this is correct!!!
  names(data_thresh) <- vars
  ret <- list(
    "marginal"    = marginal,
    "data_thresh" = data_thresh,
    "original"    = orig_dat,
    "vars"        = vars
  )
  # add evgam fit object if fitted
  if (exists("evgam_fit", envir = environment())) {
    names(evgam_fit) <- vars
    ret$evgam_fit <- evgam_fit
  }
  # make ret object of class `evc_marg`
  class(ret) <- c("cecl_marg", class(ret))
  return(ret)
}
