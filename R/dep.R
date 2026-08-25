#' @title Fit conditional extremes dependence model
#' @description Fit the conditional extremes dependence model of Heffernan
#' and Tawn (2004) to data transformed to Laplace scale.
#' @param obj Object of class `cecl_marg` containing transformed data.
#' @param cond_prob Numeric vector of quantiles to use for conditioning
#' variables.
#' @param cond_val Numeric vector of thresholds on the Laplace scale to use
#' for conditioning variables. Exactly one of `cond_prob` and `cond_val`
#' must be supplied.
#' @param vars Character vector of variable names to fit dependence model for.
#' If `NULL`, all variables are used.
#' @param cond_vars Character vector of variable names to condition on.
#' If `NULL`, all variables are used.
#' @param start Named numeric vector or list of named numeric matrices
#' with starting values for dependence parameters `a` and `b`.
#' If a list, must have one entry per location.
#' @param ncores Number of cores to use for parallel processing.
#' @param nruns Number of optimization runs with different starting values.
#' @param aLow Lower bound for dependence parameter `a`.
#' @param fixed_b Logical indicating whether to fix dependence parameter `b`
#' at start value.
#' @param fit_no_keef Logical indicating whether to fit the model without
#' Keef et al. (2012) constraints.
#' @return Object of class `cecl_dep` containing fitted dependence parameters
#' and residuals.
#' @export
cecl_dep <- \(
  obj,
  cond_prob = NULL,
  cond_val = NULL,
  # TODO Need to add more checking, can add unwanted cols (same in marg)
  vars = NULL,
  cond_vars = NULL,
  # TODO Need to be able to have vector of length vars here though!
  start = c("a" = 0.01, "b" = 0.01),
  ncores = 1,
  nruns = 1,
  aLow = -1,
  fixed_b = FALSE,
  fit_no_keef = FALSE
) {
  stopifnot(inherits(obj, "cecl_marg"))

  cond_var <- var <- name <- a <- b <- NULL

  # must have one of cond_prob or cond_val
  if (is.null(cond_prob) && is.null(cond_val) ||
    sum(!is.null(cond_prob), !is.null(cond_val)) > 1
  ) {
    stop("Must specify one of cond_prob or cond_val")
  }

  # Parallel setup
  # TODO Functionalise, used in multiple places
  apply_fun <- ifelse(ncores == 1, lapply, parallel::mclapply)
  ext_args <- NULL
  if (ncores > 1) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }

  # extract relative information from obj
  marginal_trans <- obj$transformed
  locs_keep <- names(marginal_trans)
  if (is.null(vars)) {
    vars <- obj$vars
  }
  stopifnot("vars not in data, check again" = !is.null(vars))

  # conditioning variables default to all
  if (is.null(cond_vars)) {
    cond_vars <- vars
  }
  stopifnot(
    "cond_vars not in variables, check again" = all(cond_vars %in% vars)
  )

  # only keep variables of interest from marginal object
  marginal_trans <- lapply(marginal_trans, \(x) {
    x[, unique(c(vars, cond_vars))]
  })

  # start values must either be named vector or dataframe from `coef(dep)`
  # TODO Change argument documentation above
  # TODO Expand to also allow starting values for m nd s
  is_df_start <- FALSE
  # for vector start values, check that they are named and have correct names
  if (is.vector(start)) {
    cond <- is.numeric(start) && !is.null(names(start)) &&
      all(names(start) == c("a", "b"))
    stopifnot(
      "`start` vector must have names a and b" = cond
    )
    # check for data.frame/coef.cecl_dep
    # } else if (is.data.frame(start) && !"coef.cecl_dep" %in% class(start)) {
  } else if (is.data.frame(start)) { # will pass this if `coef.cecl_dep`
    is_df_start <- TRUE
    rq_cols <- c("name", "var", "cond_var", "a", "b")
    if (!all(rq_cols %in% colnames(start))) {
      msg <- paste0(
        "`start` data.frame must have columns: ",
        paste(rq_cols, collapse = ", "),
        ", consider using `coef.cecl_dep` to extract dependence parameters",
        " from fitted model"
      )
      stop(msg)
    }
    # filter to only relevant rows of start data, if cond_vars specified
    start <- start |>
      dplyr::filter(cond_var %in% cond_vars)

    if (nrow(start) == 0L) {
      stop(
        "No rows of `start` data.frame match specified `cond_vars`; ",
        "check that `start` contains matching rows or adjust `cond_vars`.",
        call. = FALSE
      )
    }
  } else {
    msg <- paste0(
      "`start` must be either a vector with names a and b, or a data.frame",
      " with columns name, var, cond_var, a and b (e.g. from `coef.cecl_dep`)"
    )
    stop(msg)
  }

  # fit dependence model to transformed data for each location
  dependence <- loop_fun(seq_along(marginal_trans), \(i) {
    # Quantiles/thresholds are indexed by conditioning variable inside
    # ce_optim(), and the same specification is used at every location.
    cond_prob_spec <- cond_prob
    cond_val_spec <- cond_val

    # pull start values
    start_spec <- start # if a vector, just use these
    if (is_df_start) {
      start_spec <- start |>
        dplyr::filter(
          name == locs_keep[i],
          var %in% vars,
          cond_var %in% cond_vars
        ) |>
        dplyr::select(a, b, var, cond_var)

      # validate that filtering produced at least one row of start values
      if (nrow(start_spec) == 0L) {
        stop(
          sprintf(
            paste0(
              "No starting values found in `start` for location '%s' ",
              "with vars '%s' and cond_vars '%s'. ",
              "Please check that `start` specifies rows matching ",
              "`name`, `var`, and `cond_var`."
            ),
            locs_keep[i],
            paste(vars, collapse = ", "),
            paste(cond_vars, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }

    # fit dependence model
    o <- ce_optim(
      Y         = marginal_trans[[i]],
      dqu       = cond_prob_spec,
      dth       = cond_val_spec,
      cond_vars = cond_vars,
      control   = list(maxit = 1e6),
      constrain = !fit_no_keef,
      aLow      = aLow,
      fixed_b   = fixed_b,
      start     = start_spec,
      nruns     = nruns
    )
    o
  })

  # recursively pull dependence parameters and residuals out separately
  pull_element <- \(x, element) {
    if (element %in% names(x)) {
      x[[element]]
    } else {
      lapply(x, pull_element, element)
    }
  }
  ret <- lapply(c("resid", "params"), \(x) {
    stats::setNames(pull_element(dependence, x), locs_keep)
  })
  names(ret) <- c("residual", "dependence")
  # add other information to return object
  ret$call <- match.call()
  ret$transformed <- obj$transformed
  ret$start <- start

  # check that all dependence models have run successfully, message if not
  locs_fail <- locs_keep[
    vapply(ret$dependence, \(x) any(is.na(unlist(x))), logical(1))
  ]
  if (length(locs_fail) > 0) {
    message(paste0(
      length(locs_fail),
      " locations failed to fit CE model for at least one variable: ",
      paste(locs_fail, collapse = ", ")
    ))
  }

  # set class and return
  class(ret) <- c(
    "cecl_dep",
    class(ret)
  )
  return(ret)
}

# check constraints on parameters under constrained Laplace estimation
constraints_satisfied <- \(a, b, z, zpos, zneg, v) {
  C1e <- a <= min(
    1, 1 - b * min(z) * v^(b - 1), 1 - v^(b - 1) * min(z) + min(zpos) / v
  ) &
    a <= min(
      1, 1 - b * max(z) * v^(b - 1), 1 - v^(b - 1) * max(z) + max(zpos) / v
    )

  C1o <- a <= 1 &
    a > 1 - b * min(z) * v^(b - 1) &
    a > 1 - b * max(z) * v^(b - 1) &
    (1 - 1 / b) * (b * min(z))^(1 / (1 - b)) *
      (1 - a)^(-b / (1 - b)) + min(zpos) > 0 &
    (1 - 1 / b) * (b * max(z))^(1 / (1 - b)) *
      (1 - a)^(-b / (1 - b)) + max(zpos) > 0

  C2e <- -a <= min(
    1, 1 + b * v^(b - 1) * min(z), 1 + v^(b - 1) * min(z) - min(zneg) / v
  ) &
    -a <= min(
      1, 1 + b * v^(b - 1) * max(z), 1 + v^(b - 1) * max(z) - max(zneg) / v
    )

  C2o <- -a <= 1 &
    -a > 1 + b * v^(b - 1) * min(z) &
    -a > 1 + b * v^(b - 1) * max(z) &
    (1 - 1 / b) * (-b * min(z))^(1 / (1 - b)) *
      (1 + a)^(-b / (1 - b)) - min(zneg) > 0 &
    (1 - 1 / b) * (-b * max(z))^(1 / (1 - b)) *
      (1 + a)^(-b / (1 - b)) - max(zneg) > 0

  if (any(is.na(c(C1e, C1o, C2e, C2o)))) {
    message("Strayed into impossible area of parameter space")
    C1e <- C1o <- C2e <- C2o <- FALSE
  }

  (C1e | C1o) && (C2e | C2o)
}

# function to evaluate (negative log) likelihood
laplace_nll <- \(yex, ydep, a, b, m, s, constrain, v, aLow) {
  BigNumber <- 10^40
  WeeNumber <- 10^(-10)

  # give large value if parameters are out of bounds
  if (a < aLow || s < WeeNumber || a > 1 - WeeNumber || b > 1 - WeeNumber) {
    return(BigNumber)
  } else {
    # Assuming normal distribution for excesses, calculate mean & sd
    mu <- a * yex + m * yex^b
    sig <- s * yex^b

    # calculate log likelihood
    res <- sum(0.5 * log(2 * pi) + log(sig) + 0.5 * ((ydep - mu) / sig)^2)

    if (is.infinite(res)) {
      if (res < 0) {
        return(-BigNumber)
      } else {
        return(BigNumber)
      }
      warning("Infinite value of Q in mexDependence")
      # Apply Keef, Papastathopoulos constraints if specified
    } else if (constrain) {
      zpos <- range(ydep - yex) # q0 and q1
      z <- range((ydep - yex * a) / (yex^b)) # q0 and q1
      zneg <- range(ydep + yex) # q0 and q1

      if (!constraints_satisfied(a, b, z, zpos, zneg, v)) {
        return(BigNumber)
      }
    }
  }
  res
}

# function to evaluate (negative) profile (log) likelihood and optimise over
laplace_npll <- \(yex, ydep, a, b, constrain, v, aLow) {
  # first, estimate Z by rearranging the conditional extremes equation
  Z <- (ydep - yex * a) / (yex^b)
  stopifnot(
    "NaNs in Z, conditional quantile may be negative" = all(!is.nan(Z))
  )
  # estimate nuisance parameters
  m <- mean(Z)
  s <- stats::sd(Z)

  # now estimate a and b
  res <- laplace_nll(
    yex, ydep, a, b,
    m = m, s = s, constrain, v, aLow = aLow
  )
  res <- list(profLik = res, m = m, s = s)
  res
}

# function to evaluate profile likelihood and optimise over
Qpos <- \(param, yex, ydep, constrain, v, aLow) {
  a <- param[1]
  b <- param[2]

  res <- laplace_npll(yex, ydep, a, b, constrain, v, aLow)
  res$profLik
}

Qpos_fixed_b <- \(param, yex, ydep, constrain, v, aLow, b) {
  a <- param[1]
  res <- laplace_npll(yex, ydep, a, b, constrain, v, aLow)
  res$profLik
}

# fit CE for all pairs of variables
ce_optim <- \(
  Y,
  dqu = NULL,
  dth = NULL,
  cond_vars = NULL,
  start = c("a" = 0.01, "b" = 0.01),
  control = list(maxit = 1e6),
  constrain = TRUE,
  v = 10,
  aLow = -1,
  fixed_b = FALSE,
  nruns = 2
) {

  # set objects to NULL to appease R CMD check
  var <- cond_var <- a <- b <- NULL


  # must specify either dqu or dth
  if (is.null(dqu) && is.null(dth) ||
    sum(!is.null(dqu), !is.null(dth)) > 1) {
    stop("Must specify either dependence quantile (dqu) or threshold (dth)")
  }

  if (!is.data.frame(start)) {
    stopifnot(
      "`start` must be a named numeric vector with elements 'a' and 'b'" =
        is.numeric(start) &&
          is.vector(start) &&
          identical(names(start), c("a", "b"))
    )
  }

  # check that Y has names; if not give dummy names
  names_y <- colnames(Y)
  ncol_y <- ncol(Y)
  if (is.null(names_y)) {
    names_y <- paste0("var_", seq_len(ncol_y))
    colnames(Y) <- names_y
  }

  # conditioning variables default to all variables
  if (is.null(cond_vars)) {
    cond_vars <- names_y
  }
  if (!all(cond_vars %in% names_y)) {
    stop("All `cond_vars` must be columns of `Y`.", call. = FALSE)
  }
  if (anyDuplicated(cond_vars)) {
    stop("`cond_vars` must not contain duplicates.", call. = FALSE)
  }
  if (!is.null(dqu) && !length(dqu) %in% c(1L, length(cond_vars))) {
    stop(
      "`dqu` must have length 1 or one value per conditioning variable.",
      call. = FALSE
    )
  }
  if (!is.null(dth) && !length(dth) %in% c(1L, length(cond_vars))) {
    stop(
      "`dth` must have length 1 or one value per conditioning variable.",
      call. = FALSE
    )
  }

  # optimise for a single variable vs another
  single_optim <- \(yex, ydep, start, dqu, dth) {
    # browser()
    # threshold data
    thresh <- dth # if threshold value provided, use that
    # otherwise, threshold at quantile
    if (is.null(dth)) {
      thresh <- stats::quantile(yex, dqu)
    }
    wch <- yex > thresh

    # object to return if an error is found
    err_obj <- list(
      "resid" = matrix(NA, nrow = max(sum(wch), 1)), # can't be 0
      "params" = c(
        "a" = NA, "b" = NA, "m" = NA, "s" = NA, "ll" = NA, "dth" = NA
      )
    )

    if (any(is.infinite(yex))) {
      message("Inf values in Laplace transformed data, optimisation failed")
      return(err_obj)
    }

    if (sum(wch) == 0) {
      message("No exceedances above threshold, optimisation failed")
      return(err_obj)
    }

    # if fixing b, do 1D optimisation on a only
    if (fixed_b == TRUE) {
      # TODO May have to change to work with list??
      b <- start[[2]]
      start <- start[1]
      o_single <- try(stats::optim(
        par       = start,
        fn        = Qpos_fixed_b,
        method    = "Brent", # 1D optimisation
        lower     = -1,
        upper     = 1,
        control   = control,
        yex       = yex[wch],
        ydep      = ydep[wch],
        constrain = constrain,
        v         = v,
        aLow      = aLow,
        b         = b
      ), silent = TRUE)
      if (!inherits(o_single, "try-error")) {
        o_single$par <- c(o_single$par, b)
      }
      # else optimise a and b together, as normal
    } else {
      o_single <- try(stats::optim(
        par       = start,
        fn        = Qpos,
        control   = control,
        yex       = yex[wch],
        ydep      = ydep[wch],
        constrain = constrain,
        v         = v,
        aLow      = aLow
      ), silent = TRUE)
    }

    if (inherits(o_single, "try-error")) {
      message(paste("optimisation failed for constrain =", constrain))
      return(err_obj)
    }
    # set to NA if no change from (default!) starting values
    if (all(o_single$par[1:2] == start) && all(start == 0.01)) {
      message("No change from starting values, optimisation failed")
      return(err_obj)
    }
    # back-calculate residuals and nuisance parameters from a and b estimates
    if (all(!is.na(o_single$par))) {
      Z <- (ydep[wch] - yex[wch] * o_single$par[1]) /
        (yex[wch]^o_single$par[2])
      o_single$par <- c(o_single$par[1:2], mean(Z), stats::sd(Z))
    } else {
      Z <- rep(NA_real_, sum(wch))
      o_single$par <- c(o_single$par, NA, NA)
    }
    # add LL and numerical (Laplace scale) threshold, label appropriately
    o_single$par <- c(o_single$par, o_single$value, thresh[[1]])
    names(o_single$par) <- names(err_obj$params)
    # TODO Add checks afterwards on o
    list("resid" = matrix(Z), "params" = o_single$par)
  }

  # loop through RHS conditioning variables
  ret <- lapply(seq_along(cond_vars), \(i) {
    cond_var_i <- cond_vars[[i]]
    which_cond <- match(cond_var_i, names_y)
    names_conditioned <- names_y[-which_cond]

    # can have different dependence quantiles/values for each variable
    dqu_spec <- dqu
    if (length(dqu) > 1L) {
      dqu_spec <- dqu[[i]]
    }
    dth_spec <- dth
    if (length(dth) > 1L) {
      dth_spec <- dth[[i]]
    }

    # RHS of the CE equation: conditioning variable
    yex <- Y[, which_cond, drop = TRUE]

    # loop through LHS conditioned variables
    o_yex <- lapply(seq_along(names_conditioned), \(j) {
      conditioned_var_j <- names_conditioned[[j]]
      which_conditioned <- match(conditioned_var_j, names_y)
      ydep <- Y[, which_conditioned, drop = TRUE]

      # extract specific start values (if data.frame)
      start_spec <- start
      if (is.data.frame(start_spec)) {
        start_spec <- start_spec |>
          dplyr::filter(
            var == conditioned_var_j,
            cond_var == cond_var_i
          ) |>
          dplyr::select(a, b) |>
          unlist(use.names = TRUE)
        if (!identical(names(start_spec), c("a", "b"))) {
          stop(
            sprintf(
              paste0(
                "No unique starting values found for conditioned variable ",
                "'%s' and conditioning variable '%s'."
              ),
              conditioned_var_j,
              cond_var_i
            ),
            call. = FALSE
          )
        }
      }
      o <- single_optim(
        yex = yex,
        ydep = ydep,
        start = start_spec,
        dqu = dqu_spec,
        dth = dth_spec
      )

      # check if optimisation failed
      if (all(is.na(o$params))) {
        return(o)
      }

      # perform multiple runs with prev estimates as start values, if desired
      if (nruns > 1L) {
        for (run in seq_len(nruns - 1L)) {
          start_spec <- o$params[c("a", "b")]
          o <- single_optim(
            yex = yex,
            ydep = ydep,
            start = start_spec,
            dqu = dqu_spec,
            dth = dth_spec
          )
          if (all(is.na(o$params))) {
            break
          }
        }
      }
      o
    })

    # join matrices; columns are LHS conditioned variables
    o_yex <- list(
      "resid"  = do.call(cbind, lapply(o_yex, `[[`, "resid")),
      "params" = do.call(cbind, lapply(o_yex, `[[`, "params"))
    )

    colnames(o_yex$resid) <- names_conditioned
    colnames(o_yex$params) <- names_conditioned
    o_yex
  })
  # outer list elements are RHS conditioning variables
  names(ret) <- cond_vars
  ret
}

#' @title Extract dependence parameters from `cecl_dep` object
#' @description Extract dependence parameters from a fitted
#' `cecl_dep` object.
#' @param object Object of class `cecl_dep`.
#' @param ... Additional arguments (not used).
#' @return Data frame of dependence parameters for each location and
#' conditioned variable.
#' @rdname coef.cecl_dep
#' @export
#' @method coef cecl_dep
coef.cecl_dep <- \(object, ...) {
  stopifnot(inherits(object, "cecl_dep"))

  # Storage invariant:
  # dependence[[location]][[conditioning variable]][, conditioned variable]
  ret <- do.call(rbind, lapply(names(object$dependence), \(loc) {
    dep_loc <- object$dependence[[loc]]
    do.call(rbind, lapply(names(dep_loc), \(cond_var_i) {
      params <- dep_loc[[cond_var_i]]
      do.call(rbind, lapply(colnames(params), \(var_i) {
        data.frame(
          name = loc,
          var = var_i,
          cond_var = cond_var_i,
          as.list(params[, var_i, drop = TRUE]),
          row.names = NULL,
          check.names = FALSE
        )
      }))
    }))
  }))
  rownames(ret) <- NULL
  class(ret) <- c("coef.cecl_dep", class(ret))
  return(ret)
}

#' @title Extract residuals from `cecl_dep` object
#' @description Extract residuals from a fitted `cecl_dep` object.
#' @param object Object of class `cecl_dep`.
#' @param ... Additional arguments (not used).
#' @return List of matrices of residuals for each location.
#' @rdname residuals.cecl_dep
#' @export
#' @method residuals cecl_dep
residuals.cecl_dep <- \(object, ...) {
  stopifnot(inherits(object, "cecl_dep"))
  object$residual
}

#' @title Print summary of `cecl_dep` object
#' @description Print a summary of a fitted `cecl_dep` object.
#' @param x Object of class `cecl_dep`.
#' @param ... Additional arguments (not used).
#' @rdname print.cecl_dep
#' @export
#' @method print cecl_dep
print.cecl_dep <- \(x, ...) {
  stopifnot(inherits(x, "cecl_dep"))
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of locations:", length(x$dependence), "\n")
  vars <- sort(unique(unlist(lapply(
    x$dependence[[1]],
    colnames,
  ), use.names = FALSE)))
  cond_vars <- sort(names(x$dependence[[1]]))
  cat("Dependent Variables:", paste(vars, collapse = ", "), "\n")
  cat("Conditioning Variables:", paste(cond_vars, collapse = ", "), "\n")

  invisible(x)
}

#' @title Summary of `cecl_dep` object
#' @description Summarise a fitted `cecl_dep` object.
#' @param object Object of class `cecl_dep`.
#' @param ... Additional arguments (not used).
#' @return Data frame summarising dependence parameters for each location
#' and conditioned variable.
#' @rdname summary.cecl_dep
#' @export
#' @method summary cecl_dep
summary.cecl_dep <- \(object, ...) {
  stopifnot(inherits(object, "cecl_dep"))
  dep_params <- stats::coef(object)
  dep_params
}

#' @title Generic scatter plot function
#' @description Generic scatter plot function for different object classes.
#' @param x Object to plot.
#' @param ... Additional arguments passed to methods.
#' @return Plot of object.
#' @rdname plot_scatter
#' @keywords internal
plot_scatter <- \(x, ...) {
  UseMethod("plot_scatter")
}

#' @title Plot scatter plot from `cecl_dep` object
#' @description Plot scatter plot of dependence parameters from a fitted
#' `cecl_dep` object.
#' @param x Object of class `cecl_dep`.
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
#' @method plot_scatter cecl_dep
plot_scatter.cecl_dep <- \(
  x, var, cond_var, labels = NULL, type = c("ggplot", "plot"), ...
) {
  stopifnot(inherits(x, "cecl_dep"))
  type <- match.arg(type)

  a <- b <- name <- NULL

  # pull dependence parameters for all locations
  dep_params <- stats::coef(x)
  # pull for specific var/cond_var
  dep_params_spec <- dep_params[
    dep_params$var == var & dep_params$cond_var == cond_var,
  ]

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
      ggplot2::ggplot(ggplot2::aes(x = a, y = b)) +
      ggplot2::geom_point(...) +
      ggplot2::facet_wrap(~facet_lab) +
      cecl_theme() +
      ggplot2::labs(
        x = expression(a),
        y = expression(b),
      )

    # add labels if ggrepel is installed
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      plot <- plot +
        ggrepel::geom_text_repel(ggplot2::aes(label = name))
    } else {
      plot <- plot +
        ggplot2::geom_text(ggplot2::aes(label = name), vjust = -0.5)
    }

    return(plot)
  } else {
    plot(
      dep_params_spec$a,
      dep_params_spec$b,
      xlab = expression(alpha),
      ylab = expression(beta),
      main = paste0(var_lab, " | ", cond_var_lab),
      pch = 16,
      col = grDevices::rgb(0, 0, 0, 0.5),
      xlim = c(-1, 1),
      ylim = c(min(dep_params_spec$b) - 0.1, max(dep_params_spec$b) + 0.1),
      ...
    )

    graphics::text(
      dep_params_spec$a,
      dep_params_spec$b,
      labels = dep_params_spec$name,
      pos = 3
    )
  }
  return(invisible(NULL))
}

#' @title Plot residuals from `cecl_dep` object
#' @description Plot residuals from a fitted `cecl_dep` object.
#' @param obj Object of class `cecl_dep`.
#' @param loc Location name to plot residuals for.
#' @param var Conditioned variable name to plot residuals for.
#' @param cond_var Conditioning variable name to plot residuals against.
#' @param labels List mapping variable names to plot labels, e.g.,
#' `list("rain" = "Precipitation", "wind" = "Wind Speed")`, for use in axis
#' labels.
#' Default is `NULL`, which uses variable names as is.
#' @return ggplot object of residuals plot.
#' @rdname plot_resid
#' @keywords internal
plot_resid <- \(
  obj, loc, var, cond_var, labels = NULL, type = c("ggplot", "plot")
) {
  stopifnot(inherits(obj, "cecl_dep"))
  type <- match.arg(type)
  # pull specific residuals
  Z <- obj$residual[[loc]][[cond_var]][, var, drop = FALSE]
  # pull specific dependence quantile
  dqu <- obj$dependence[[loc]][[cond_var]]["dth", var, drop = TRUE]

  # For plotting, tidy up variable names
  cond_var_plt <- cond_var
  lhs_var_plt <- var
  if (!is.null(labels)) {
    cond_var_plt <- labels[[cond_var]]
    lhs_var_plt <- labels[[var]]
  }

  if (all(is.na(Z))) {
    return(NA)
  }

  n <- length(Z)

  # TODO Fix this!!! Should use dqu, but actually uses dth!!!
  # p <- seq(dqu, 1 - (1 / n), length = n)
  p <- seq(0.9, 1 - (1 / n), length = n)

  if (type == "ggplot") {
    plot <- data.frame(p, "resid" = Z) |>
      ggplot2::ggplot(ggplot2::aes(x = p, y = Z)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_smooth() +
      cecl_theme() +
      ggplot2::labs(
        x = paste0("F(", cond_var_plt, ")"),
        y = paste0("Z ", lhs_var_plt, " | ", cond_var_plt)
      )
    return(plot)
  } else {
    plot(
      p, Z,
      xlab = paste0("F(", cond_var_plt, ")"),
      ylab = paste0("Z ", lhs_var_plt, " | ", cond_var_plt),
      # TODO Do I want to use a different colour set?
      pch = 16, col = grDevices::rgb(0, 0, 0, 0.5)
    )
    graphics::lines(
      stats::loess.smooth(p, Z),
      col = "blue", lwd = 2
    )
  }
  return(invisible(NULL))
}

# plot quantiles of conditional expectation at single location (for single var)
plot_quantile <- \(
  obj,
  loc,
  var,
  cond_var,
  quantiles = seq(0.1, by = 0.2, len = 5),
  labels = NULL,
  type = c("ggplot", "plot")
) {
  type <- match.arg(type)
  stopifnot(inherits(obj, "cecl_dep"))

  x <- y <- NULL # to appease R CMD check

  # take out data for one location
  dep_fit_spec <- list(
    "residual" = obj$residual[[loc]][[cond_var]][, var, drop = FALSE],
    "dependence" = obj$dependence[[loc]][[cond_var]][, var, drop = FALSE],
    "transformed" = obj$transformed[[loc]][, c(var, cond_var), drop = FALSE]
  )
  if (any(vapply(dep_fit_spec, length, numeric(1)) == 0)) {
    stop("No data for specified loc, var, and cond_var in `cecl_dep` object")
  }

  n <- nrow(dep_fit_spec$residual)

  # dependence parameters for conditioning variables
  dep <- dep_fit_spec$dependence
  # dependence quantiles and thresholds (on laplace scale)
  dth <- dep["dth", ]
  # calculate dependence quantile by reversing dth for transformed data
  # TODO Maybe supply from outside???
  dqu <- stats::ecdf(dep_fit_spec$transformed[, cond_var, drop = TRUE])(dth)
  # Determine x-axis values to estimate CE quantiles at along conditioned var
  xmax <- max(dep_fit_spec$transformed[, cond_var])
  dif <- xmax - dth
  xlim <- c(dth - 0.1 * dif, dth + 1.5 * dif)

  # Upper limit of x-axis
  plim <- 1
  # CDF probabilities to plot at
  p <- seq(dqu, 1 - 1 / n, length = n)
  # take out largest point to avoid Inf in CDF transform
  len <- 501
  plotp <- seq(dqu, plim, len = len)[-len]
  # transform to Laplace scale; these will be x-values in plot
  plotx <- as.vector(dlaplace(plotp))

  # convert probs to Laplace scale (these are values to calculate CE line at)
  xq <- dlaplace(plotp)

  # pull dependence coefficients and quantiles of residuals
  co <- dep_fit_spec$dependence
  zq <- stats::quantile(dep_fit_spec$residual, quantiles)

  # calculates regression lines from quantiles of residuals
  yq <- sapply(zq, \(z, xq) {
    (co["a", ] * xq) + ((xq^co["b", ]) * z)
  }, xq = xq)

  # Previously transformed to original margins; now keep on Laplace scale
  ploty <- yq
  dth_spec <- dth

  # For plotting, tidy up variable names
  if (!is.null(labels)) {
    labels <- labels(c(cond_var, var))
  } else {
    labels <- c(cond_var, var)
  }

  if (type == "ggplot") {
    base_plot <- data.frame(dep_fit_spec$transformed) |>
      stats::setNames(c("x", "y")) |>
      ggplot2::ggplot(ggplot2::aes(x, y)) +
      ggplot2::geom_point() +
      cecl_theme() +
      # add vertical line at threshold
      ggplot2::geom_vline(xintercept = dth) +
      ggplot2::labs(x = labels[[1]], y = labels[[2]])

    # plot CE quantiles recursively
    add_line_ggplot <- \(p, ploty) {
      if (length(ploty) == 0) {
        p
      } else {
        add_line_ggplot(
          p +
            ggplot2::geom_line(
              data     = data.frame(x = plotx, y = ploty[, 1]),
              mapping  = ggplot2::aes(x = plotx, y = ploty[, 1]),
              linetype = 2,
              col      = "blue"
            ),
          ploty[, -1, drop = FALSE]
        )
      }
    }
    p <- add_line_ggplot(base_plot, ploty)
    return(p)
  } else {
    plot(
      dep_fit_spec$transformed[, cond_var, drop = TRUE],
      dep_fit_spec$transformed[, var, drop = TRUE],
      xlab = labels[[1]],
      ylab = labels[[2]],
      pch = 16,
      col = grDevices::rgb(0, 0, 0, 0.5)
    )
    graphics::abline(v = dth, lty = 2)

    # plot CE quantiles recursively
    add_line_plot <- \(ploty) {
      if (length(ploty) == 0) {
        NULL
      } else {
        graphics::lines(
          plotx,
          ploty[, 1],
          lty = 2,
          col = "blue"
        )
        add_line_plot(ploty[, -1, drop = FALSE])
      }
    }
    add_line_plot(ploty)
  }
  return(invisible(NULL))
}

#' @title Plot from `cecl_dep` object
#' @description Plot residuals or conditional quantiles from a fitted
#' `cecl_dep` object.
#' @param x Object of class `cecl_dep`.
#' @param which Character string specifying which plot to produce.
#' Either `"residual"` for residuals plot, `"quantile"` for conditional
#' quantiles plot, or `"scatter"` for dependence parameters scatter plot.
#' @param loc Location name to plot for.
#' @param var Conditioned variable name to plot for.
#' @param cond_var Conditioning variable name to plot against.
#' @param quantiles Numeric vector of quantiles to plot for conditional
#' quantiles plot. Default is `seq(0.1, by = 0.2, len = 5)`.
#' @param labels List mapping variable names to plot labels, e.g.,
#' `list("rain" = "Precipitation", "wind" = "Wind Speed")`, for use in axis
#' labels.
#' Default is `NULL`, which uses variable names as is.
#' @param ... Additional arguments passed to plotting functions.
#' @return ggplot object of specified plot.
#' @rdname plot.cecl_dep
#' @export
#' @method plot cecl_dep
plot.cecl_dep <- \(
  x,
  which = c("residual", "quantile", "scatter"),
  loc,
  var,
  cond_var,
  quantiles = seq(0.1, by = 0.2, len = 5),
  labels = NULL,
  ...
) {
  stopifnot(inherits(x, "cecl_dep"))
  which <- match.arg(which)

  if (which == "residual") {
    plot_resid(
      obj      = x,
      loc      = loc,
      var      = var,
      cond_var = cond_var,
      labels   = labels,
      type     = "plot"
    )
  } else if (which == "quantile") {
    plot_quantile(
      obj       = x,
      loc       = loc,
      var       = var,
      cond_var  = cond_var,
      quantiles = quantiles,
      labels    = labels,
      type      = "plot"
    )
  } else if (which == "scatter") {
    plot_scatter(
      x        = x,
      var      = var,
      cond_var = cond_var,
      labels   = labels,
      type     = "plot"
    )
  }
}

#' @title ggplot from `cecl_dep` object
#' @description Create ggplot of residuals or conditional quantiles from a
#' fitted `cecl_dep` object.
#' @param data Object of class `cecl_dep`.
#' @param mapping Not used.
#' @param which Character string specifying which plot to produce.
#' Either `"residual"` for residuals plot or `"quantile"` for conditional
#' quantiles plot.
#' @inheritParams plot.cecl_dep
#' @param ... Additional arguments passed to plotting functions.
#' @param environment Not used.
#' @return ggplot object of specified plot.
#' @rdname ggplot.cecl_dep
#' @export
#' @method ggplot cecl_dep
ggplot.cecl_dep <- \(
  data = NULL,
  mapping = ggplot2::aes(),
  which = c("residual", "quantile", "scatter"),
  loc,
  var,
  cond_var,
  quantiles = seq(0.1, by = 0.2, len = 5),
  labels = NULL,
  ...,
  environment = parent.frame()
) {
  stopifnot(inherits(data, "cecl_dep"))
  which <- match.arg(which)

  if (which == "residual") {
    p <- plot_resid(
      obj      = data,
      loc      = loc,
      var      = var,
      cond_var = cond_var,
      labels   = labels,
      type     = "ggplot",
      ...
    )
  } else if (which == "quantile") {
    p <- plot_quantile(
      obj       = data,
      loc       = loc,
      var       = var,
      cond_var  = cond_var,
      quantiles = quantiles,
      labels    = labels,
      type      = "ggplot",
      ...
    )
  } else if (which == "scatter") {
    p <- plot_scatter(
      x        = data,
      var      = var,
      cond_var = cond_var,
      labels   = labels,
      type     = "ggplot",
      ...
    )
  }
  return(p)
}

# TODO Be more detailed with description, make it clear that the user can
# supply their own method as long as it has a, b, m, s and dth values
#' @title Generic function to convert to `cecl_dep` object
#' @description Generic function to convert an object to class `cecl_dep`.
#' @param x Object to convert.
#' @param ... Additional arguments passed to methods.
#' @return Object of class `cecl_dep`.
#' @rdname as_cecl_dep
#' @export
as_cecl_dep <- \(x, ...) {
  UseMethod("as_cecl_dep")
}

#' @title Convert data.frame to `cecl_dep` object
#' @description Convert a data.frame of dependence parameters to class
#' `cecl_dep`.
#' @param x Data frame of dependence parameters.
#' @param name_col Name of column in `x` containing location names.
#' Default is `"name"`.
#' @param add_obj Additional named list object to add to returned
#' `cecl_dep` object. Default is `NULL`.
#' @param ... Additional arguments (not used).
#' @return Object of class `cecl_dep`.
#' @rdname as_cecl_dep
#' @export
#' @method as_cecl_dep data.frame
as_cecl_dep.data.frame <- \(
  x,
  name_col = "name",
  add_obj = NULL,
  ...
) {
  stopifnot(inherits(x, "data.frame"))
  stopifnot(
    "name_col, var and cond_var columns must be present in data.frame" =
      all(c(name_col, "var", "cond_var") %in% colnames(x))
  )
  stopifnot(
    all(c("a", "b", "m", "s", "dth") %in% colnames(x))
  )

  # first split into a list with one element per group/location
  x_lst <- x |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(name_col), \(y) factor(y, levels = unique(x[[name_col]]))
    )) |>
    dplyr::group_split(.data[[name_col]], .keep = FALSE)

  # next, convert each to a matrix of parameters
  vars <- unique(x$var)
  ret <- lapply(x_lst, \(y) {
    # create list objects for each variable
    params_list <- lapply(vars, \(v) {
      y_var <- y[y$var == v, ]
      # add dummy ll if not in data
      if (!"ll" %in% colnames(y_var)) {
        y_var$ll <- NA
      }
      params_mat <- t(as.matrix(
        y_var[, c("a", "b", "m", "s", "ll", "dth")]
      ))
      colnames(params_mat) <- y_var$cond_var
      params_mat
    })
    names(params_list) <- vars
    params_list
  })
  names(ret) <- unique(x[[name_col]])
  ret <- list("dependence" = ret)
  if (!is.null(add_obj)) {
    ret <- c(ret, add_obj)
  }
  class(ret) <- "cecl_dep"
  ret
}


#' @title Convert list to `cecl_dep` object
#' @description Convert a list of data.frames of dependence parameters to class
#' `cecl_dep`.
#' @param x List of data frames of dependence parameters.
#' @param name_col Name of column in each data.frame in `x` containing
#' location names. Default is `"name"`.
#' @param add_obj Additional named list object to add to returned
#' `cecl_dep` object. Default is `NULL`.
#' @param ... Additional arguments (not used).
#' @return Object of class `cecl_dep`.
#' @rdname as_cecl_dep
#' @export
#' @method as_cecl_dep list
as_cecl_dep.list <- \(
  x,
  name_col = "name",
  add_obj = NULL,
  ...
) {
  stopifnot(inherits(x, "list"))
  stopifnot(
    all(sapply(x, \(y) {
      inherits(y, "data.frame")
    }))
  )

  # if x is unnamed, assume name column is present in each data.frame
  if (is.null(names(x))) {
    x_df <- dplyr::bind_rows(x)
  } else {
    x_df <- dplyr::bind_rows(x, .id = name_col)
  }

  # call data.frame method
  as_cecl_dep.data.frame(
    x = x_df,
    name_col = name_col,
    add_obj = add_obj
  )
}
