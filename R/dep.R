#' @title Fit conditional extremes dependence model
#' @description Fit the conditional extremes dependence model of Heffernan
#' and Tawn (2004) to data transformed to Laplace scale.
#' @param obj Object of class `cecl_marg` containing transformed data.
#' @param cond_prob Numeric vector of quantiles to use for conditioning
#' variables.
#' @param vars Character vector of variable names to fit dependence model for.
#' If `NULL`, all variables are used.
#' @param cond_var Character vector of variable names to condition on.
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
cecl_dep <- \(
  obj,
  cond_prob,
  vars = NULL,
  cond_var = NULL,
  start = c("a" = 0.01, "b" = 0.01),
  ncores = 1,
  nruns = 1,
  aLow = -1,
  fixed_b = FALSE,
  fit_no_keef = FALSE
) {
  stopifnot(inherits(obj, "cecl_marg"))

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
  if (is.null(cond_var)) {
    cond_var <- vars
  }
  stopifnot("cond_var not in variables, check again" = cond_var %in% vars)

  # test_start <- \(start, marginal) {
  #   # check same locations
  #   test_loc <- length(start) == length(marginal)
  #   # check same variables
  #   test_var <- all(unlist(lapply(seq_along(start), \(i) {
  #     length(start[[i]]) == length(marginal[[i]]) &&
  #       all(names(start[[i]]) == names(marginal[[i]]))
  #   })))
  #   # check start values for each variable against each conditioning variable
  #   test_dim <- all(unlist(lapply(seq_along(start), \(i) {
  #     lapply(seq_along(start[[i]]), \(j) {
  #       all(colnames(start[[i]][[j]]) == names(marginal[[i]])[-j]) &&
  #         nrow(start[[i]][[j]]) == 2
  #     })
  #   })))
  #   stopifnot(
  #     "Start values not correct" = all(c(test_loc, test_var, test_dim))
  #   )
  # }

  # Test start values for dependence parameters (a and b)
  test_start <- \(start, cond_var, vars) {
    if (is.list(start)) {
      test_loc <- length(start) > 0
      test_var <- all(unlist(lapply(seq_along(start), \(i) {
        all(names(start[[i]]) %in% cond_var)
      })))
      test_params <- all(unlist(lapply(seq_along(start), \(i) {
        all(unlist(lapply(start[[i]], \(x) {
          all(rownames(x) == c("a", "b"))
        })))
      })))
      stopifnot(
        "Start values for dependence parameters not correct" =
          all(c(test_loc, test_var, test_params))
      )
    } else {
      stopifnot(
        "Start values for dependence parameters must be named vector" =
          is.numeric(start) && !is.null(names(start)) &&
            all(c("a", "b") %in% names(start))
      )
    }
  }


  # check if start values for dependence is a list
  is_list_start <- FALSE
  if (is.list(start)) {
    is_list_start <- TRUE
    test_start(start, cond_var, vars)
    # keep start values only for conditioned variables
    if (all(cond_var == vars) == FALSE) {
      start <- lapply(start, \(x) {
        x[names(x) %in% cond_var]
      })
    }
  }

  # fit dependence model to transformed data for each location
  dependence <- loop_fun(seq_along(marginal_trans), \(i) {
    start_spec <- start
    if (is_list_start) {
      start_spec <- start[[i]] # pull specific start vals, if required
    }
    # if multiple dependence quantiles are specified
    cond_prob_spec <- cond_prob
    if (length(cond_prob) == length(cond_var)) {
      cond_prob_spec <- cond_prob[i]
    }

    # fit dependence model
    o <- ce_optim(
      Y         = marginal_trans[[i]],
      dqu       = cond_prob_spec,
      cond_var  = cond_var,
      control   = list(maxit = 1e6),
      constrain = !fit_no_keef,
      aLow      = aLow,
      start     = start_spec,
      nruns     = nruns,
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
  ret$call <- match.call()
  ret$transformed <- obj$transformed

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
  class(ret) <- class(ret) <- c(
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
laplace_npll <- function(yex, ydep, a, b, constrain, v, aLow) {
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
Qpos <- function(param, yex, ydep, constrain, v, aLow) {
  a <- param[1]
  b <- param[2]

  res <- laplace_npll(yex, ydep, a, b, constrain, v, aLow)
  res$profLik
}

Qpos_fixed_b <- function(param, yex, ydep, constrain, v, aLow, b) {
  a <- param[1]
  res <- laplace_npll(yex, ydep, a, b, constrain, v, aLow)
  res$profLik
}

# fit CE for all pairs of variables
ce_optim <- \(
  Y,
  dqu,
  cond_var = NULL,
  start = c("a" = 0.01, "b" = 0.01),
  control = list(maxit = 1e6),
  constrain = TRUE,
  v = 10,
  aLow = -1,
  fixed_b = FALSE,
  nruns = 2
) {
  # check if start is a list (of start values for each location and variable)
  is_list_start <- is.list(start)

  # check that Y has names; if not give dummy names
  names_y <- colnames(Y)
  ncol_y <- ncol(Y)
  if (is.null(names_y)) {
    names_y <- paste0("var_", seq_len(ncol_y))
    colnames(Y) <- names_y
  }

  # check that dqu is a single value or vector
  if (is.null(cond_var)) {
    cond_var <- names_y
  }

  # optimise for a single variable vs another
  single_optim <- \(yex, ydep, start, dqu) {
    # threshold data
    thresh <- stats::quantile(yex, dqu)
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
      o_single$par <- c(o_single$par, NA, NA)
    }
    # add LL and numerical (Laplace scale) threshold, label appropriately
    o_single$par <- c(o_single$par, o_single$value, thresh[[1]])
    names(o_single$par) <- names(err_obj$params)
    # TODO Add checks afterwards on o
    list("resid" = matrix(Z), "params" = o_single$par)
  }

  # loop through variables, fit CE model against other variables
  ret <- lapply(seq_along(cond_var), \(i) {
    # can have different depenendence quantiles for each variable
    dqu_spec <- dqu
    if (length(dqu) > 1) {
      dqu_spec <- dqu[i]
    }

    # loop through conditioning variables
    o_yex <- lapply(seq_len(ncol_y - 1), \(j) {
      # conditioning variable (Y_{i}/LHS in CE model)
      yex <- Y[, which(colnames(Y) == cond_var[i])]
      # j'th conditioned variable (single vec in Y_{-i}/RHS of model)
      ydep <- Y[, -which(colnames(Y) == cond_var[i]), drop = FALSE][
        , j,
        drop = FALSE
      ]

      # extract specific start values
      start_spec <- start
      if (is_list_start) {
        start_spec <- start[[i]][, j, drop = TRUE]
      }
      o <- single_optim(yex, ydep, start_spec, dqu_spec)

      # check if optimisation failed
      if (all(is.na(o$params))) {
        return(o)
      }

      # perform multiple runs with prev estimates as start values, if desired
      if (nruns > 1) {
        for (i in seq_len(nruns) - 1) {
          start_spec <- o$params[1:2]
          o <- single_optim(yex, ydep, start_spec, dqu_spec)
        }
      }
      o
    })

    # join matrices from lists (each column will be for each conditioning var)
    o_yex <- list(
      "resid"  = do.call(cbind, lapply(o_yex, `[[`, "resid")),
      "params" = do.call(cbind, lapply(o_yex, `[[`, "params"))
    )

    if (!is.null(names_y)) {
      names_other <- names_y[names_y != cond_var[[i]]]
      colnames(o_yex$resid) <- names_other
      colnames(o_yex$params) <- names_other
    }
    o_yex
  })
  names(ret) <- cond_var
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
coef.cecl_dep <- \(object, ...) {
  stopifnot(inherits(object, "cecl_dep"))

  name <- var <- cond_var <- NULL # to appease R CMD check

  # extract dependence parameters
  dep_params <- object$dependence
  # convert to data.frame for easier viewing
  dep_params_df <- do.call(rbind, lapply(names(dep_params), \(loc) {
    params_loc <- as.data.frame(dep_params[[loc]])
    params_loc$parameter <- rownames(params_loc)
    params_loc$name <- loc
    rownames(params_loc) <- NULL

    # convert to wide to match coef.cecl_marg output
    params_loc_wide <- tidyr::pivot_longer(
      params_loc,
      cols = -c("name", "parameter"),
      names_to = "cond_var",
      values_to = "value"
    ) |>
      tidyr::pivot_wider(
        names_from = "parameter",
        values_from = "value"
      )

    # split column name, if required
    if (ncol(params_loc) > 4) {
      var_names <- stringr::str_split(params_loc_wide$cond_var, "\\.", n = 2)
      params_loc_wide$var <- vapply(
        var_names,
        \(x) x[[1]],
        character(1)
      )
      params_loc_wide$cond_var <- vapply(
        var_names,
        \(x) x[[2]],
        character(1)
      )
    } else {
      # opposite to cond_var
      params_loc_wide$var <- names(dep_params[[loc]])
    }

    params_loc_wide
  }))

  ret <- dplyr::relocate(
    dep_params_df, name, var, cond_var, dplyr::everything()
  )
  as.data.frame(ret)
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
  dep_params <- coef.cecl_dep(object)
  dep_params
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
  Z <- obj$residual[[loc]][[var]][, cond_var, drop = FALSE]
  # pull specific dependence quantile
  dqu <- obj$dependence[[loc]][[var]]["dth", cond_var, drop = TRUE]

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

  p <- seq(dqu, 1 - (1 / n), length = n)

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
    "residual" = obj$residual[[loc]][[var]][, cond_var, drop = FALSE],
    "dependence" = obj$dependence[[loc]][[var]][, cond_var, drop = FALSE],
    "transformed" = obj$transformed[[loc]][, c(var, cond_var), drop = FALSE]
  )

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
}

#' @title Plot from `cecl_dep` object
#' @description Plot residuals or conditional quantiles from a fitted
#' `cecl_dep` object.
#' @param x Object of class `cecl_dep`.
#' @param which Character string specifying which plot to produce.
#' Either `"residual"` for residuals plot or `"quantile"` for conditional
#' quantiles plot.
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
  which = c("residual", "quantile"),
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
  which = c("residual", "quantile"),
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
  }
  return(p)
}
