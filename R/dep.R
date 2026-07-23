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
  } else if (is.data.frame(start)) {
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
    
    # <<<<<<< Updated upstream
    # =======
    #     # filter to only relevant rows of start data, if cond_var specified
    #     start <- start |>
    #       dplyr::filter(cond_var %in% cond_vars)
    #
    #     if (nrow(start) == 0) {
    #       stop(
    #         "No rows of `start` data.frame match specified `cond_vars`, check that",
    #         " start has rows matching `cond_vars` or adjust `cond_vars` argument"
    #       )
    #     }
    # >>>>>>> Stashed changes
    
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
    # if location-specific dependence quantiles are specified
    cond_prob_spec <- cond_prob
    if (length(cond_prob) == length(locs_keep)) {
      cond_prob_spec <- cond_prob[[i]]
    }
    
    # same for cond_val
    cond_val_spec <- cond_val
    if (length(cond_val) == length(locs_keep)) {
      cond_val_spec <- cond_val[[i]]
    }
    
    # pull start values
    start_spec <- start
    if (is_df_start) {
      start_spec <- start |>
        dplyr::filter(
          name == locs_keep[i],
          var %in% vars,
          cond_var %in% cond_vars
        ) |>
        dplyr::select(a, b, var, cond_var)
      
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
    
    # <<<<<<< Updated upstream
    #     o <- ce_optim(
    #       Y         = marginal_trans[[i]],
    #       dqu       = cond_prob_spec,
    #       dth       = cond_val_spec,
    #       cond_var  = cond_var,
    # =======
    #     o <- ce_optim(
    #       Y         = marginal_trans[[i]],
    #       # dqu       = cond_prob,
    #       # dth       = cond_val,
    #       cond_vars = cond_vars,
    # >>>>>>> Stashed changes
    
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
    } else if (constrain) {
      zpos <- range(ydep - yex)
      z <- range((ydep - yex * a) / (yex^b))
      zneg <- range(ydep + yex)
      
      if (!constraints_satisfied(a, b, z, zpos, zneg, v)) {
        return(BigNumber)
      }
    }
  }
  res
}

# function to evaluate (negative) profile (log) likelihood and optimise over
laplace_npll <- \(yex, ydep, a, b, constrain, v, aLow) {
  Z <- (ydep - yex * a) / (yex^b)
  stopifnot(
    "NaNs in Z, conditional quantile may be negative" = all(!is.nan(Z))
  )
  
  m <- mean(Z)
  s <- stats::sd(Z)
  
  res <- laplace_nll(
    yex, ydep, a, b,
    m = m, s = s, constrain, v, aLow = aLow
  )
  list(profLik = res, m = m, s = s)
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
  # must specify either dqu or dth
  if (is.null(dqu) && is.null(dth) ||
      sum(!is.null(dqu), !is.null(dth)) > 1
  ) {
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
  
  # <<<<<<< Updated upstream
  #   # check that dqu is a single value or vector
  #   if (is.null(cond_var)) {
  #     cond_var <- names_y
  # =======
  #   var <- names_y
  #   if (is.null(cond_vars)) {
  #     cond_vars <- names_y # by default, condition on all variables
  # >>>>>>> Stashed changes
  #   }
  
  # conditioning variables default to all variables
  if (is.null(cond_vars)) {
    cond_vars <- names_y
  }
  
  if (!all(cond_vars %in% names_y)) {
    stop(
      "All `cond_vars` must be columns of `Y`.",
      call. = FALSE
    )
  }
  
  if (anyDuplicated(cond_vars)) {
    stop(
      "`cond_vars` must not contain duplicates.",
      call. = FALSE
    )
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
    # threshold data
    thresh <- dth
    if (is.null(dth)) {
      thresh <- stats::quantile(yex, dqu)
    }
    wch <- yex > thresh
    
    err_obj <- list(
      "resid" = matrix(NA, nrow = max(sum(wch), 1)),
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
    
    if (fixed_b == TRUE) {
      b <- start[[2]]
      start <- start[1]
      o_single <- try(stats::optim(
        par       = start,
        fn        = Qpos_fixed_b,
        method    = "Brent",
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
    
    if (all(o_single$par[1:2] == start) && all(start == 0.01)) {
      message("No change from starting values, optimisation failed")
      return(err_obj)
    }
    
    if (all(!is.na(o_single$par))) {
      Z <- (ydep[wch] - yex[wch] * o_single$par[1]) /
        (yex[wch]^o_single$par[2])
      o_single$par <- c(o_single$par[1:2], mean(Z), stats::sd(Z))
    } else {
      o_single$par <- c(o_single$par, NA, NA)
    }
    
    o_single$par <- c(o_single$par, o_single$value, thresh[[1]])
    names(o_single$par) <- names(err_obj$params)
    
    list("resid" = matrix(Z), "params" = o_single$par)
  }
  
  # <<<<<<< Updated upstream
  #   # loop through variables, fit CE model against other variables
  #   ret <- lapply(seq_along(cond_var), \(i) {
  # =======
  #   # loop through conditioning variables, fit CE model against other variables
  #   # TODO Can we fit models with multiple conditioning and/or conditioned vars?
  #   ret <- lapply(seq_along(cond_vars), \(i) {
  #     # ret <- lapply(seq_along(var), \(i) {
  # >>>>>>> Stashed changes
  
  # loop through RHS conditioning variables
  ret <- lapply(seq_along(cond_vars), \(i)) {
    cond_var_i <- cond_vars[[i]]
    which_cond <- match(cond_var_i, names_y)
    names_conditioned <- names_y[-which_cond]
    
    # can have different dependence quantiles/values for each conditioning var
    dqu_spec <- dqu
    if (length(dqu) > 1L) {
      dqu_spec <- dqu[[i]]
    }
    
    dth_spec <- dth
    if (length(dth) > 1L) {
      dth_spec <- dth[[i]]
    }
    
    # <<<<<<< Updated upstream
    #     # loop through conditioning/dependent variables
    #     # TODO Very hard to follow which variable is which! change `cond_var` name
    #     o_yex <- lapply(seq_len(ncol_y - 1), \(j) {
    #       # conditioning variable (Y_{i}/LHS in CE model)
    #       which_cond <- which(colnames(Y) == cond_var[i])
    #       yex <- Y[, which_cond, drop = TRUE]
    #       # j'th conditioned variable (single vec in Y_{-i}/RHS of model)
    #       ydep <- Y[, -which_cond, drop = FALSE][
    #         , j,
    #         drop = FALSE
    #       ]
    # =======
    #     # i'th conditioning variable (single vec in Y_{-i}/RHS of model)
    #     which_cond <- which(names_y == cond_vars[[i]])
    #     ydep <- Y[, which_cond, drop = FALSE]
    #
    #     # Loop through conditioned variables
    #     # TODO Very hard to follow which variable is which! change `cond_vars` name
    #     # TODO Does this logic still follow for > 2 variables?
    #     o_yex <- lapply(seq_len(ncol_y - 1), \(j) {
    #       # conditioned variable (Y_{j}/LHS in CE model)
    #       # which_cond <- which(names_y == cond_vars[i])
    #       # which_var <- which(names_y == var[i])
    #       # yex <- Y[, which_cond, drop = TRUE]
    #       # yex <- Y[, which_var, drop = TRUE]
    #       yex <- Y[, -which_cond, drop = FALSE][, j, drop = TRUE]
    #
    #       # j'th conditioning variable (single vec in Y_{-i}/RHS of model)
    #       # ydep <- Y[, -which_cond, drop = FALSE][
    #       #   , j,
    #       #   drop = FALSE
    #       # ]
    #       # ydep <- Y[, -which_var, drop = FALSE][
    #       #   , j,
    #       #   drop = FALSE
    #       # ]
    # >>>>>>> Stashed changes
    
    # RHS of the CE equation: conditioning variable
    yex <- Y[, which_cond, drop = TRUE]
    
    # loop through LHS conditioned variables
    o_yex <- lapply(seq_along(names_conditioned), \(j) {
      conditioned_var_j <- names_conditioned[[j]]
      which_conditioned <- match(conditioned_var_j, names_y)
      
      # LHS of the CE equation: conditioned variable
      ydep <- Y[, which_conditioned, drop = TRUE]
      
      # extract specific start values (if data.frame)
      start_spec <- start
      if (is.data.frame(start_spec)) {
        # <<<<<<< Updated upstream
        #         start_spec <- start_spec |>
        #           dplyr::filter(
        #             var == !!cond_var[[which_cond]],
        #             cond_var == !!cond_var[-which_cond][j]
        # =======
        #         start_spec <- start_spec |>
        #           dplyr::filter(
        #             # var == !!cond_vars[[which_cond]],
        #             # cond_var == !!cond_vars[-which_cond][j]
        #             var == names_y[-which_cond][j],
        #             cond_var == names_y[which_cond]
        # >>>>>>> Stashed changes
        #           ) |>
        #           dplyr::select(a, b) |>
        #           unlist()
        
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
      
      # perform multiple runs with previous estimates as start values
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
    
    # join matrices from lists
    o_yex <- list(
      "resid" = do.call(cbind, lapply(o_yex, `[[`, "resid")),
      "params" = do.call(cbind, lapply(o_yex, `[[`, "params"))
    )
    
    # <<<<<<< Updated upstream
    #     if (!is.null(names_y)) {
    #       names_other <- names_y[names_y != cond_var[[i]]]
    #       colnames(o_yex$resid) <- names_other
    #       colnames(o_yex$params) <- names_other
    # =======
    #     if (!is.null(names_y)) {
    #       names_ex <- names_y[names_y != cond_vars[[i]]]
    #       # names_ex <- names_y[names_y != var[[i]]]
    #       colnames(o_yex$resid) <- names_ex
    #       colnames(o_yex$params) <- names_ex
    # >>>>>>> Stashed changes
    #     }
    
    # Columns are named after the LHS conditioned variables.
    colnames(o_yex$resid) <- names_conditioned
    colnames(o_yex$params) <- names_conditioned
    
    o_yex
  })

# Outer list elements are named after RHS conditioning variables.
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
  
  name <- var <- cond_var <- NULL
  
  dep_params <- object$dependence
  dep_params_df <- do.call(rbind, lapply(names(dep_params), \(loc) {
    params_loc <- as.data.frame(dep_params[[loc]])
    params_loc$parameter <- rownames(params_loc)
    params_loc$name <- loc
    rownames(params_loc) <- NULL
    
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
      params_loc_wide$var <- names(dep_params[[loc]])
    }
    
    params_loc_wide
  }))
  
  ret <- as.data.frame(dplyr::relocate(
    dep_params_df, name, var, cond_var, dplyr::everything()
  ))
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
  
  # <<<<<<< Updated upstream
  #   n_vars <- colnames(x$transformed[[1]])
  #   # fail safe for user-specified dependence object with no trans data
  #   if (is.null(n_vars)) {
  #     n_vars <- names(x$dependence[[1]])
  #   }
  #   cat(
  #     "Variables:", paste(n_vars, collapse = ", "), "\n"
  #   )
  # =======
  #   vars <- sort(unname(colnames(x$dependence[[1]][[1]])))
  #   cond_vars <- sort(names(x$dependence[[1]]))
  #   cat("Dependent Variables:", paste(vars, collapse = ", "), "\n")
  #   cat("Conditioning Variables:", paste(cond_vars, collapse = ", "), "\n")
  # >>>>>>> Stashed changes
  
  vars <- sort(unique(unlist(lapply(
    x$dependence[[1]],
    colnames,
    use.names = FALSE
  ))))
  cond_vars <- sort(names(x$dependence[[1]]))
  
  cat("Dependent Variables:", paste(vars, collapse = ", "), "\n")
  cat("Conditioning Variables:", paste(cond_vars, collapse = ", "), "\n")
  
  invisible(x)
}