#### ?? ####

#### Libs ####

devtools::load_all(".")
library(copula)
library(dplyr)

#### Metadata ####

cor <- c(0.5)
n_vars <- 2
n_locs <- 10
df <- 3
n <- 10000
stopifnot(n %% n_locs == 0)

#### Generate Data ####

# Generate t copula data with student-t marginals
cop_t <- copula::tCopula(cor, dim = n_vars, df = df, dispstr = "un")
u <- copula::rCopula(n, cop_t)
df <- data.frame(qt(u, df = df))
df$name <- rep(paste0("loc_", 1:n_locs), each = n / n_locs)

#### Threshold ####

# TODO Anything else to test??

# 1: Specific value
thresh_val <- cecl_marg(
  df,
  thresh_method = "value",
  thresh_args = quantile(df$X1, 0.9), # TODO can be length 1 or length n
  thresh_only = TRUE,
  ncores = 1
)

# 2. Quantile
thresh_q <- cecl_marg(
  df,
  thresh_method = "quantile",
  thresh_args = 0.9, # TODO can be length 1 or length n
  thresh_only = TRUE,
  ncores = 1
)

# 3. Regression
# debugonce(cecl_marg)
# debugonce(marg_thresh)
thresh_reg <- cecl_marg(
  df,
  thresh_method = "regression",
  # thresh_args = list(formula = X1 ~ name),
  # TODO Add more checking for these arguments
  thresh_args = list(
    f      = list("response ~ name", "~ name"),
    qu     = .9,
    jitter = TRUE
  ),
  thresh_only = TRUE,
  ncores = 1
)

# 4. None
data_input <- lapply(thresh_reg$data_thresh, function(x) bind_rows(unname(x)))
thresh_reg <- cecl_marg(
  thresh_method = "none",
  thresh_only = TRUE,
  ncores = 1
)



#### Fit marginal model ####

# 1: ECDF

# 2: ismev

# 3: evgam
