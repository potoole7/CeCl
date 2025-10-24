#### ?? ####

# TODO Need to add functions for choosing thresholds

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
set.seed(123)
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
thresh_reg <- cecl_marg(
  df,
  thresh_method = "regression",
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
# TODO Fix
# data_input <- lapply(thresh_reg$data_thresh, function(x) bind_rows(unname(x)))
# thresh_reg <- cecl_marg(
#   thresh_method = "none",
#   thresh_only = TRUE,
#   ncores = 1
# )

#### Fit marginal model ####

# 1: ECDF
marg_ecdf <- cecl_marg(
  df,
  thresh_method = "quantile",
  thresh_args = 0.9,
  marg_method = "ecdf",
  ncores = 1
)

# 2: ismev
marg_ismev <- cecl_marg(
  df,
  thresh_method = "quantile",
  thresh_args = 0.9,
  marg_method = "ismev",
  ncores = 1
)

# 3: evgam
marg_evgam <- cecl_marg(
  df,
  thresh_method = "quantile",
  thresh_args = 0.9,
  marg_method = "evgam",
  marg_args = list(f = list("excess ~ name", "~ name")),
  ncores = 1
)

#### Marginal methods ####

devtools::document()
devtools::load_all(".")

# coef method
tryCatch(
  coef(marg_ecdf),
  error = function(e) message(e$message)
)
coef(marg_ismev)
coef(marg_evgam)

# print method
print(marg_ecdf)
print(marg_ismev)
print(marg_evgam)

# summary method
tryCatch(
  summary(marg_ecdf),
  error = function(e) message(e$message)
)
summary(marg_ismev, n = 10)
summary(marg_evgam)

# plot method
tryCatch(
  plot(marg_ecdf),
  error = function(e) message(e$message)
)
plot(marg_ismev, which = "pp", loc = "loc_1", var = "X1")
plot(marg_ismev, which = "qq", loc = "loc_1", var = "X1")
plot(marg_ismev, which = "hist", loc = "loc_1", var = "X1")
plot(marg_ismev, which = "return", loc = "loc_1", var = "X1")

# ggplot method
ggplot(marg_ismev, which = "pp", loc = "loc_1", var = "X1")
# TODO investigate why this is so bad! Lol
ggplot(marg_ismev, which = "qq", loc = "loc_1", var = "X1")
ggplot(marg_ismev, which = "hist", loc = "loc_1", var = "X1")
ggplot(marg_ismev, which = "return", loc = "loc_1", var = "X1")
