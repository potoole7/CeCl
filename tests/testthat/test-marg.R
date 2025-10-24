library(testthat)
library(dplyr)
library(copula)
devtools::load_all(".")

#### Setup data ####
set.seed(123)
n_vars <- 2
n_locs <- 10
df_t <- 3
n <- 10000
stopifnot(n %% n_locs == 0)

cop_t <- tCopula(0.7, dim = n_vars, df = df_t, dispstr = "un")
u <- rCopula(n, cop_t)
df <- data.frame(qt(u, df = df_t))
df$name <- rep(paste0("loc_", 1:n_locs), each = n / n_locs)

vars <- c("X1", "X2")

#### Thresholding tests ####
test_that("cecl_marg thresholds correctly", {
  # Value threshold
  val_thresh <- cecl_marg(df,
    thresh_method = "value",
    thresh_args = quantile(df$X1, 0.9),
    thresh_only = TRUE, ncores = 1
  )
  expect_s3_class(val_thresh, "cecl_thresh")
  expect_true(all(names(val_thresh$data_thresh) %in% vars))

  # Quantile threshold
  q_thresh <- cecl_marg(df,
    thresh_method = "quantile",
    thresh_args = 0.9,
    thresh_only = TRUE
  )
  expect_s3_class(q_thresh, "cecl_thresh")

  # Regression threshold
  reg_thresh <- cecl_marg(df,
    thresh_method = "regression",
    # thresh_args = list(f = list("X1 ~ name", "~ name"),
    thresh_args = list(
      f = list("response ~ name", "~ name"),
      qu = 0.9,
      jitter = TRUE
    ),
    thresh_only = TRUE
  )
  expect_s3_class(reg_thresh, "cecl_thresh")
})

#### Marginal model tests ####
test_that("cecl_marg fits different marginal methods", {
  # ECDF
  m_ecdf <- cecl_marg(df,
    thresh_method = "quantile", thresh_args = 0.9,
    marg_method = "ecdf"
  )
  expect_s3_class(m_ecdf, "cecl_marg")
  expect_null(m_ecdf$marginal)

  # ISMEV
  m_ismev <- cecl_marg(df,
    thresh_method = "quantile", thresh_args = 0.9,
    marg_method = "ismev"
  )
  expect_s3_class(m_ismev, "cecl_marg")
  expect_true(!is.null(m_ismev$marginal))
  expect_true(all(names(m_ismev$marginal) %in% paste0("loc_", 1:n_locs)))

  # EVGAM
  m_evgam <- cecl_marg(df,
    thresh_method = "quantile", thresh_args = 0.9,
    marg_method = "evgam",
    # marg_args = list(f = list("X1 ~ name", "~ name"))
    marg_args = list(f = list("excess ~ name", "~ name")),
  )
  expect_s3_class(m_evgam, "cecl_marg")
  expect_true(!is.null(m_evgam$evgam_fit))
})

#### S3 methods tests ####
test_that("cecl_marg S3 methods work", {
  m_ismev <- cecl_marg(df,
    thresh_method = "quantile", thresh_args = 0.9,
    marg_method = "ismev"
  )

  expect_s3_class(coef(m_ismev), "data.frame")
  expect_output(print(m_ismev))
  expect_output(summary(m_ismev, n = 5))
})

#### Base R plotting tests ####
test_that("cecl_marg base plotting works", {
  m_ismev <- cecl_marg(df,
    thresh_method = "quantile", thresh_args = 0.9,
    marg_method = "ismev"
  )

  pdf(NULL)
  expect_error(plot(m_ismev, which = "qq", loc = "loc_1", var = "X1"), NA)
  expect_error(plot(m_ismev, which = "pp", loc = "loc_1", var = "X1"), NA)
  expect_error(plot(m_ismev, which = "hist", loc = "loc_1", var = "X1"), NA)
  expect_error(plot(m_ismev, which = "return", loc = "loc_1", var = "X1"), NA)
  dev.off()
})

#### ggplot tests ####
test_that("cecl_marg ggplot works", {
  m_ismev <- cecl_marg(df,
    thresh_method = "quantile", thresh_args = 0.9,
    marg_method = "ismev"
  )

  # ggplot objects
  pdf(NULL)
  expect_s3_class(ggplot(m_ismev, which = "qq", loc = "loc_1", var = "X1"), "gg")
  expect_s3_class(ggplot(m_ismev, which = "pp", loc = "loc_1", var = "X1"), "gg")
  expect_s3_class(ggplot(m_ismev, which = "hist", loc = "loc_1", var = "X1"), "gg")
  expect_s3_class(ggplot(m_ismev, which = "return", loc = "loc_1", var = "X1"), "gg")
  dev.off()
})
