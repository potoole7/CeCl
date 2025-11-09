#### Testing of CeCl functions ####

# TODO Need to add functions for choosing thresholds (i.e. helping choice)
# TODO For thresholding, must allow either threshold for either variable, or
# a list for each variable of thresholds for each location
# Can't have the same thresholds for different variables!!!

#### Libs ####

devtools::load_all(".")
library(copula)
library(dplyr)

#### Metadata ####

cor <- c(0.7)
n_vars <- 2
n_locs <- 10
df_t <- 3
n <- 10000
stopifnot(n %% n_locs == 0)

cor_t3 <- c(0.3, 0.6, 0.9) # correlation values for 3 clusters
n_t3 <- 9000
n_locs_t3 <- 9
clust_mem <- rep(1:3, each = n_t3 / 3 / (n_t3 / n_locs_t3))

#### Generate Data ####

# Generate t copula data with student-t marginals
set.seed(123)
gen_dat <- function(cor, n_vars, df_t, n, n_locs, start_loc_n = 1) {
  cop_t <- copula::tCopula(param = cor, dim = n_vars, df = df_t, dispstr = "ex")
  u <- copula::rCopula(n, cop_t)
  data <- data.frame(apply(u, 2, qt, df = df_t))
  name_col <- rep(
    paste0("loc_", start_loc_n:(n_locs + start_loc_n - 1)),
    each = floor(n / n_locs)
  )
  # ensure name_col is same length as data; if not
  if (length(name_col) != nrow(data)) {
    name_col <- c(
      name_col, rep(last(name_col), nrow(data) - length(name_col))
    )
  }
  data$name <- name_col
  return(data)
}

df <- gen_dat(cor, n_vars, df_t, n, n_locs)
df3 <- gen_dat(cor, 3, df_t, n, n_locs) # for testing with 3 variables

# generate 3 clusters with different t-copula correlations
df_clust <- bind_rows(
  gen_dat(cor_t3[1], 3, df_t, n_t3 / 3, n_locs_t3 / 3, start_loc_n = 1),
  gen_dat(cor_t3[2], 3, df_t, n_t3 / 3, n_locs_t3 / 3, start_loc_n = (n_locs_t3 / 3) + 1),
  gen_dat(cor_t3[3], 3, df_t, n_t3 / 3, n_locs_t3 / 3, start_loc_n = 2 * (n_locs_t3 / 3) + 1)
)


#### Threshold ####

# TODO Anything else to test here??

# 1: Specific value
thresh_val <- cecl_marg(
  df,
  thresh_method = "value",
  # thresh_args = quantile(df$X1, 0.9), # TODO can be length 1 or length n
  thresh_args = apply(df[, c("X1", "X2")], 2, quantile, probs = 0.9),
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
  thresh_method = "qgam",
  # TODO Add more checking for these arguments
  thresh_args = list(
    f      = list("response ~ name", "~ name"),
    qu     = .9,
    jitter = TRUE
  ),
  thresh_only = TRUE,
  ret_obj = TRUE,
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
  ncores = 1,
  # ret_obj = FALSE
  ret_obj = TRUE
)

# fit for 3 variables
marg_ismev3 <- cecl_marg(
  df3,
  thresh_method = "quantile",
  thresh_args = 0.9,
  marg_method = "ismev",
  ncores = 1,
  # ret_obj = FALSE
  ret_obj = TRUE
)

# 3: evgam
marg_evgam <- cecl_marg(
  df,
  thresh_method = "quantile",
  thresh_args = 0.9,
  marg_method = "evgam",
  marg_args = list(f = list("excess ~ name", "~ name")),
  ncores = 1,
  # ret_obj = FALSE
  ret_obj = TRUE
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
coef(marg_ismev3)
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
# TODO Should return levels be on log scale? See ismev::gpd.diag
# TODO Should also provide global plots across all locations?
tryCatch(
  plot(marg_ecdf),
  error = function(e) message(e$message)
)
plot(marg_ismev, which = "pp", loc = "loc_1", var = "X1")
# TODO should QQ plot show more negative quantiles?
# TODO Also doesn't seem to match, investigate
# rnfit <- gpd.fit(df[df$name == "loc_1", ]$X1, quantile(df[df$name == "loc_1", ]$X1, 0.9))
# gpd.diag(rnfit)
plot(marg_ismev, which = "qq", loc = "loc_1", var = "X1")
plot(marg_ismev, which = "hist", loc = "loc_1", var = "X1")
plot(marg_ismev, which = "return", loc = "loc_1", var = "X1")

# ggplot method
ggplot(marg_ismev, which = "pp", loc = "loc_1", var = "X1")
ggplot(marg_ismev, which = "qq", loc = "loc_1", var = "X1")
ggplot(marg_ismev, which = "hist", loc = "loc_1", var = "X1")
ggplot(marg_ismev, which = "return", loc = "loc_1", var = "X1")

#### Dependence modelling ####

devtools::load_all()
devtools::document()

# TODO Also output marginal model from this
# TODO Add method for also fitting marginal model within dependence?
dep <- cecl_dep(
  obj = marg_ismev,
  cond_prob = 0.9 # TODO Could also have dependence value?? i.e. dth not dqu
)
dep3 <- cecl_dep(
  obj = marg_ismev3,
  cond_prob = 0.9
)

coef(dep)
coef(dep3)

print(dep)
print(dep3)

summary(dep)
summary(dep3)
summary.cecl_dep(dep)
summary.cecl_dep(dep3)

# TODO Add argument checks for loc, var, cond_var
plot(dep, which = "residual", var = "X1", cond_var = "X2", loc = "loc_1")
plot(dep, which = "quantile", var = "X1", cond_var = "X2", loc = "loc_1")
plot(dep, which = "scatter", var = "X1", cond_var = "X2")

ggplot(dep, which = "residual", var = "X1", cond_var = "X2", loc = "loc_1")
ggplot(dep, which = "quantile", var = "X1", cond_var = "X2", loc = "loc_1")
ggplot(dep, which = "scatter", var = "X1", cond_var = "X2")


#### Clustering ####

devtools::load_all()

# prep data by running marginal and dependence modelling
clust_marg <- cecl_marg(
  df_clust,
  thresh_method = "quantile",
  thresh_args = 0.9,
  marg_method = "ecdf",
  ncores = 1,
  ret_obj = TRUE
)
clust_dep <- cecl_dep(
  obj = clust_marg,
  cond_prob = 0.9
)

# calculate divergence matrix
# debugonce(cecl_dist)
dist <- cecl_dist(clust_dep, clust_marg, seed = 123)

# methods:
print(dist)
summary(dist)
summary(dist, var = "X3")

# TODO Is order of labels wrong??
# TODO Allow distance binning for colours
plot(dist, which = "image")
plot(dist, which = "scree")
plot(dist, which = "scree", var = "X2")

# TODO Use my theme in this plot (and all others!)
ggplot(dist, which = "image") # TODO Change colour scheme, don't like
ggplot(dist, which = "scree")
ggplot(dist, which = "scree", var = "X2")

# cluster
# TODO Add lambda (weights) argument for clustering
# method for pre-computed distance
set.seed(123)
clust1 <- cecl_clust(dist, k = 3, cluster_mem = clust_mem)
# method for dep and marg objects
set.seed(123)
clust2 <- cecl_clust(
  clust_dep,
  clust_marg,
  k = 3,
  cluster_mem = clust_mem,
  seed = 123
)

# check clustering results are the same
norm(as.matrix(clust1$dist_mat - clust2$dist_mat))

clust <- clust1

print(clust) # TODO Print which variables used
summary(clust) # TODO Can we somehow remove silhouette? Not applicable here..

# TODO Test this works with locs mixed up (i.e. not loc 1 -> loc 3 in clust 1)
plot(clust, which = "image")
plot(clust, which = "scatter", var = "X1", cond_var = "X2")

# TODO Optionally bin values (and in other heatmaps)
# TODO White squares across diagonal here, but not for other images?? Standardise!!
# TODO See if cluster colouring works when not ordering by cluster
ggplot(clust, which = "image")
ggplot(clust, which = "scatter", var = "X1", cond_var = "X2")
