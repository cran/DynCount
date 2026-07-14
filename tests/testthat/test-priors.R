test_that("dynamic_prior has documented defaults", {
  p <- dynamic_prior()
  expect_s3_class(p, "dynamic_prior")
  expect_equal(p$var_shape, 2.5)
  expect_equal(p$var_rate, 0.5)
  expect_equal(p$df_min, 3)
  expect_equal(p$df_mean_excess, 6)
  expect_equal(p$mix_components, 2L)
  expect_equal(p$mix_concentration, 1)
  expect_equal(p$zi_open_a, 1)
  expect_equal(p$zi_open_b, 1)
  expect_equal(p$ar_rho_mean, 0)
  expect_equal(p$ar_rho_sd, 1)
  expect_equal(p$mu_mean, 0)
  expect_equal(p$mu_sd, 1)
  expect_equal(p$init_mean, 0)
  expect_equal(p$init_var, 100)
})

test_that("dynamic_prior validates the rho prior", {
  expect_error(dynamic_prior(ar_rho_sd = 0))
  expect_error(dynamic_prior(ar_rho_sd = -1))
  expect_error(dynamic_prior(ar_rho_mean = Inf))
})

test_that("dynamic_prior validates its arguments", {
  expect_error(dynamic_prior(var_shape = -1))
  expect_error(dynamic_prior(var_rate = 0))
  expect_error(dynamic_prior(df_min = 0))
  expect_error(dynamic_prior(mix_components = 0))
  expect_error(dynamic_prior(mix_components = 2.7))
  expect_error(dynamic_prior(zi_open_a = -2))
})

test_that("a custom prior flows through to the fit spec", {
  sim <- simulate_dynamic_poisson(n = 40, sigma = 0.2, seed = 1)
  pr <- dynamic_prior(var_shape = 5, var_rate = 1)
  fit <- fit_dynamic_model(sim$y, family = "poisson", prior = pr,
                           nsave = 150, nburn = 75, seed = 1)
  expect_equal(fit$spec$prior$var_shape, 5)
  expect_equal(fit$spec$prior$var_rate, 1)
})

test_that("a tighter innovation-variance prior shrinks the posterior", {
  sim <- simulate_dynamic_poisson(n = 80, sigma = 0.05, log_rate0 = 2, seed = 21)
  loose <- fit_dynamic_model(sim$y, family = "poisson",
                             prior = dynamic_prior(var_shape = 2.5, var_rate = 0.5),
                             nsave = 450, nburn = 225, seed = 21)
  tight <- fit_dynamic_model(sim$y, family = "poisson",
                             prior = dynamic_prior(var_shape = 20, var_rate = 0.05),
                             nsave = 450, nburn = 225, seed = 21)
  expect_lt(mean(tight$draws$innov_var), mean(loose$draws$innov_var) * 1.5)
})
