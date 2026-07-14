test_that("fit_dynamic_model returns a dynamic_fit with the expected structure", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 1)
  fit <- fit_dynamic_model(sim$y, family = "poisson",
                           nsave = 300, nburn = 150, seed = 1)
  expect_s3_class(fit, "dynamic_fit")
  expect_named(fit, c("draws", "data", "spec"))
  expect_equal(fit$spec$family, "poisson")
  expect_equal(fit$spec$innovations, "gaussian")
  expect_equal(fit$data$n, 60)
  # posterior draw matrices have the requested number of rows
  expect_equal(nrow(fit$draws$z), 300)
  # latent state vector is the leading pad + n observations (n + 1 columns)
  expect_equal(ncol(fit$draws$z), 61)
  expect_equal(ncol(fit$draws$fitted), 60)
  # draws holds only posterior draws -- fixed inputs live elsewhere
  expect_null(fit$draws$offset)
  expect_equal(fit$data$offset, rep(0, 60))
})

test_that("Poisson model recovers the latent rate reasonably", {
  sim <- simulate_dynamic_poisson(n = 120, sigma = 0.15, log_rate0 = 3, seed = 11)
  fit <- fit_dynamic_model(sim$y, family = "poisson",
                           nsave = 750, nburn = 375, seed = 11)
  fitted_mean <- colMeans(fit$draws$fitted)
  expect_gt(cor(fitted_mean, sim$rate), 0.8)
})

test_that("all innovation structures run and summarise", {
  sim <- simulate_dynamic_poisson(n = 50, sigma = 0.2, log_rate0 = 2, seed = 5)
  for (innov in c("gaussian", "t", "mixture")) {
    fit <- fit_dynamic_model(sim$y, family = "poisson", innovations = innov,
                             nsave = 225, nburn = 112, seed = 5)
    expect_s3_class(fit, "dynamic_fit")
    expect_equal(fit$spec$innovations, innov)
    s <- summary(fit)
    expect_s3_class(s, "summary.dynamic_fit")
  }
})

test_that("the t innovation samples degrees of freedom", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 6)
  fit <- fit_dynamic_model(sim$y, family = "poisson", innovations = "t",
                           nsave = 300, nburn = 150, seed = 6)
  expect_true(!is.null(fit$draws$nu))
  expect_true(all(fit$draws$nu > 0))
})

test_that("invalid count input is rejected", {
  expect_error(fit_dynamic_model(c(1, 2, -3), family = "poisson"))
  expect_error(fit_dynamic_model(c(1.5, 2, 3), family = "poisson"))
  expect_error(fit_dynamic_model(c(1, 2), family = "poisson"))  # n < 3
})

test_that("nburn = 0 is allowed", {
  sim <- simulate_dynamic_poisson(n = 40, sigma = 0.2, log_rate0 = 2, seed = 1)
  fit <- fit_dynamic_model(sim$y, nsave = 75, nburn = 0, seed = 1)
  expect_s3_class(fit, "dynamic_fit")
  expect_equal(nrow(fit$draws$z), 75)
})

test_that("sv innovations run when stochvol is available", {
  skip_if_not_installed("stochvol")
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 1)
  fit <- fit_dynamic_model(sim$y, family = "poisson", innovations = "sv",
                           nsave = 150, nburn = 75, seed = 1)
  expect_s3_class(fit, "dynamic_fit")
  expect_equal(fit$spec$innovations, "sv")
  expect_true(all(is.finite(fit$draws$innov_var)))
})

test_that("the mixture innovation runs with three components", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 13)
  fit <- fit_dynamic_model(sim$y, family = "poisson", innovations = "mixture",
                           prior = dynamic_prior(mix_components = 3),
                           nsave = 225, nburn = 112, seed = 13)
  expect_s3_class(fit, "dynamic_fit")
  expect_equal(fit$spec$prior$mix_components, 3L)
  expect_true(all(is.finite(fit$draws$innov_var)))
  expect_true(all(fit$draws$innov_var > 0))
})

test_that("non-finite offsets are rejected", {
  sim <- simulate_dynamic_poisson(n = 20, sigma = 0.2, log_rate0 = 2, seed = 1)
  expect_error(fit_dynamic_model(sim$y, offset = NA_real_), "finite")
  expect_error(fit_dynamic_model(sim$y, offset = Inf), "finite")
  expect_error(
    fit_dynamic_model(sim$y, horizon = 3, forecast_offset = "a"), "finite"
  )
})
