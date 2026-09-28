test_that("the stored fit-time forecast is what forecast() returns by default", {
  sim <- simulate_dynamic_poisson(n = 50, sigma = 0.2, log_rate0 = 2, seed = 1)
  fit <- fit_dynamic_model(sim$y, nsave = 200, nburn = 100, horizon = 4, seed = 1)
  fc <- forecast(fit)
  expect_identical(fc$draws, fit$draws$forecast_y)
  expect_identical(fc$latent_draws, fit$draws$forecast_z)
  # a new horizon triggers a fresh forward simulation
  fc8 <- forecast(fit, horizon = 8, seed = 3)
  expect_equal(dim(fc8$draws), c(200, 8))
})

test_that("post-hoc forecasts are reproducible with a seed", {
  sim <- simulate_dynamic_poisson(n = 50, sigma = 0.2, log_rate0 = 2, seed = 2)
  fit <- fit_dynamic_model(sim$y, nsave = 150, nburn = 100, seed = 2)
  a <- forecast(fit, horizon = 6, seed = 10)
  b <- forecast(fit, horizon = 6, seed = 10)
  expect_identical(a$draws, b$draws)
})

test_that("Gaussian random-walk forecasts have the implied mean and variance", {
  # given a draw, z_{n+H} = z_n + N(0, H * sigma^2), so across draws
  # E[z_{n+H}] = E[z_n] and Var[z_{n+H}] = Var[z_n] + H * E[sigma^2]
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.15, log_rate0 = 3, seed = 4)
  fit <- fit_dynamic_model(sim$y, nsave = 3000, nburn = 500, seed = 4)
  H <- 10
  fc <- forecast(fit, horizon = H, seed = 5)
  zn <- fit$draws$z[, 60]
  zf <- fc$latent_draws[, H]
  expect_lt(abs(mean(zf) - mean(zn)), 4 * sqrt(var(zf) / 3000))
  v_theory <- var(zn) + H * mean(fit$draws$scale)
  expect_lt(abs(var(zf) / v_theory - 1), 0.1)
})

test_that("AR(1) and drift forecasts follow the state equation", {
  sim <- simulate_dynamic_poisson(n = 80, sigma = 0.1, log_rate0 = 2, mu = 0.05, seed = 6)
  fit <- fit_dynamic_model(sim$y, include_mu = TRUE, nsave = 1000, nburn = 300, seed = 6)
  fc <- forecast(fit, horizon = 20, seed = 7)
  # expected latent path grows with the drift: E[z_{n+H}] = E[z_n + H mu]
  expect_lt(abs(mean(fc$latent_draws[, 20]) -
                  mean(fit$draws$z[, 80] + 20 * fit$draws$mu)), 0.1)
})

test_that("post-hoc forecasts run for every innovation structure", {
  sim <- simulate_dynamic_poisson(n = 50, sigma = 0.2, log_rate0 = 2, seed = 8)
  innov <- c("gaussian", "t", "mixture")
  if (requireNamespace("stochvol", quietly = TRUE)) innov <- c(innov, "sv")
  for (inn in innov) {
    fit <- fit_dynamic_model(sim$y, innovations = inn, nsave = 100, nburn = 100, seed = 8)
    fc <- forecast(fit, horizon = 3, seed = 1)
    expect_true(all(is.finite(fc$latent_draws)))
    expect_true(all(fc$draws >= 0))
  }
})

test_that("forecast inputs are validated and defaulted", {
  sim <- simulate_dynamic_binomial(n = 40, sigma = 0.1, trials = 30, seed = 9)
  fit <- fit_dynamic_model(sim$y, family = "binomial", trials = 30,
                           nsave = 100, nburn = 50, seed = 9)
  fc <- forecast(fit, horizon = 3, forecast_trials = c(10, 20, 30), seed = 1)
  expect_true(all(sweep(fc$draws, 2, c(10, 20, 30), "<=")))
  expect_equal(forecast(fit, horizon = 2, seed = 1)$forecast_trials, c(30, 30))
  expect_error(forecast(fit, horizon = 3, forecast_trials = c(1, 2)), "length")
  expect_error(forecast(fit, horizon = 0), "horizon")
  # a model with an offset warns when no forecast offset is given
  simp <- simulate_dynamic_poisson(n = 40, sigma = 0.1, log_rate0 = -2,
                                   offset = log(100), seed = 9)
  fitp <- fit_dynamic_model(simp$y, offset = log(100), nsave = 100, nburn = 50, seed = 9)
  expect_warning(forecast(fitp, horizon = 2), "forecast_offset")
  expect_silent(forecast(fitp, horizon = 2, forecast_offset = log(100)))
})

test_that("multinomial forecast totals default to the last non-zero row total", {
  sim <- simulate_dynamic_multinomial(n = 30, sigma = 0.1, trials = 60, seed = 10)
  y <- sim$y; y[30, ] <- 0; y[29, ] <- c(10, 20, 30)
  fit <- fit_dynamic_model(y, family = "multinomial", nsave = 50, nburn = 30, seed = 10)
  fc <- forecast(fit, horizon = 2, seed = 1)
  expect_equal(fc$forecast_trials, c(60, 60))
  expect_equal(unname(apply(fc$draws, c(1, 2), sum)), matrix(60, 50, 2))
})

test_that("forecast() dispatches through the generics generic", {
  sim <- simulate_dynamic_poisson(n = 30, sigma = 0.2, log_rate0 = 2, seed = 11)
  fit <- fit_dynamic_model(sim$y, nsave = 50, nburn = 30, seed = 11)
  expect_s3_class(generics::forecast(fit, horizon = 2), "dynamic_forecast")
})
