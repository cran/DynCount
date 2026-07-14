test_that("predict.dynamic_fit returns fitted summaries aligned to data", {
  sim <- simulate_dynamic_poisson(n = 50, sigma = 0.2, log_rate0 = 2, seed = 1)
  fit <- fit_dynamic_model(sim$y, family = "poisson",
                           nsave = 300, nburn = 150, seed = 1)
  pr <- predict(fit, type = "mean")
  expect_true(is.list(pr))
  expect_equal(nrow(pr$summary), 50)
  expect_equal(pr$summary$observed, sim$y)
  expect_true(all(c("mean", "sd") %in% names(pr$summary)))
})

test_that("predict response draws are non-negative integers", {
  sim <- simulate_dynamic_poisson(n = 40, sigma = 0.2, seed = 2)
  fit <- fit_dynamic_model(sim$y, family = "poisson",
                           nsave = 225, nburn = 112, seed = 2)
  pr <- predict(fit, type = "response")
  expect_true(all(pr$draws >= 0))
  expect_true(all(pr$draws == round(pr$draws)))
})

test_that("forecast produces a dynamic_forecast of the requested horizon", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 3)
  fit <- fit_dynamic_model(sim$y, family = "poisson",
                           nsave = 300, nburn = 150, horizon = 8, seed = 3)
  fc <- forecast(fit)
  expect_s3_class(fc, "dynamic_forecast")
  expect_equal(fc$horizon, 8)
  expect_equal(nrow(fc$summary), 8)
  expect_equal(nrow(fc$latent), 8)
  expect_equal(fc$family, "poisson")
  expect_true(all(fc$summary[["q2.5"]] <= fc$summary[["q97.5"]]))
})

test_that("binomial predict response draws are integers within [0, trials]", {
  sim <- simulate_dynamic_binomial(n = 50, sigma = 0.12, trials = 40, seed = 5)
  fit <- fit_dynamic_model(sim$y, family = "binomial", trials = sim$trials,
                           nsave = 225, nburn = 112, seed = 5)
  pr <- predict(fit, type = "response")
  expect_true(all(pr$draws == round(pr$draws)))
  expect_true(all(pr$draws >= 0))
  # each column is bounded by that observation's trial count
  expect_true(all(sweep(pr$draws, 2L, sim$trials, `<=`)))
})

test_that("binomial predict mean lies within [0, trials]", {
  sim <- simulate_dynamic_binomial(n = 40, sigma = 0.12, trials = 30, seed = 6)
  fit <- fit_dynamic_model(sim$y, family = "binomial", trials = 30,
                           nsave = 225, nburn = 112, seed = 6)
  pr <- predict(fit, type = "mean")
  expect_true(all(pr$summary$mean >= 0 & pr$summary$mean <= 30))
})

test_that("binomial forecasts respect new trial sizes", {
  sim <- simulate_dynamic_binomial(n = 50, sigma = 0.12, trials = 40, seed = 4)
  fit <- fit_dynamic_model(sim$y, family = "binomial", trials = 40,
                           nsave = 300, nburn = 150, horizon = 5,
                           forecast_trials = 100, seed = 4)
  fc <- forecast(fit)
  expect_equal(fc$family, "binomial")
  expect_true(all(fc$draws <= 100))
})
