test_that("binomial fit requires trials", {
  sim <- simulate_dynamic_binomial(n = 50, sigma = 0.15, trials = 40, seed = 1)
  expect_error(fit_dynamic_model(sim$y, family = "binomial"), "trials")
})

test_that("binomial model recovers the latent probability", {
  sim <- simulate_dynamic_binomial(n = 120, sigma = 0.12, trials = 50, seed = 4)
  fit <- fit_dynamic_model(sim$y, family = "binomial", trials = sim$trials,
                           nsave = 750, nburn = 375, seed = 4)
  expect_s3_class(fit, "dynamic_fit")
  expect_equal(fit$spec$family, "binomial")
  fitted_prob <- colMeans(fit$draws$fitted) / sim$trials
  expect_gt(cor(fitted_prob, sim$prob), 0.8)
})

test_that("trials are recycled from a single value", {
  sim <- simulate_dynamic_binomial(n = 40, sigma = 0.1, trials = 25, seed = 8)
  fit <- fit_dynamic_model(sim$y, family = "binomial", trials = 25,
                           nsave = 225, nburn = 112, seed = 8)
  expect_equal(fit$data$trials, rep(25, 40))
})

test_that("y greater than trials is rejected", {
  y <- c(5, 12, 3, 8)
  expect_error(
    fit_dynamic_model(y, family = "binomial", trials = 10),
    "trials"
  )
})

test_that("zero-inflated binomial runs and exposes a gate", {
  sim <- simulate_dynamic_binomial(n = 80, sigma = 0.12, trials = 30,
                                   zero_inflation = 0.2, seed = 2)
  fit <- fit_dynamic_model(sim$y, family = "binomial", trials = 30,
                           zeros = "inflated", nsave = 375, nburn = 188, seed = 2)
  expect_equal(fit$spec$zeros, "inflated")
  expect_true(!is.null(fit$draws$gate))
  expect_true(!is.null(fit$draws$pi_open))
})

test_that("structural_zero_prob works for the binomial family", {
  set.seed(123)
  z <- 1.5 + cumsum(c(0, rnorm(59, 0, 0.05)))   # logit ~ 1.5 => p ~ 0.82
  trials <- rep(30, 60)
  y <- rbinom(60, trials, plogis(z))
  y[c(10, 30, 50)] <- 0                          # inject structural zeros
  fit <- fit_dynamic_model(y, family = "binomial", trials = trials,
                           zeros = "inflated", nsave = 525, nburn = 262, seed = 123)
  sz <- structural_zero_prob(fit, zeros_only = TRUE)
  expect_true(all(sz$observed == 0))
  injected <- sz$p_structural[sz$time %in% c(10, 30, 50)]
  expect_true(all(injected > 0.5))
})

test_that("zero_inflation = TRUE alias works for binomial", {
  sim <- simulate_dynamic_binomial(n = 50, sigma = 0.1, trials = 25,
                                   zero_inflation = 0.15, seed = 5)
  fit <- fit_dynamic_model(sim$y, family = "binomial", trials = 25,
                           zero_inflation = TRUE, nsave = 225, nburn = 112, seed = 5)
  expect_equal(fit$spec$zeros, "inflated")
})
