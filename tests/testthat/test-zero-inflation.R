test_that("zero_inflation = TRUE turns on the inflated branch", {
  sim <- simulate_dynamic_poisson(n = 80, sigma = 0.2, log_rate0 = 2,
                                  zero_inflation = 0.2, seed = 3)
  fit <- fit_dynamic_model(sim$y, family = "poisson", zero_inflation = TRUE,
                           nsave = 375, nburn = 188, seed = 3)
  expect_equal(fit$spec$zeros, "inflated")
  expect_true(!is.null(fit$draws$gate))
  expect_true(!is.null(fit$draws$pi_open))
})

test_that("structural_zero_prob returns probabilities for observed zeros", {
  sim <- simulate_dynamic_poisson(n = 100, sigma = 0.2, log_rate0 = 3,
                                  zero_inflation = 0.3, seed = 12)
  fit <- fit_dynamic_model(sim$y, family = "poisson", zeros = "inflated",
                           nsave = 450, nburn = 225, seed = 12)
  sp <- structural_zero_prob(fit, zeros_only = TRUE)
  expect_true(all(sp$p_structural >= 0 & sp$p_structural <= 1))
  expect_true(all(sp$observed == 0))
})

test_that("structural zeros are flagged with high probability", {
  set.seed(99)
  z <- 4 + cumsum(c(0, rnorm(59, 0, 0.05)))
  y <- rpois(60, exp(z))
  y[c(10, 30, 50)] <- 0                      # inject obvious structural zeros
  fit <- fit_dynamic_model(y, family = "poisson", zeros = "inflated",
                           nsave = 600, nburn = 300, seed = 99)
  sp <- structural_zero_prob(fit, zeros_only = TRUE)
  injected <- sp$p_structural[sp$time %in% c(10, 30, 50)]
  expect_true(all(injected > 0.5))
})

test_that("structural_zero_prob errors without an inflated fit", {
  sim <- simulate_dynamic_poisson(n = 40, sigma = 0.2, seed = 1)
  fit <- fit_dynamic_model(sim$y, family = "poisson", zeros = "none",
                           nsave = 150, nburn = 75, seed = 1)
  expect_error(structural_zero_prob(fit), "inflat")
})

test_that("zeros = 'missing' runs and treats observed zeros as missing", {
  sim <- simulate_dynamic_poisson(n = 80, sigma = 0.2, log_rate0 = 1,
                                  zero_inflation = 0.3, seed = 3)
  fit <- fit_dynamic_model(sim$y, family = "poisson", zeros = "missing",
                           nsave = 150, nburn = 75, seed = 3)
  expect_equal(fit$spec$zeros, "missing")
  expect_true(all(is.finite(fit$draws$z)))
  # structural-zero diagnostics are only defined for the inflated model
  expect_error(structural_zero_prob(fit), "inflat")
})

# Conditional vs unconditional replicates and fitted values --------------------

test_that("yrep includes the gate while yrep_open is conditional on open", {
  sim <- simulate_dynamic_poisson(n = 120, sigma = 0.15, log_rate0 = 3,
                                  zero_inflation = 0.3, seed = 21)
  fit <- fit_dynamic_model(sim$y, zeros = "inflated",
                           nsave = 450, nburn = 225, seed = 21)
  p0_obs  <- mean(fit$data$y == 0)
  p0_rep  <- mean(fit$draws$yrep == 0)       # unconditional: gate applied
  p0_open <- mean(fit$draws$yrep_open == 0)  # conditional on gate open
  # the unconditional replicate reproduces the zero share of the data ...
  expect_lt(abs(p0_rep - p0_obs), 0.1)
  # ... while the conditional replicate has systematically too few zeros
  expect_lt(p0_open, p0_rep - 0.1)
})

test_that("fitted = pi_open * fitted_open exactly, per draw, under ZI", {
  sim <- simulate_dynamic_poisson(n = 80, sigma = 0.2, log_rate0 = 2,
                                  zero_inflation = 0.25, seed = 22)
  fit <- fit_dynamic_model(sim$y, zeros = "inflated",
                           nsave = 225, nburn = 150, seed = 22)
  expect_equal(fit$draws$fitted,
               fit$draws$fitted_open * fit$draws$pi_open)
})

test_that("without zero inflation the conditional/unconditional pairs coincide", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 23)
  fit <- fit_dynamic_model(sim$y, zeros = "none",
                           nsave = 150, nburn = 75, seed = 23)
  expect_identical(fit$draws$yrep,   fit$draws$yrep_open)
  expect_identical(fit$draws$fitted, fit$draws$fitted_open)
})

test_that("predict() exposes both flavours via the conditional argument", {
  sim <- simulate_dynamic_poisson(n = 100, sigma = 0.15, log_rate0 = 3,
                                  zero_inflation = 0.3, seed = 24)
  fit <- fit_dynamic_model(sim$y, zeros = "inflated",
                           nsave = 300, nburn = 150, seed = 24)
  pr_u <- predict(fit, type = "response")                      # unconditional
  pr_c <- predict(fit, type = "response", conditional = TRUE)  # gate-open
  expect_identical(pr_u$draws, fit$draws$yrep)
  expect_identical(pr_c$draws, fit$draws$yrep_open)
  expect_gt(mean(pr_u$draws == 0), mean(pr_c$draws == 0))
  pm_u <- predict(fit, type = "mean")
  pm_c <- predict(fit, type = "mean", conditional = TRUE)
  expect_true(all(pm_u$summary$mean <= pm_c$summary$mean + 1e-12))
})

test_that("forecast_y applies the gate: forecast zeros match the ZI mixture", {
  sim <- simulate_dynamic_poisson(n = 120, sigma = 0.1, log_rate0 = 3,
                                  zero_inflation = 0.4, seed = 25)
  fit <- fit_dynamic_model(sim$y, zeros = "inflated", horizon = 6,
                           nsave = 450, nburn = 225, seed = 25)
  # Under the full ZI model, P(Y = 0) per draw and horizon is
  #   (1 - pi_open) + pi_open * P(Poisson zero | forecast state),
  # so the empirical zero share of forecast_y must match its average.
  p0_pois  <- dpois(0, exp(fit$draws$forecast_z))       # draws x H
  p0_model <- mean((1 - fit$draws$pi_open) + fit$draws$pi_open * p0_pois)
  p0_fc    <- mean(fit$draws$forecast_y == 0)
  expect_lt(abs(p0_fc - p0_model), 0.03)
  # and the gate genuinely adds zeros beyond the pure observation model
  expect_gt(p0_fc, mean(p0_pois) + 0.1)
})
