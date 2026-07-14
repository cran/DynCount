test_that("RW is the default and fixes rho = 1 without sampling it", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 1)
  fit <- fit_dynamic_model(sim$y, nsave = 225, nburn = 112, seed = 1)
  expect_equal(fit$spec$latent_dynamics, "rw")          # default
  expect_true(!is.null(fit$draws$rho))
  expect_equal(unique(fit$draws$rho), 1)                # never sampled
})

test_that("AR(1) samples rho (it varies across draws)", {
  sim <- simulate_dynamic_poisson(n = 120, sigma = 0.2, log_rate0 = 3, seed = 2)
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1", include_mu = TRUE,
                           nsave = 450, nburn = 300, seed = 2)
  expect_equal(fit$spec$latent_dynamics, "ar1")
  expect_gt(stats::sd(fit$draws$rho), 0)                # genuinely sampled
  expect_true(all(is.finite(fit$draws$rho)))
  s <- summary(fit)
  expect_true("ar1_rho" %in% rownames(s$params))        # reported in summary
})

test_that("rho is constrained to the stationary region (-1, 1)", {
  # a tight prior centred outside (-1, 1) would, without a constraint, pull the
  # posterior past 1; the stationarity truncation keeps every draw inside (-1, 1)
  sim <- simulate_dynamic_poisson(n = 80, sigma = 0.2, log_rate0 = 3, seed = 4)
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1", include_mu = TRUE,
                           prior = dynamic_prior(ar_rho_mean = 1.5, ar_rho_sd = 0.05),
                           nsave = 375, nburn = 225, seed = 4)
  expect_true(all(fit$draws$rho > -1 & fit$draws$rho < 1))  # never escapes (-1, 1)
  expect_gt(mean(fit$draws$rho), 0.5)                       # but pulled toward 1
})

test_that("AR(1) recovers rho in a high-signal regime", {
  sim <- simulate_dynamic_poisson(n = 200, sigma = 0.2, log_rate0 = 4,
                                  rho = 0.95, mu = 0.2, seed = 3)   # stationary mean 4
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1", include_mu = TRUE,
                           nsave = 600, nburn = 450, seed = 3)
  expect_gt(mean(fit$draws$rho), 0.8)                   # near the true 0.95
  expect_lt(mean(fit$draws$rho), 1)                     # strictly stationary
})

test_that("the GMRF difference operator generalises the random walk", {
  D1 <- make_gmrf_diff(5, rho = 1)
  expect_equal(D1, diff(diag(5)))                       # rho = 1 -> first diff
  Dr <- make_gmrf_diff(5, rho = 0.7)
  expect_equal(diag(Dr[, -ncol(Dr)]), rep(-0.7, 4))     # -rho on the subdiagonal
})

test_that("arbitrary forecast horizons work and expose the final step", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 5)
  for (h in c(1, 7, 30)) {
    fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 150, horizon = h, seed = 5)
    fc <- forecast(fit)
    expect_equal(fc$horizon, h)
    expect_equal(nrow(fc$summary), h)
    expect_equal(nrow(fc$final), 1L)
    expect_equal(fc$final$horizon, h)
    expect_length(fc$final_draws, nrow(fit$draws$z))
  }
})

test_that("forecast errors when the model was fitted without a horizon", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 5)
  fit <- fit_dynamic_model(sim$y, nsave = 225, nburn = 112, seed = 5)  # horizon = 0
  expect_error(forecast(fit), "horizon")
})

test_that("AR(1) forecasts use rho in the propagation", {
  sim <- simulate_dynamic_poisson(n = 120, sigma = 0.2, log_rate0 = 3, seed = 6)
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1", include_mu = TRUE,
                           nsave = 300, nburn = 225, horizon = 5, seed = 6)
  fc <- forecast(fit)
  expect_equal(nrow(fc$summary), 5)
  expect_true(all(is.finite(fc$latent[["mean"]])))
})

test_that("horizon is honoured for multi-step forecasts", {
  sim <- simulate_dynamic_poisson(n = 71, sigma = 0.2, log_rate0 = 2, seed = 7)
  fit1 <- fit_dynamic_model(sim$y[1:60], horizon = 1,
                            nsave = 300, nburn = 150, seed = 7)
  fit5 <- fit_dynamic_model(sim$y[1:60], horizon = 5,
                            nsave = 300, nburn = 150, seed = 7)
  expect_equal(fit1$spec$horizon, 1)
  expect_equal(fit5$spec$horizon, 5)
  expect_equal(nrow(forecast(fit5)$summary), 5)
})

test_that("plot_forecast works for multi-step horizons", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 9)
  fit <- fit_dynamic_model(sim$y, nsave = 225, nburn = 112, horizon = 25, seed = 9)
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  expect_silent(plot_forecast(fit))
  dev.off()
  expect_true(file.exists(tmp))
  unlink(tmp)
})
