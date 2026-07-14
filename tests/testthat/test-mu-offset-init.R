test_that("RW without drift keeps mu fixed at zero", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 1)
  fit <- fit_dynamic_model(sim$y, nsave = 225, nburn = 112, seed = 1)
  expect_false(isTRUE(fit$spec$include_mu))
  expect_equal(unique(fit$draws$mu), 0)        # never sampled
})

test_that("RW with drift fixes rho = 1 and samples mu", {
  sim <- simulate_dynamic_poisson(n = 150, sigma = 0.15, log_rate0 = 2,
                                  mu = 0.03, seed = 1)
  fit <- fit_dynamic_model(sim$y, include_mu = TRUE, nsave = 450, nburn = 300,
                           seed = 1)
  expect_true(fit$spec$include_mu)
  expect_equal(unique(fit$draws$rho), 1)        # rho fixed
  expect_gt(stats::sd(fit$draws$mu), 0)         # mu genuinely sampled
  expect_gt(mean(fit$draws$mu), 0)              # positive drift detected
  s <- summary(fit)
  expect_true("drift_mu" %in% rownames(s$params))
})

test_that("AR(1) silently defaults to including an intercept", {
  sim <- simulate_dynamic_poisson(n = 120, sigma = 0.2, log_rate0 = 3,
                                  rho = 0.7, mu = 0.9, seed = 2)   # stationary mean 3
  # include_mu is left at its default (FALSE) but AR(1) enables it automatically
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1",
                           nsave = 375, nburn = 225, seed = 2)
  expect_true(fit$spec$include_mu)              # auto-enabled
  expect_gt(stats::sd(fit$draws$mu), 0)         # intercept genuinely sampled
  expect_gt(stats::sd(fit$draws$rho), 0)        # rho sampled
})

test_that("AR(1) with intercept jointly updates (mu, rho)", {
  sim <- simulate_dynamic_poisson(n = 200, sigma = 0.2, log_rate0 = 4,
                                  rho = 0.8, mu = 0.5, seed = 2)
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1", include_mu = TRUE,
                           nsave = 600, nburn = 450, seed = 2)
  expect_gt(stats::sd(fit$draws$mu), 0)
  expect_gt(stats::sd(fit$draws$rho), 0)
  # the implied stationary mean mu/(1-rho) is well identified (~2.5)
  smean <- mean(fit$draws$mu / (1 - fit$draws$rho))
  expect_gt(smean, 1.5)
  expect_lt(smean, 3.5)
  s <- summary(fit)
  expect_true(all(c("ar1_rho", "intercept_mu") %in% rownames(s$params)))
  # the (rho, mu) update is a conjugate Gibbs draw, not Metropolis-Hastings:
  # consecutive stored draws (thin = 1) must never repeat (no rejections)
  expect_true(all(diff(fit$draws$rho) != 0))
  expect_true(all(diff(fit$draws$mu)  != 0))
})

test_that("rho stays within the stationary region with an intercept", {
  sim <- simulate_dynamic_poisson(n = 80, sigma = 0.2, log_rate0 = 3, seed = 4)
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1", include_mu = TRUE,
                           prior = dynamic_prior(ar_rho_mean = 1.5, ar_rho_sd = 0.05),
                           nsave = 375, nburn = 225, seed = 4)
  expect_true(all(fit$draws$rho > -1 & fit$draws$rho < 1))  # constrained
})

test_that("the GMRF stays sparse: D_rho is bidiagonal so D'WD is tridiagonal", {
  D <- make_gmrf_diff(6, rho = 0.7)
  W <- diag(runif(5) + 0.1)
  P <- crossprod(D, W %*% D)
  offband <- P[abs(row(P) - col(P)) > 1]
  expect_true(all(abs(offband) < 1e-12))
})

test_that("a drift is recovered and its forecast runs", {
  sim <- simulate_dynamic_poisson(n = 120, sigma = 0.1, log_rate0 = 1,
                                  mu = 0.05, seed = 5)
  fit <- fit_dynamic_model(sim$y, include_mu = TRUE, nsave = 375, nburn = 225,
                           horizon = 20, seed = 5)
  # the positive drift is recovered in the posterior of mu
  expect_gt(mean(fit$draws$mu), 0)
  # and a forecast is produced over the full horizon
  fc <- forecast(fit)
  expect_equal(nrow(fc$latent), 20)
  expect_true(all(is.finite(fc$latent[["mean"]])))
})

test_that("plot_forecast works for a drift model", {
  sim <- simulate_dynamic_poisson(n = 80, sigma = 0.1, log_rate0 = 1, mu = 0.04,
                                  seed = 6)
  fit <- fit_dynamic_model(sim$y, include_mu = TRUE, nsave = 225, nburn = 112,
                           horizon = 15, seed = 6)
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  expect_silent(plot_forecast(fit))
  dev.off(); expect_true(file.exists(tmp)); unlink(tmp)
})

test_that("a Poisson offset is absorbed into the mean", {
  off <- rep(log(100), 120)
  sim <- simulate_dynamic_poisson(120, 0.1, log_rate0 = 0, offset = off, seed = 3)
  expect_gt(mean(sim$y), 50)                    # offset lifts the counts
  fit <- fit_dynamic_model(sim$y, offset = off, nsave = 225, nburn = 150, seed = 3)
  fm <- colMeans(fit$draws$fitted)
  expect_gt(mean(fm), 40)
  expect_lt(mean(fm), 250)
})

test_that("offset validation and forecasting with forecast_offset", {
  sim <- simulate_dynamic_poisson(60, 0.15, log_rate0 = 1, seed = 7)
  expect_error(fit_dynamic_model(sim$y, offset = c(1, 2)), "length")
  fit <- fit_dynamic_model(sim$y, offset = 0.5, nsave = 225, nburn = 112,
                           horizon = 6, forecast_offset = 0.5, seed = 7)
  fc <- forecast(fit)
  expect_equal(nrow(fc$summary), 6)
  expect_error(
    fit_dynamic_model(sim$y, offset = 0.5, nsave = 150, nburn = 75,
                      horizon = 6, forecast_offset = c(1, 2), seed = 7),
    "length"
  )
})

test_that("the initial-state prior anchors the first latent state (RW)", {
  sim <- simulate_dynamic_poisson(60, 0.2, log_rate0 = 2, seed = 8)
  fit <- fit_dynamic_model(sim$y,
                           prior = dynamic_prior(init_mean = 5, init_var = 1e-3),
                           nsave = 300, nburn = 150, seed = 8)
  expect_lt(abs(mean(fit$draws$z[, 1]) - 5), 0.5)
})

test_that("AR(1) with the fixed diffuse initial prior runs and gives finite draws", {
  sim <- simulate_dynamic_poisson(120, 0.2, log_rate0 = 3, rho = 0.7,
                                  mu = 0.9, seed = 9)   # stationary mean 3
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1", include_mu = TRUE,
                           nsave = 300, nburn = 225, seed = 9)
  expect_true(all(is.finite(fit$draws$z)))
  expect_true(all(is.finite(fit$draws$rho)))
})

test_that("include_mu forecasts run at the final h-step", {
  sim <- simulate_dynamic_poisson(71, 0.15, log_rate0 = 2, mu = 0.02, seed = 10)
  fit <- fit_dynamic_model(sim$y[1:60], include_mu = TRUE,
                           horizon = 5, nsave = 300, nburn = 150, seed = 10)
  expect_equal(fit$spec$horizon, 5)
  expect_true(all(is.finite(forecast(fit)$summary$mean)))
})
