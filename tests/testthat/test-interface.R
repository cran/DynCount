test_that("zero_inflation must be a single logical and must not conflict", {
  sim <- simulate_dynamic_poisson(n = 30, sigma = 0.2, log_rate0 = 2, seed = 1)
  expect_error(fit_dynamic_model(sim$y, zero_inflation = 0.2), "TRUE or FALSE")
  expect_error(fit_dynamic_model(sim$y, zero_inflation = NA), "TRUE or FALSE")
  expect_error(fit_dynamic_model(sim$y, zero_inflation = TRUE, zeros = "none"),
               "conflicts")
})

test_that("ignored arguments trigger warnings", {
  sim <- simulate_dynamic_poisson(n = 30, sigma = 0.2, log_rate0 = 2, seed = 2)
  expect_warning(fit_dynamic_model(sim$y, trials = 10, nsave = 10, nburn = 5),
                 "ignored for the Poisson")
  expect_warning(fit_dynamic_model(sim$y, baseline = 2, nsave = 10, nburn = 5),
                 "multinomial")
  expect_warning(fit_dynamic_model(sim$y, forecast_offset = 1, nsave = 10, nburn = 5),
                 "horizon = 0")
})

test_that("a seed does not alter the global random number stream", {
  sim <- simulate_dynamic_poisson(n = 30, sigma = 0.2, log_rate0 = 2, seed = 3)
  set.seed(99); a <- stats::runif(3)
  set.seed(99)
  fit <- fit_dynamic_model(sim$y, nsave = 20, nburn = 5, seed = 1)
  invisible(forecast(fit, horizon = 2, seed = 2))
  invisible(simulate_dynamic_binomial(10, 0.1, trials = 5, seed = 4))
  b <- stats::runif(3)
  expect_identical(a, b)
  # and the same seed gives the same fit
  fit2 <- fit_dynamic_model(sim$y, nsave = 20, nburn = 5, seed = 1)
  expect_identical(fit$draws$z, fit2$draws$z)
})

test_that("plot arguments in ... override the defaults", {
  sim <- simulate_dynamic_poisson(n = 30, sigma = 0.2, log_rate0 = 2, seed = 5)
  fit <- fit_dynamic_model(sim$y, nsave = 50, nburn = 20, seed = 5)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot_fitted(fit, main = "my title", ylab = "y"))
  expect_no_error(plot_latent(fit, main = "latent", ylim = c(-5, 10)))
  expect_no_error(plot_forecast(fit, horizon = 3, main = "forecast"))
})

test_that("summary reports SV parameters but no per-component mixture rows", {
  sim <- simulate_dynamic_poisson(n = 40, sigma = 0.2, log_rate0 = 2, seed = 6)
  fit <- fit_dynamic_model(sim$y, innovations = "mixture", nsave = 60, nburn = 30, seed = 6)
  expect_equal(rownames(summary(fit)$params), "innov_sd")
  # the component draws are stored, in order of increasing variance per draw
  expect_true(all(fit$draws$mix_var[, 1] <= fit$draws$mix_var[, 2]))
  expect_equal(rowSums(fit$draws$mix_weight), rep(1, 60))
  skip_if_not_installed("stochvol")
  fsv <- fit_dynamic_model(sim$y, innovations = "sv", nsave = 40, nburn = 20, seed = 6)
  expect_true(all(c("sv_mu", "sv_phi", "sv_sigma") %in% rownames(summary(fsv)$params)))
})

test_that("the vectorised multinomial log-likelihood is exact", {
  set.seed(1)
  Y <- matrix(c(3, 5, 1, 0, 7, 2), 2, 3)        # 2 observations, 3 categories
  Ntot <- rowSums(Y)
  base <- 3
  Ynb <- Y[, -base, drop = FALSE]
  Z <- matrix(rnorm(6), 3, 2)                   # padded latent matrix (N = n + 1)
  off <- matrix(0.1, 2, 2)
  zz <- c(0.4, -0.3)
  ll <- DynCount:::obs_loglik("multinomial", 1:2, 1L, zz, Z, Ynb, Ntot, off)
  for (j in 1:2) {
    eta <- off[j, ] + Z[j + 1, ]; eta[1] <- off[j, 1] + zz[j]
    p <- DynCount:::alr_to_prob(eta, base)
    expect_equal(ll[j], stats::dmultinom(Y[j, ], prob = p[1, ], log = TRUE))
  }
})

test_that("the multinomial simulator accepts a baseline name and a per-category offset", {
  sim <- simulate_dynamic_multinomial(n = 20, sigma = 0.1, trials = 50,
                                      alr0 = c(0, 0), offset = c(1, -1),
                                      categories = c("a", "b", "c"),
                                      baseline = "a", seed = 7)
  expect_equal(sim$baseline, 1L)
  expect_equal(sim$offset, matrix(c(1, -1), 20, 2, byrow = TRUE))
  expect_equal(colnames(sim$alr), c("b", "c"))
})

test_that("the example data carry a date column", {
  expect_s3_class(uk_weekly$date, "Date")
  expect_true(all(diff(med_weekly$date) == 7))
})
