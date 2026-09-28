# Multinomial (ALR) family -----------------------------------------------------

test_that("alr_to_prob is a softmax with the baseline fixed at zero", {
  eta <- matrix(c(0.5, -1, 2, 0), 2, 2)
  P <- alr_to_prob(eta, baseline = 2)
  expect_equal(dim(P), c(2, 3))
  expect_equal(rowSums(P), c(1, 1))
  full <- cbind(eta[, 1], 0, eta[, 2])
  expect_equal(P, exp(full) / rowSums(exp(full)))
  # numerically stable for large predictors
  expect_true(all(is.finite(alr_to_prob(matrix(c(800, -800), 1), 1))))
})

test_that("simulate_dynamic_multinomial returns a coherent object", {
  sim <- simulate_dynamic_multinomial(n = 30, sigma = c(0.1, 0.2), trials = 50,
                                      alr0 = c(-1, 1), baseline = 2, seed = 1)
  expect_equal(dim(sim$y), c(30, 3))
  expect_equal(rowSums(sim$y), rep(50, 30))
  expect_equal(dim(sim$alr), c(30, 2))
  expect_equal(rowSums(sim$prob), rep(1, 30))
  expect_equal(sim$baseline, 2L)
  expect_equal(colnames(sim$y), c("cat1", "cat2", "cat3"))
  expect_error(simulate_dynamic_multinomial(30, sigma = c(0.1, 0.2, 0.3), trials = 5,
                                            alr0 = c(0, 0)), "length")
})

test_that("multinomial input validation", {
  sim <- simulate_dynamic_multinomial(n = 20, sigma = 0.1, trials = 40, seed = 2)
  expect_error(fit_dynamic_model(sim$y[, 1], family = "multinomial"), "matrix")
  expect_error(fit_dynamic_model(sim$y, family = "multinomial", trials = 40), "trials")
  expect_error(fit_dynamic_model(sim$y, family = "multinomial", zeros = "inflated"),
               "multinomial")
  expect_error(fit_dynamic_model(sim$y, family = "multinomial", zero_inflation = TRUE),
               "multinomial")
  expect_error(fit_dynamic_model(sim$y, family = "multinomial", baseline = "nope"),
               "column")
  expect_error(fit_dynamic_model(sim$y, family = "multinomial", baseline = 7), "1..K")
  expect_error(fit_dynamic_model(sim$y, family = "multinomial",
                                 offset = matrix(0, 20, 3)), "K - 1")
  ybad <- sim$y; ybad[2, 1] <- -1
  expect_error(fit_dynamic_model(ybad, family = "multinomial"), "non-negative")
})

test_that("the baseline defaults to the largest category and can be chosen", {
  sim <- simulate_dynamic_multinomial(n = 30, sigma = 0.1, trials = 60,
                                      alr0 = c(-1, 0.5), baseline = 2,
                                      categories = c("A", "B", "C"), seed = 3)
  largest <- which.max(colSums(sim$y))
  fit <- fit_dynamic_model(sim$y, family = "multinomial", nsave = 40, nburn = 20, seed = 3)
  expect_equal(fit$data$baseline, largest)
  expect_equal(fit$data$baseline_name, colnames(sim$y)[largest])
  expect_equal(fit$data$categories, c("A", "B", "C"))
  fit_a <- fit_dynamic_model(sim$y, family = "multinomial", baseline = "A",
                             nsave = 40, nburn = 20, seed = 3)
  expect_equal(fit_a$data$baseline, 1L)
  fit_3 <- fit_dynamic_model(sim$y, family = "multinomial", baseline = 3,
                             nsave = 40, nburn = 20, seed = 3)
  expect_equal(fit_3$data$baseline_name, "C")
  # unnamed input gets default labels; data frames are accepted
  fit_u <- fit_dynamic_model(unname(sim$y), family = "multinomial",
                             nsave = 20, nburn = 10, seed = 3)
  expect_equal(fit_u$data$categories, c("cat1", "cat2", "cat3"))
  fit_d <- fit_dynamic_model(as.data.frame(sim$y), family = "multinomial",
                             nsave = 20, nburn = 10, seed = 3)
  expect_equal(fit_d$data$categories, c("A", "B", "C"))
})

test_that("multinomial draws have the documented array layout", {
  sim <- simulate_dynamic_multinomial(n = 25, sigma = 0.1, trials = 80,
                                      alr0 = c(0, 0, 0), seed = 4)     # K = 4
  fit <- fit_dynamic_model(sim$y, family = "multinomial", nsave = 30, nburn = 10,
                           horizon = 3, seed = 4)
  d <- fit$draws
  lat <- fit$data$categories[-fit$data$baseline]
  expect_equal(dim(d$z), c(30, 25, 3))
  expect_equal(dimnames(d$z)[[3]], lat)
  expect_equal(dim(d$sig2), c(30, 25, 3))
  expect_equal(dim(d$z0), c(30, 3))
  expect_equal(dim(d$fitted), c(30, 25, 4))
  expect_equal(dimnames(d$fitted)[[3]], fit$data$categories)
  expect_equal(dim(d$yrep), c(30, 25, 4))
  expect_equal(dim(d$fitted_prob), c(30, 25, 4))
  expect_equal(dim(d$rho), c(30, 3)); expect_equal(dim(d$mu), c(30, 3))
  expect_equal(dim(d$innov_var), c(30, 3))
  expect_equal(dim(d$forecast_z), c(30, 3, 3))
  expect_equal(dim(d$forecast_y), c(30, 3, 4))
  expect_equal(dim(d$forecast_prob), c(30, 3, 4))
  expect_null(d$fitted_open); expect_null(d$yrep_open); expect_null(d$gate)
  # shares sum to one, expected counts to the row totals, replicates to trials
  expect_equal(unname(apply(d$fitted_prob, c(1, 2), sum)), matrix(1, 30, 25))
  expect_equal(unname(apply(d$fitted, c(1, 2), sum)),
               matrix(rep(sim$trials, each = 30), 30, 25))
  expect_equal(unname(apply(d$yrep, c(1, 2), sum)),
               matrix(rep(sim$trials, each = 30), 30, 25))
  expect_equal(unname(apply(d$forecast_y, c(1, 2), sum)), matrix(80, 30, 3))
  expect_equal(fit$spec$forecast_trials, rep(80, 3))
  expect_equal(fit$data$trials, sim$trials)
})

test_that("the univariate families keep their matrix/vector layout", {
  sim <- simulate_dynamic_poisson(n = 30, sigma = 0.2, seed = 5)
  fit <- fit_dynamic_model(sim$y, nsave = 20, nburn = 10, horizon = 2, seed = 5)
  expect_true(is.matrix(fit$draws$z)); expect_null(dim(fit$draws$rho))
  expect_true(is.matrix(fit$draws$forecast_y))
  expect_null(fit$draws$fitted_prob); expect_null(fit$draws$forecast_prob)
})

test_that("multinomial model recovers the category shares", {
  sim <- simulate_dynamic_multinomial(n = 100, sigma = c(0.15, 0.1), trials = 300,
                                      alr0 = c(-1, 0.5), baseline = 2, seed = 6)
  fit <- fit_dynamic_model(sim$y, family = "multinomial", baseline = 2,
                           nsave = 500, nburn = 300, seed = 6)
  phat <- apply(fit$draws$fitted_prob, c(2, 3), mean)
  for (k in 1:3) expect_gt(cor(phat[, k], sim$prob[, k]), 0.8)
  expect_lt(max(abs(phat - sim$prob)), 0.1)
  # latent ALR paths are recovered on the right scale
  n <- 100
  zhat <- apply(fit$draws$z, c(2, 3), mean)
  expect_gt(cor(zhat[, 1], sim$alr[, 1]), 0.8)
  expect_gt(cor(zhat[, 2], sim$alr[, 2]), 0.8)
})

test_that("K = 2 multinomial with the second column as baseline is the binomial model", {
  sim <- simulate_dynamic_binomial(n = 80, sigma = 0.12, trials = 60, logit0 = 0.5,
                                   seed = 7)
  Y <- cbind(success = sim$y, failure = sim$trials - sim$y)
  fit_m <- fit_dynamic_model(Y, family = "multinomial", baseline = "failure",
                             nsave = 500, nburn = 300, seed = 7)
  fit_b <- fit_dynamic_model(sim$y, family = "binomial", trials = sim$trials,
                             nsave = 500, nburn = 300, seed = 7)
  # same latent series (the ALR is the logit): posterior means agree closely
  zm <- colMeans(fit_m$draws$z[, , 1])
  zb <- colMeans(fit_b$draws$z)
  expect_gt(cor(zm, zb), 0.98)
  expect_lt(mean(abs(zm - zb)), 0.1)
  expect_lt(abs(mean(fit_m$draws$innov_var[, 1]) - mean(fit_b$draws$innov_var)) /
              mean(fit_b$draws$innov_var), 0.5)
  # fitted counts of the success category match the binomial fitted values
  fm <- colMeans(fit_m$draws$fitted[, , "success"])
  fb <- colMeans(fit_b$draws$fitted)
  expect_lt(max(abs(fm - fb)) / max(fb), 0.1)
})

test_that("parameters are not shared across categories", {
  # one very smooth and one very rough ALR series (relative to baseline 3);
  # the baseline must be fixed to the simulated one, because ALR series
  # relative to a different (rough) baseline are all rough
  sim <- simulate_dynamic_multinomial(n = 120, sigma = c(0.02, 0.4), trials = 400,
                                      alr0 = c(0, 0), baseline = 3, seed = 8)
  fit <- fit_dynamic_model(sim$y, family = "multinomial", baseline = 3,
                           nsave = 500, nburn = 300, seed = 8)
  sds <- sqrt(colMeans(fit$draws$innov_var))
  expect_lt(sds[1], 0.1)
  expect_lt(sds[1], sds[2] / 3)
  expect_lt(abs(sds[2] - 0.4), 0.15)
})

test_that("all dynamics / innovation structures run for the multinomial family", {
  sim <- simulate_dynamic_multinomial(n = 40, sigma = 0.15, trials = 100,
                                      alr0 = c(-0.5, 0.5), rho = 0.8, mu = 0, seed = 9)
  for (innov in c("gaussian", "t", "mixture")) {
    for (dyn in c("rw", "ar1")) {
      fit <- fit_dynamic_model(sim$y, family = "multinomial", innovations = innov,
                               latent_dynamics = dyn, include_mu = TRUE,
                               nsave = 60, nburn = 40, seed = 9)
      expect_s3_class(fit, "dynamic_fit")
      expect_true(all(is.finite(fit$draws$z)))
      s <- summary(fit)
      expect_true(all(grepl("\\[", rownames(s$params))))
      if (dyn == "ar1") {
        expect_true(all(fit$draws$rho > -1 & fit$draws$rho < 1))
        expect_true(any(grepl("^ar1_rho\\[", rownames(s$params))))
      }
      if (innov == "t") expect_true(any(grepl("^t_df\\[", rownames(s$params))))
    }
  }
})

test_that("multinomial offsets shift the ALR predictor", {
  sim <- simulate_dynamic_multinomial(n = 60, sigma = 0.1, trials = 200,
                                      alr0 = c(0, 0), baseline = 3, seed = 10)
  off <- cbind(rep(1, 60), rep(-1, 60))
  fit <- fit_dynamic_model(sim$y, family = "multinomial", baseline = 3, offset = off,
                           nsave = 150, nburn = 100, horizon = 3,
                           forecast_offset = c(1, -1), seed = 10)
  expect_equal(fit$data$offset, off)
  expect_equal(fit$spec$forecast_offset, matrix(c(1, -1), 3, 2, byrow = TRUE))
  # the latent state absorbs the offset: z ~ alr - offset
  zhat <- apply(fit$draws$z, c(2, 3), mean)
  expect_lt(abs(mean(zhat[, 1]) - mean(sim$alr[, 1]) + 1), 0.4)
  expect_lt(abs(mean(zhat[, 2]) - mean(sim$alr[, 2]) - 1), 0.4)
  # a scalar offset is recycled, and so is a per-category vector
  fit0 <- fit_dynamic_model(sim$y, family = "multinomial", offset = 0.5,
                            nsave = 20, nburn = 10, seed = 10)
  expect_equal(fit0$data$offset, matrix(0.5, 60, 2))
  fit1 <- fit_dynamic_model(sim$y, family = "multinomial", offset = c(1, -1),
                            nsave = 20, nburn = 10, seed = 10)
  expect_equal(fit1$data$offset, matrix(c(1, -1), 60, 2, byrow = TRUE))
})

test_that("zero-total rows are allowed and treated as missing", {
  sim <- simulate_dynamic_multinomial(n = 40, sigma = 0.1, trials = 50, seed = 11)
  y <- sim$y; y[c(5, 20), ] <- 0
  fit <- fit_dynamic_model(y, family = "multinomial", nsave = 60, nburn = 30, seed = 11)
  expect_equal(fit$data$trials[c(5, 20)], c(0, 0))
  expect_true(all(is.finite(fit$draws$z)))
  expect_true(all(fit$draws$yrep[, 5, ] == 0))
  expect_equal(summary(fit)$n_zero, 2L)
  expect_error(fit_dynamic_model(0 * y, family = "multinomial"), "zero")
})

test_that("predict, forecast and plots work for multinomial fits", {
  sim <- simulate_dynamic_multinomial(n = 30, sigma = 0.1, trials = 70,
                                      alr0 = c(-0.5, 0.5),
                                      categories = c("A", "B", "C"), seed = 12)
  fit <- fit_dynamic_model(sim$y, family = "multinomial", nsave = 80, nburn = 40,
                           horizon = 4, forecast_trials = 90, seed = 12)
  p <- predict(fit)
  expect_equal(nrow(p$summary), 30 * 3)
  expect_named(p$summary[, 1:3], c("time", "category", "observed"))
  expect_equal(p$summary$observed[p$summary$category == "B"], unname(sim$y[, "B"]))
  pp <- predict(fit, type = "prob")
  expect_true(all(pp$summary$mean >= 0 & pp$summary$mean <= 1))
  expect_equal(predict(fit, type = "response", conditional = TRUE)$draws, fit$draws$yrep)
  expect_error(predict(fit_dynamic_model(sim$y[, 1], nsave = 10, nburn = 5),
                       type = "prob"), "multinomial")

  fc <- forecast(fit)
  expect_s3_class(fc, "dynamic_forecast")
  expect_equal(nrow(fc$summary), 4 * 3)
  expect_equal(nrow(fc$latent), 4 * 2)
  expect_equal(nrow(fc$prob), 4 * 3)
  expect_equal(nrow(fc$final), 3)
  expect_equal(dim(fc$final_draws), c(80, 3))
  expect_equal(unname(apply(fc$draws, c(1, 2), sum)), matrix(90, 80, 4))
  expect_output(print(fc), "multinomial")
  expect_output(print(fit), "baseline")
  expect_output(print(summary(fit)), "categories")

  expect_error(structural_zero_prob(fit), "multinomial")
  expect_error(plot_zero_inflation(fit), "multinomial")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot_latent(fit))
  expect_silent(plot_fitted(fit, category = c("A", "C")))
  expect_silent(plot_forecast(fit, category = 2))
  expect_error(plot_fitted(fit, category = "Z"), "Unknown")
  expect_silent(plot(fit, "latent", category = "A"))
})
