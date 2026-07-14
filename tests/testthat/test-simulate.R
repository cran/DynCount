test_that("simulate_dynamic_poisson returns a well-formed list", {
  sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 1.5, seed = 1)
  expect_named(sim, c("y", "log_rate", "rate", "offset", "structural"))
  expect_length(sim$y, 60)
  expect_length(sim$log_rate, 60)
  expect_true(all(sim$y >= 0))
  expect_true(all(sim$y == round(sim$y)))      # integer counts
  expect_true(all(sim$rate > 0))
  expect_equal(sim$rate, exp(sim$log_rate))
})

test_that("simulate_dynamic_poisson is reproducible with a seed", {
  a <- simulate_dynamic_poisson(n = 40, sigma = 0.3, seed = 42)
  b <- simulate_dynamic_poisson(n = 40, sigma = 0.3, seed = 42)
  expect_identical(a$y, b$y)
})

test_that("zero inflation injects structural zeros", {
  sim <- simulate_dynamic_poisson(n = 500, sigma = 0.1, log_rate0 = 3,
                                  zero_inflation = 0.5, seed = 7)
  expect_true(sum(sim$structural) > 0)
  expect_true(all(sim$y[sim$structural] == 0))
  expect_gt(mean(sim$structural), 0.3)
  expect_lt(mean(sim$structural), 0.7)
})

test_that("simulate_dynamic_binomial respects the trial ceiling", {
  sim <- simulate_dynamic_binomial(n = 80, sigma = 0.15, trials = 30, seed = 3)
  expect_named(sim, c("y", "trials", "logit", "prob", "offset", "structural"))
  expect_true(all(sim$y <= sim$trials))
  expect_true(all(sim$y >= 0))
  expect_true(all(sim$prob >= 0 & sim$prob <= 1))
  expect_equal(sim$prob, plogis(sim$logit))
})

test_that("simulate_dynamic_binomial accepts a vector of trials", {
  tr <- sample(20:60, 50, replace = TRUE)
  sim <- simulate_dynamic_binomial(n = 50, sigma = 0.1, trials = tr, seed = 9)
  expect_identical(sim$trials, tr)
  expect_true(all(sim$y <= tr))
})

test_that("simulate_dynamic_binomial injects structural zeros", {
  sim <- simulate_dynamic_binomial(n = 400, sigma = 0.08, trials = 40, logit0 = 2,
                                   zero_inflation = 0.4, seed = 11)
  expect_true(all(sim$y[sim$structural] == 0))
  expect_gt(mean(sim$structural), 0.25)
  expect_lt(mean(sim$structural), 0.55)
})

test_that("invalid simulation arguments error", {
  expect_error(simulate_dynamic_poisson(n = 10, sigma = -1))
  expect_error(simulate_dynamic_poisson(n = 10, sigma = 0.1, zero_inflation = 1.5))
  expect_error(simulate_dynamic_binomial(n = 10, sigma = 0.1, trials = c(5, 6)))
})
