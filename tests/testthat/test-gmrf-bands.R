# gmrf_bands(): the tridiagonal precision bands and the linear term b ---------
#
# These tests pin down the behaviour flagged in review: expanding the joint
# density of the GMRF plus the proper N(m0, v0) prior on the first state, the
# linear term must satisfy b = mu * (D_rho' w) + (m0 / v0) e_1. In particular
# b[1] contains m0 / v0 REGARDLESS of mu, so "b = 0 when mu = 0" only holds if
# the initial-state prior mean is itself zero.

dense_reference <- function(w, rho, mu, v0, m0, N) {
  D <- make_gmrf_diff(N, rho)
  P <- crossprod(D, diag(w, length(w)) %*% D)
  P[1, 1] <- P[1, 1] + 1 / v0
  b <- as.vector(mu * crossprod(D, w)) # mu * D_rho' w
  b[1] <- b[1] + m0 / v0               # init-prior contribution
  list(P = P, b = b)
}

test_that("gmrf_bands matches the dense P and b (rho != 1, mu != 0, m0 != 0)", {
  set.seed(1)
  N <- 7; w <- runif(N - 1, 0.5, 2)
  rho <- 0.6; mu <- 0.8; v0 <- 4; m0 <- -2.5
  bd  <- gmrf_bands(w, rho, mu, v0, m0, N)
  ref <- dense_reference(w, rho, mu, v0, m0, N)
  expect_equal(bd$diag, diag(ref$P))
  expect_equal(bd$off,  ref$P[cbind(1:(N - 1), 2:N)])
  expect_equal(bd$b,    ref$b)
})

test_that("with mu = 0 but init_mean != 0, b is NOT zero: b[1] = m0 / v0", {
  set.seed(2)
  N <- 6; w <- runif(N - 1, 0.5, 2)
  v0 <- 10; m0 <- 5
  for (rho in c(1, 0.7)) {            # random walk and AR(1)
    bd <- gmrf_bands(w, rho, mu = 0, v0 = v0, m0 = m0, N = N)
    expect_equal(bd$b[1], m0 / v0)
    expect_equal(bd$b[-1], rep(0, N - 1))
    # and it agrees with the dense expansion of the joint density
    ref <- dense_reference(w, rho, 0, v0, m0, N)
    expect_equal(bd$b, ref$b)
  }
})

test_that("b = 0 everywhere only when both mu = 0 and m0 = 0", {
  w <- rep(1, 4)
  bd <- gmrf_bands(w, rho = 1, mu = 0, v0 = 100, m0 = 0, N = 5)
  expect_equal(bd$b, rep(0, 5))
})

test_that("the z[1] full conditional uses the init-prior mean (mu = 0, RW)", {
  # Single-site full conditional of z[1]: mean = (b[1] - P[1,2] z[2]) / P[1,1].
  # With a tight N(5, 1e-4) init prior it must sit essentially at 5, which the
  # sampler-level test in test-mu-offset-init.R exercises end-to-end; here we
  # check the algebra directly.
  w <- rep(2, 4); v0 <- 1e-4; m0 <- 5
  bd <- gmrf_bands(w, rho = 1, mu = 0, v0 = v0, m0 = m0, N = 5)
  z2 <- 1.7
  cm <- (bd$b[1] - bd$off[1] * z2) / bd$diag[1]
  expect_lt(abs(cm - 5), 0.01)
})

test_that("a nonzero init_mean anchors z[1] under AR(1) dynamics too", {
  # End-to-end: with the fixed (non-stationary) diffuse init prior, a tight
  # N(5, 1e-3) prior must pin the first latent state under ar1 as well.
  sim <- simulate_dynamic_poisson(80, 0.2, log_rate0 = 2, rho = 0.8,
                                  mu = 0.4, seed = 11)
  fit <- fit_dynamic_model(sim$y, latent_dynamics = "ar1",
                           prior = dynamic_prior(init_mean = 5, init_var = 1e-3),
                           nsave = 300, nburn = 225, seed = 11)
  expect_lt(abs(mean(fit$draws$z0) - 5), 0.5)
})
