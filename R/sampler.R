# Core sampler ----------------------------------------------------------------
#
# A single Metropolis-within-Gibbs sampler that handles both observation
# families, all zero treatments, and both latent dynamics. The latent log-rate
# / logit trajectory `z` follows a GMRF state evolution
# z_t = mu + rho z_{t-1} + eps_t, encoded with precision
# P = D_rho' diag(1/sig2) D_rho  (D_rho being the generalised first-difference
# matrix). rho = 1 gives the random walk ("rw"), where mu is a drift; for AR(1)
# ("ar1") rho and the intercept mu are drawn jointly from their exact conjugate
# Gaussian full conditional (a Gibbs step), with rho truncated to the
# stationary region (-1, 1). Conjugacy holds because the first latent state
# carries a fixed, diffuse N(init_mean, init_var) prior (default N(0, 100))
# under BOTH dynamics; since that prior does not involve (mu, rho), the
# (mu, rho) full conditional is exactly the (truncated) Gaussian implied by the
# transition likelihood and the Gaussian priors, and no Metropolis-Hastings
# correction is needed. Each latent state is updated by a single-site adaptive
# random-walk Metropolis step using its GMRF full conditional; the innovation
# variance is updated by one of the four innovation structures.
#
# State vector layout (length N = n + 1 + H, where H = horizon):
#   z[1]            leading pad         (always missing -> prior only)
#   z[2..n+1]       latent states aligned with observations y[1..n]
#   z[n+2..n+1+H]   forecast states     (always missing -> prior only)
#
# Forecasting is done by appending H states to the latent vector and sampling
# them as missing values, exactly like the initial value, an observed zero treated
# as missing, or a structural zero under zero inflation: they make no likelihood
# contribution and are updated from their (Gaussian) GMRF full conditionals. Forecast responses are
# drawn from the observation model at the sampled forecast states.

#' Internal sampler (not exported)
#' @keywords internal
#' @noRd
dynamic_sampler <- function(y,
                            family,
                            trials = NULL,
                            innovations = "gaussian",
                            zeros = "none",
                            latent_dynamics = "rw",
                            include_mu = FALSE,
                            prior = dynamic_prior(),
                            nsave = 4000,
                            nburn = 1000,
                            thin = 1,
                            horizon = 0L,
                            forecast_trials = NULL,
                            forecast_offset = NULL,
                            offset = NULL,
                            verbose = FALSE,
                            seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  n   <- length(y)
  H   <- as.integer(horizon)
  N   <- n + 1L + H            # leading pad + n observations + H forecast states
  # `nsave` is the number of draws KEPT. We run `thin` sampling iterations per
  # kept draw, so the post-burn-in phase is nsave * thin iterations long.
  ntot <- nburn + nsave * thin
  n_store <- nsave

  if (is.null(offset)) offset <- rep(0, n)
  if (length(offset) == 1L) offset <- rep(offset, n)
  if (length(offset) != n) stop("`offset` must have length 1 or n.", call. = FALSE)

  # forecast-period offsets / trials (length H) -------------------------------
  if (H >= 1L) {
    if (is.null(forecast_offset)) forecast_offset <- rep(0, H)
    if (length(forecast_offset) == 1L) forecast_offset <- rep(forecast_offset, H)
    if (length(forecast_offset) != H)
      stop("`forecast_offset` must have length 1 or horizon.", call. = FALSE)
    if (family == "binomial") {
      if (is.null(forecast_trials)) forecast_trials <- rep(trials[n], H)
      if (length(forecast_trials) == 1L) forecast_trials <- rep(forecast_trials, H)
      if (length(forecast_trials) != H)
        stop("`forecast_trials` must have length 1 or horizon.", call. = FALSE)
    }
  }

  is_pois <- family == "poisson"
  zero_inflated <- isTRUE(zeros == "inflated")
  is_ar1 <- isTRUE(latent_dynamics == "ar1")
  rho <- if (is_ar1) 0.9 else 1               # RW fixes rho = 1 (never sampled)
  mu  <- 0                                    # drift (rw) / intercept (ar1); 0 if disabled

  # link, inverse link and per-observation log-likelihood ----------------------
  # The linear predictor is offset[j] + z; the offset is a known per-observation
  # term (e.g. log-exposure for the Poisson log link).
  if (is_pois) {
    z_init <- log(y + 0.5) - offset
    loglik_i <- function(j, zz) stats::dpois(y[j], exp(offset[j] + zz), log = TRUE)
  } else {
    phat <- (y + 0.5) / (trials + 1)
    z_init <- stats::qlogis(phat) - offset
    loglik_i <- function(j, zz) stats::dbinom(y[j], size = trials[j],
                                              prob = stats::plogis(offset[j] + zz),
                                              log = TRUE)
  }

  # latent state and missingness ----------------------------------------------
  # The H forecast states are appended after the n observation-aligned states
  # and, like the leading pad, are always missing (no likelihood contribution).
  z <- c(z_init[1], z_init, rep(z_init[n], H))
  fc_missing <- rep(TRUE, H)
  M <- switch(
    zeros,
    none     = c(TRUE, rep(FALSE, n), fc_missing),  # leading pad + forecast
    missing  = c(TRUE, y == 0,        fc_missing),  # + zeros treated missing
    inflated = c(TRUE, y == 0,        fc_missing)   # gate overwrites obs part
  )

  # GMRF precision: P = D_rho' diag(1/sig2) D_rho (rho = 1 -> random walk).
  # GMRF precision P = D_rho' diag(1/sig2) D_rho (rho = 1 -> random walk),
  # represented by its tridiagonal bands (Pdiag, Poff) rather than a dense
  # matrix; the drift/intercept mu enters only as a linear term b (a mean/offset
  # in the latent full conditionals), NOT as a change to P. The bands and b are
  # rebuilt in O(N) each iteration as the increment variances change.
  inno <- init_innovation_state(innovations, N - 1L, prior)

  # Proper initial-state prior on z[1], anchoring the otherwise-improper GMRF:
  # a fixed, diffuse N(init_mean, init_var) under BOTH dynamics. Keeping it
  # free of (mu, rho) is what makes step 2b conjugate (see the file header).
  v0 <- prior$init_var
  m0 <- prior$init_mean
  bands <- gmrf_bands(1 / inno$sig2, rho, mu, v0, m0, N)
  Pdiag <- bands$diag       # P[i, i]      (length N)
  Poff  <- bands$off        # P[i, i+1]    (length N - 1), symmetric
  drift_b <- bands$b        # linear term  (length N)

  # truncated-normal draw via the inverse CDF (base R, no dependencies). The
  # clamp guards against the degenerate case in which essentially all proposal
  # mass lies outside (a, b) on one side, where pnorm(a) == pnorm(b) in floating
  # point and qnorm() returns an infinite value; the draw is then pinned just
  # inside the nearer bound instead of corrupting the chain.
  rtnorm <- function(m, sd, a, b) {
    r <- stats::qnorm(stats::runif(1, stats::pnorm(a, m, sd),
                                   stats::pnorm(b, m, sd)), m, sd)
    min(max(r, a + 1e-12), b - 1e-12)
  }

  # adaptive MH scales and zero-inflation parameters ---------------------------
  tau_z   <- rep(0.05, N)
  pi_open <- 0.5                       # P(gate open) for zero inflation
  v <- ifelse(y == 0, 0, 1)            # gate indicator (1 = open / sampling zero)

  # storage --------------------------------------------------------------------
  # Conditional versus unconditional in-sample quantities (relevant only under
  # zero inflation; the pairs coincide exactly for zeros = "none"/"missing"):
  #   fitted      : E[y_t | draw]   under the FULL model, i.e. the gate-open
  #                 probability times the latent-implied mean (unconditional).
  #   fitted_open : latent-implied mean of y_t given the gate is open
  #                 (conditional on v_t = 1); exp(eta) / trials * plogis(eta).
  #   yrep        : posterior predictive replicate of y_t under the FULL model,
  #                 including the zero-inflation gate (a fresh Bernoulli(pi_open)
  #                 gate is drawn per draw and observation, marginalising the
  #                 imputed gate states). Use THIS for posterior predictive
  #                 checks of a zero-inflated fit.
  #   yrep_open   : replicate drawn straight from the observation model, i.e.
  #                 conditional on the gate being open. Under zero inflation it
  #                 has systematically too few zeros for PPC of the full model.
  z_store     <- matrix(NA_real_, n_store, N)
  sig2_store  <- matrix(NA_real_, n_store, N - 1L)
  fitted_store      <- matrix(NA_real_, n_store, n)   # unconditional mean of y
  fitted_open_store <- matrix(NA_real_, n_store, n)   # gate-open (latent-implied) mean
  yrep_store        <- matrix(NA_real_, n_store, n)   # unconditional PP replicate
  yrep_open_store   <- matrix(NA_real_, n_store, n)   # gate-open replicate
  v_store     <- matrix(NA_real_, n_store, n)
  pi_open_store   <- numeric(n_store)
  nu_store    <- numeric(n_store)
  rho_store   <- numeric(n_store)       # AR(1) coefficient (1 for the RW)
  mu_store    <- numeric(n_store)       # drift/intercept (0 if disabled)
  innov_var_store <- numeric(n_store)   # representative increment variance

  has_fc <- H >= 1L
  fc_z_store <- if (has_fc) matrix(NA_real_, n_store, H) else NULL  # latent forecast
  fc_y_store <- if (has_fc) matrix(NA_real_, n_store, H) else NULL  # response forecast

  tick <- make_progress(ntot, verbose)

  for (irep in seq_len(ntot)) {

    # --- 1. update latent states by single-site adaptive RW Metropolis --------
    # The GMRF full conditional of z[i] depends only on its two neighbours
    # (P is tridiagonal), so the conditional mean uses P[i,i-1] z[i-1] and
    # P[i,i+1] z[i+1] only -- an O(1) computation per site (O(N) per sweep).
    for (i in seq_len(N)) {
      cvar <- 1 / Pdiag[i]
      nb <- 0                                     # sum_{j != i} P[i,j] z[j]
      if (i > 1L) nb <- nb + Poff[i - 1L] * z[i - 1L]   # P[i, i-1] = Poff[i-1]
      if (i < N)  nb <- nb + Poff[i]      * z[i + 1L]   # P[i, i+1] = Poff[i]
      cm   <- cvar * (drift_b[i] - nb)
      # NOTE: drift_b is NOT zero when mu = 0. gmrf_bands() folds the
      # initial-state prior into the linear term, drift_b[1] = mu*(D_rho'w)[1]
      # + m0/v0, so the z[1] full conditional is pulled towards the prior mean
      # m0 = init_mean even without a drift/intercept (see test-gmrf-bands.R).
      if (M[i]) {
        # Missing state (leading pad, structural zero, or forecast step): the
        # full conditional is exactly the Gaussian N(cm, cvar) -- no likelihood
        # term -- so draw it directly (Gibbs). This is exact and mixes far
        # better than a random-walk proposal, which matters especially for the
        # appended forecast states.
        z[i] <- stats::rnorm(1, cm, sqrt(cvar))
      } else {
        # Observed state: adaptive random-walk Metropolis targeting a 0.44
        # acceptance rate (the univariate single-site optimum).
        j <- i - 1L                               # observation index
        z_prop  <- stats::rnorm(1, z[i], sqrt(tau_z[i]))
        ll_diff <- loglik_i(j, z_prop) - loglik_i(j, z[i])
        pr_diff <- stats::dnorm(z_prop, cm, sqrt(cvar), log = TRUE) -
                   stats::dnorm(z[i],   cm, sqrt(cvar), log = TRUE)
        alpha <- min(1, exp(pr_diff + ll_diff))
        # A non-finite ratio (e.g. -Inf - -Inf from an overflowing linear
        # predictor) must count as a rejection; letting NaN reach the
        # adaptation below would freeze tau_z[i] (and the site) permanently.
        if (!is.finite(alpha)) alpha <- 0
        if (stats::runif(1) < alpha) z[i] <- z_prop
        tau_z[i] <- exp(log(tau_z[i]) + irep^(-0.6) * (alpha - 0.44))
      }
    }

    # --- 2. update innovation variance (one of four structures) ---------------
    # All N - 1 increments are treated identically; the forecast-period
    # increments are ordinary increments of the sampled latent path.
    dz <- z[-1] - rho * z[-N] - mu             # increments e_t = z_t - rho z_{t-1} - mu
    inno <- update_innovations(innovations, inno, dz, prior, irep)

    # --- 2b. update (mu, rho): exact conjugate Gibbs draw ----------------------
    # The state equation is linear in (mu, rho) and the fixed z[1] prior does
    # not involve them, so the conjugate Gaussian law is the full conditional
    # (see the file header): draw rho from its marginal truncated to (-1, 1),
    # then mu | rho exactly. Since the truncation involves rho only, the pair
    # is an exact sample from the joint truncated full conditional.
    #   ar1 (+ mu) : joint truncated-conjugate Gaussian Gibbs draw of (rho, mu)
    #   rw  + mu   : exact conjugate Gaussian draw of mu       (rho = 1 fixed)
    #   rw  only   : nothing                                   (rho = 1, mu = 0)
    w <- 1 / inno$sig2

    if (is_ar1) {                               # AR(1) always includes intercept
      av <- z[-N]; yv <- z[-1]
      Lam0 <- diag(c(1 / prior$mu_sd^2, 1 / prior$ar_rho_sd^2))
      m0p  <- c(prior$mu_mean, prior$ar_rho_mean)
      XtWX <- matrix(c(sum(w), sum(w * av), sum(w * av), sum(w * av^2)), 2L, 2L)
      XtWy <- c(sum(w * yv), sum(w * av * yv))
      Vpost <- solve(Lam0 + XtWX)
      mpost <- as.vector(Vpost %*% (Lam0 %*% m0p + XtWy))
      rho   <- rtnorm(mpost[2L], sqrt(Vpost[2L, 2L]), -1, 1)
      mu_cm <- mpost[1L] + Vpost[1L, 2L] / Vpost[2L, 2L] * (rho - mpost[2L])
      mu_cv <- Vpost[1L, 1L] - Vpost[1L, 2L]^2 / Vpost[2L, 2L]
      mu    <- stats::rnorm(1, mu_cm, sqrt(max(mu_cv, 0)))
    } else if (include_mu) {                    # RW drift: rho = 1 fixed, exact
      prec_post <- 1 / prior$mu_sd^2 + sum(w)
      mean_post <- (prior$mu_mean / prior$mu_sd^2 + sum(w * (z[-1] - z[-N]))) / prec_post
      mu <- stats::rnorm(1, mean_post, sqrt(1 / prec_post))
    }

    # rebuild the GMRF precision bands and drift linear term (O(N)) -------------
    bands <- gmrf_bands(w, rho, mu, v0, m0, N)
    Pdiag <- bands$diag
    Poff  <- bands$off
    drift_b <- bands$b

    # --- 3. zero inflation: gate indicators and gate-open probability ---------
    if (zero_inflated) {
      z_y <- z[2:(n + 1L)]                         # latent states at observations
      v <- ifelse(y == 0, 0, 1)
      idx0 <- which(y == 0)
      if (length(idx0)) {
        # P(Y = 0 | gate open) under the observation model
        if (is_pois) {
          p0 <- stats::dpois(0, exp(offset[idx0] + z_y[idx0]))
        } else {
          p0 <- stats::dbinom(0, size = trials[idx0],
                              prob = stats::plogis(offset[idx0] + z_y[idx0]))
        }
        # P(gate open | y = 0) = pi * P0 / (pi * P0 + (1 - pi))
        post_open <- pi_open * p0 / (pi_open * p0 + (1 - pi_open))
        v[idx0] <- as.numeric(post_open > stats::runif(length(idx0)))
      }
      M <- c(TRUE, v == 0, fc_missing)             # structural zeros + forecast missing
      pi_open <- stats::rbeta(1, prior$zi_open_a + sum(v),
                              prior$zi_open_b + n - sum(v))
    }

    # --- 4. store -------------------------------------------------------------
    if (irep > nburn && ((irep - nburn) %% thin == 0)) {
      pos <- (irep - nburn) %/% thin
      z_store[pos, ]    <- z
      sig2_store[pos, ] <- inno$sig2
      z_obs <- z[2:(n + 1L)]                        # observation-aligned states
      eta   <- offset + z_obs                       # linear predictor incl. offset
      # gate-open (conditional) mean and replicate from the observation model
      if (is_pois) {
        m_open <- exp(eta)
        y_open <- stats::rpois(n, m_open)
      } else {
        m_open <- trials * stats::plogis(eta)
        y_open <- stats::rbinom(n, size = trials, prob = stats::plogis(eta))
      }
      fitted_open_store[pos, ] <- m_open
      yrep_open_store[pos, ]   <- y_open
      if (zero_inflated) {
        # unconditional versions include the zero-inflation gate: the mean is
        # scaled by pi_open, and the replicate applies a fresh Bernoulli(pi_open)
        # gate (marginalising the imputed gate states v), mirroring the gate
        # applied to the forecast draws below.
        fitted_store[pos, ] <- pi_open * m_open
        open_rep <- stats::runif(n) < pi_open
        yrep_store[pos, ] <- ifelse(open_rep, y_open, 0)
      } else {
        fitted_store[pos, ] <- m_open               # identical without a gate
        yrep_store[pos, ]   <- y_open
      }
      v_store[pos, ]  <- v
      pi_open_store[pos]  <- if (zero_inflated) pi_open else 1
      nu_store[pos]   <- if (innovations == "t") inno$nu else NA_real_
      rho_store[pos]  <- rho
      mu_store[pos]   <- mu
      # representative increment variance for the innovation-SD summary:
      # the estimated constant variance (gaussian), the marginal variance
      # implied by the scale and mixing distribution (t, mixture), or the
      # average per-increment variance over the series (sv).
      if (innovations == "t") {
        innov_var_store[pos] <- inno$scale *
          (if (inno$nu > 2) inno$nu / (inno$nu - 2) else 3)
      } else if (innovations == "mixture") {
        innov_var_store[pos] <- inno$scale * sum(inno$eta * inno$sig2_h)
      } else if (innovations == "sv") {
        innov_var_store[pos] <- mean(inno$sig2)
      } else {                       # gaussian
        innov_var_store[pos] <- inno$scale
      }

      # --- forecast: the appended states z[n+2 .. N] are the sampled latent
      # forecast path; responses are drawn from the observation model at them.
      # Under zero inflation the gate IS applied to the response forecast, so
      # forecast_y is an UNCONDITIONAL predictive draw (consistent with `yrep`,
      # not with `yrep_open`); the latent forecast path forecast_z is gate-free.
      if (has_fc) {
        zf <- z[(n + 2L):N]
        fc_z_store[pos, ] <- zf
        eta_f <- forecast_offset + zf
        if (is_pois) {
          yf <- stats::rpois(H, exp(eta_f))
        } else {
          yf <- stats::rbinom(H, size = forecast_trials, prob = stats::plogis(eta_f))
        }
        if (zero_inflated) {
          closed <- stats::runif(H) > pi_open       # gate closed -> structural zero
          yf[closed] <- 0
        }
        fc_y_store[pos, ] <- yf
      }
    }
    tick(irep)
  }

  list(
    draws = list(
      z           = z_store,
      sig2        = sig2_store,
      fitted      = fitted_store,       # unconditional mean (includes ZI gate)
      fitted_open = fitted_open_store,  # mean conditional on gate open
      yrep        = yrep_store,         # unconditional posterior predictive
      yrep_open   = yrep_open_store,    # replicate conditional on gate open
      gate        = v_store,
      pi_open     = pi_open_store,
      nu          = nu_store,
      rho         = rho_store,
      mu          = mu_store,
      innov_var   = innov_var_store,
      forecast_z  = fc_z_store,        # NULL when horizon = 0
      forecast_y  = fc_y_store         # NULL when horizon = 0; includes ZI gate
    ),
    offset          = offset,
    forecast_offset = if (has_fc) forecast_offset else NULL,
    forecast_trials = if (has_fc && family == "binomial") forecast_trials else NULL
  )
}
