# Core sampler ----------------------------------------------------------------
#
# A single Metropolis-within-Gibbs sampler that handles all observation
# families, all zero treatments, and both latent dynamics. The latent state is
# an N x J matrix Z with one column per latent series. The Poisson and binomial
# families have a single latent series (J = 1, the log-rate / logit). The
# multinomial family with K categories has J = K - 1 additive-log-ratio (ALR)
# series, one per non-baseline category. Every column follows its own GMRF
# state evolution
#   z_{t,k} = mu_k + rho_k z_{t-1,k} + eps_{t,k},
# encoded with precision P_k = D_rho' diag(1/sig2_k) D_rho (D_rho being the
# generalised first-difference matrix). rho = 1 gives the random walk ("rw"),
# where mu is a drift; for AR(1) ("ar1") rho and the intercept mu are drawn
# jointly from their exact conjugate Gaussian full conditional (a Gibbs
# step), with rho truncated to the stationary region (-1, 1). Conjugacy holds
# because the first latent state carries a fixed, diffuse
# N(init_mean, init_var) prior (default N(0, 100)) under BOTH dynamics. Since
# that prior does not involve (mu, rho), the (mu, rho) full conditional is
# exactly the (truncated) Gaussian implied by the transition likelihood and
# the Gaussian priors, and no Metropolis-Hastings correction is needed.
#
# The columns are a priori independent and share NO parameters. Each has its
# own innovation state (sig2, and for "t"/"mixture"/"sv" the auxiliary
# quantities), its own (mu, rho), its own adaptive MH scales, and its own
# copy of the prior. They are coupled only through the likelihood. In the
# multinomial family the category probabilities at time t are the softmax of
# all K - 1 ALR states at t (with the baseline fixed at 0), so the
# log-likelihood of z_{t,k} depends on the other columns' values at t.
#
# Latent-state update (checkerboard sweep). Because P_k is tridiagonal, the
# full conditional of Z[i, k] depends on the rest of the model only through
# its neighbours Z[i - 1, k] and Z[i + 1, k] and, for the multinomial family,
# through the other columns at the same time point, which are held fixed
# while column k is updated. States with odd index are therefore conditionally
# independent given the states with even index (and vice versa), so all odd
# sites of a column can be updated simultaneously, followed by all even sites.
# Each site gets its own adaptive random-walk Metropolis step (or an exact
# Gibbs draw if it is missing). This is the same single-site Metropolis-within-Gibbs
# scheme as a sequential sweep, visited in red-black order, and it lets every
# half-sweep be computed with vectorised operations.
#
# State vector layout (length N = n + 1), per column:
#   z[1]        leading pad (the initial state; always missing -> prior only)
#   z[2..n+1]   latent states aligned with observations y[1..n]
#
# Forecasting is NOT part of the sampler. Future states carry no likelihood,
# so their joint posterior given the parameters and z_n is the forward
# transition law. Forecast draws are therefore produced after sampling by
# forward simulation from the stored parameter draws (see forecast.R). This
# targets exactly the same posterior predictive distribution as appending
# the future states to the latent vector and sampling them as missing values.

#' Log-likelihood of observations under the three observation families
#'
#' Vectorised over observations: `j` are observation indices, `zz` the values
#' of latent column `k` at those observations (all other columns are read
#' from `Z`, the padded `N x J` latent matrix, whose row `j + 1` is aligned
#' with observation `j`). Terms that do not depend on `zz` are kept, so the
#' values are exact log-probabilities.
#' @keywords internal
#' @noRd
obs_loglik <- function(family, j, k, zz, Z, Y, Ntot, offset) {
  switch(family,
    poisson = stats::dpois(Y[j, 1L], exp(offset[j, 1L] + zz), log = TRUE),
    binomial = stats::dbinom(Y[j, 1L], size = Ntot[j],
                             prob = stats::plogis(offset[j, 1L] + zz), log = TRUE),
    multinomial = {
      # sum_k y_k eta_k - N log(1 + sum_k exp(eta_k)) + log multinomial coef,
      # with eta_b = 0 for the baseline
      eta <- offset[j, , drop = FALSE] + Z[j + 1L, , drop = FALSE]
      eta[, k] <- offset[j, k] + zz
      Yj <- Y[j, , drop = FALSE]
      rowSums(Yj * eta) - Ntot[j] * row_log_sum_exp(cbind(0, eta)) +
        lgamma(Ntot[j] + 1) - rowSums(lgamma(Yj + 1)) - lgamma(Ntot[j] - rowSums(Yj) + 1)
    })
}

#' Internal sampler (not exported)
#'
#' @param y numeric vector (poisson/binomial) or n x K count matrix
#'   (multinomial, columns in the user's order).
#' @param baseline integer column index of the baseline category
#'   (multinomial only).
#' @param offset for the multinomial family an n x J matrix on the ALR scale
#'   (J = K - 1, columns aligned with the non-baseline categories in their
#'   original order); otherwise a vector.
#' @keywords internal
#' @noRd
dynamic_sampler <- function(y,
                            family,
                            trials = NULL,
                            baseline = NULL,
                            innovations = "gaussian",
                            zeros = "none",
                            latent_dynamics = "rw",
                            include_mu = FALSE,
                            prior = dynamic_prior(),
                            nsave = 4000,
                            nburn = 1000,
                            thin = 1,
                            offset = NULL,
                            verbose = FALSE) {

  is_pois  <- family == "poisson"
  is_binom <- family == "binomial"
  is_multi <- family == "multinomial"

  # observation layout ---------------------------------------------------------
  # For the multinomial family, Y holds the J = K - 1 non-baseline columns
  # (in their original order), yb the baseline counts and Ntot the row totals
  # (the known "trials" of the multinomial). For the univariate families
  # J = 1 and Y is y as a one-column matrix.
  if (is_multi) {
    y <- as.matrix(y)
    n  <- nrow(y)
    K  <- ncol(y)
    J  <- K - 1L
    nb_idx <- setdiff(seq_len(K), baseline)     # non-baseline columns
    Y    <- y[, nb_idx, drop = FALSE]
    yb   <- y[, baseline]
    Ntot <- rowSums(y)
    cat_names <- colnames(y)
  } else {
    n <- length(y)
    K <- 1L
    J <- 1L
    Y <- matrix(y, n, 1L)
    Ntot <- trials
  }
  N   <- n + 1L                # leading pad + n observations
  # `nsave` is the number of draws KEPT. We run `thin` sampling iterations per
  # kept draw, so the post-burn-in phase is nsave * thin iterations long.
  ntot <- nburn + nsave * thin
  n_store <- nsave

  # offsets: always held internally as an n x J matrix ------------------------
  if (is.null(offset)) offset <- 0
  if (is_multi) {
    if (length(offset) == 1L) offset <- matrix(offset, n, J)
    offset <- as.matrix(offset)
    if (!identical(dim(offset), c(n, J)))
      stop("`offset` must be a scalar or an n x (K - 1) matrix.", call. = FALSE)
  } else {
    if (length(offset) == 1L) offset <- rep(offset, n)
    if (length(offset) != n) stop("`offset` must have length 1 or n.", call. = FALSE)
    offset <- matrix(offset, n, 1L)
  }

  zero_inflated <- isTRUE(zeros == "inflated")
  is_ar1 <- isTRUE(latent_dynamics == "ar1")
  rho <- rep(if (is_ar1) 0.9 else 1, J)       # RW fixes rho = 1 (never sampled)
  mu  <- rep(0, J)                            # drift (rw) / intercept (ar1); 0 if disabled

  # initial values of the latent states ----------------------------------------
  if (is_pois) {
    Z_init <- log(Y + 0.5) - offset
  } else if (is_binom) {
    phat <- (Y + 0.5) / (Ntot + 1)
    Z_init <- stats::qlogis(phat) - offset
  } else {
    # ALR initial values: log((y_k + 0.5) / (y_b + 0.5)) - offset
    Z_init <- log(Y + 0.5) - log(yb + 0.5) - offset
  }

  # latent state and missingness ----------------------------------------------
  # Missingness is a property of the time index and is shared by all columns.
  Z <- rbind(Z_init[1L, , drop = FALSE], Z_init)
  if (is_multi) {
    # A row with zero total count carries no information about the shares:
    # its latent states are drawn from their GMRF full conditionals.
    M <- c(TRUE, Ntot == 0)
  } else {
    M <- switch(
      zeros,
      none     = c(TRUE, rep(FALSE, n)),   # leading pad only
      missing  = c(TRUE, y == 0),          # + zeros treated missing
      inflated = c(TRUE, y == 0)           # gate overwrites obs part
    )
  }

  # GMRF precision P = D_rho' diag(1/sig2) D_rho (rho = 1 -> random walk),
  # represented by its tridiagonal bands (Pdiag, Poff) rather than a dense
  # matrix; the drift/intercept mu enters only as a linear term b (a mean/offset
  # in the latent full conditionals), NOT as a change to P. The bands and b are
  # rebuilt in O(N) each iteration as the increment variances change. One
  # innovation state and one set of bands per latent column.
  inno <- lapply(seq_len(J), function(k) init_innovation_state(innovations, N - 1L, prior))

  # Proper initial-state prior on z[1], anchoring the otherwise-improper GMRF:
  # a fixed, diffuse N(init_mean, init_var) under BOTH dynamics. Keeping it
  # free of (mu, rho) is what makes step 2b conjugate (see the file header).
  v0 <- prior$init_var
  m0 <- prior$init_mean
  Pdiag   <- matrix(NA_real_, N, J)        # P[i, i]     per column
  Poff    <- matrix(NA_real_, N - 1L, J)   # P[i, i+1]   per column, symmetric
  drift_b <- matrix(NA_real_, N, J)        # linear term per column
  for (k in seq_len(J)) {
    bands <- gmrf_bands(1 / inno[[k]]$sig2, rho[k], mu[k], v0, m0, N)
    Pdiag[, k] <- bands$diag; Poff[, k] <- bands$off; drift_b[, k] <- bands$b
  }

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
  tau_z   <- matrix(0.05, N, J)
  pi_open <- 0.5                       # P(gate open) for zero inflation
  v <- if (is_multi) rep(1, n) else ifelse(y == 0, 0, 1)  # gate indicator (1 = open)

  # the two halves of the checkerboard sweep
  parity <- list(seq.int(1L, N, by = 2L), seq.int(2L, N, by = 2L))

  # storage --------------------------------------------------------------------
  # Conditional versus unconditional in-sample quantities (relevant only under
  # zero inflation; the pairs coincide exactly for zeros = "none"/"missing"):
  #   fitted      : E[y_t | draw]   under the FULL model, i.e. the gate-open
  #                 probability times the latent-implied mean (unconditional).
  #   fitted_open : latent-implied mean of y_t given the gate is open
  #                 (conditional on v_t = 1).
  #   yrep        : posterior predictive replicate of y_t under the FULL model,
  #                 including the zero-inflation gate (a fresh Bernoulli(pi_open)
  #                 gate is drawn per draw and observation, marginalising the
  #                 imputed gate states). Use THIS for posterior predictive
  #                 checks of a zero-inflated fit.
  #   yrep_open   : replicate drawn straight from the observation model, i.e.
  #                 conditional on the gate being open. Under zero inflation it
  #                 has systematically too few zeros for PPC of the full model.
  # Response-scale quantities have one slice per OBSERVED category (K slices,
  # baseline included, in the user's column order); latent quantities have
  # one slice per latent column (J). Univariate families are dropped back to
  # matrices/vectors on exit.
  Hm <- prior$mix_components
  z_store     <- array(NA_real_, c(n_store, n, J))   # observation-aligned states
  z0_store    <- matrix(NA_real_, n_store, J)        # initial (pad) state
  sig2_store  <- array(NA_real_, c(n_store, n, J))   # variance of the increment into z_t
  fitted_store      <- array(NA_real_, c(n_store, n, K))   # unconditional mean of y
  yrep_store        <- array(NA_real_, c(n_store, n, K))   # unconditional PP replicate
  if (is_multi) {
    prob_store <- array(NA_real_, c(n_store, n, K))         # category shares
    fitted_open_store <- yrep_open_store <- NULL
  } else {
    fitted_open_store <- matrix(NA_real_, n_store, n)   # gate-open (latent-implied) mean
    yrep_open_store   <- matrix(NA_real_, n_store, n)   # gate-open replicate
    prob_store <- NULL
  }
  v_store       <- if (zero_inflated) matrix(NA_real_, n_store, n) else NULL
  pi_open_store <- if (zero_inflated) numeric(n_store) else NULL
  rho_store   <- matrix(NA_real_, n_store, J)   # AR(1) coefficient (1 for the RW)
  mu_store    <- matrix(NA_real_, n_store, J)   # drift/intercept (0 if disabled)
  innov_var_store <- matrix(NA_real_, n_store, J)  # representative increment variance
  scale_store <- if (innovations != "sv") matrix(NA_real_, n_store, J) else NULL
  nu_store    <- if (innovations == "t") matrix(NA_real_, n_store, J) else NULL
  if (innovations == "mixture") {
    mixw_store <- array(NA_real_, c(n_store, Hm, J))  # weights, sorted by variance
    mixv_store <- array(NA_real_, c(n_store, Hm, J))  # component variances, ascending
  } else {
    mixw_store <- mixv_store <- NULL
  }
  if (innovations == "sv") {
    svmu_store  <- matrix(NA_real_, n_store, J)
    svphi_store <- matrix(NA_real_, n_store, J)
    svsig_store <- matrix(NA_real_, n_store, J)
  } else {
    svmu_store <- svphi_store <- svsig_store <- NULL
  }

  obs_rows <- 2:(n + 1L)                        # observation-aligned states
  tick <- make_progress(ntot, verbose)

  for (irep in seq_len(ntot)) {

    # --- 1. update latent states: checkerboard sweep of single-site steps ----
    # See the file header for why the odd and the even sites can each be
    # updated simultaneously. Within a half-sweep the neighbour terms use the
    # current values of the other half.
    for (k in seq_len(J)) {
      for (idx in parity) {
        zk <- Z[, k]
        # sum_{j != i} P[i, j] z[j] = P[i, i-1] z[i-1] + P[i, i+1] z[i+1]
        nb <- c(0, Poff[, k] * zk[-N])[idx] + c(Poff[, k] * zk[-1L], 0)[idx]
        cvar <- 1 / Pdiag[idx, k]
        cm   <- cvar * (drift_b[idx, k] - nb)
        # NOTE: drift_b is NOT zero when mu = 0. gmrf_bands() folds the
        # initial-state prior into the linear term, drift_b[1] = mu*(D_rho'w)[1]
        # + m0/v0, so the z[1] full conditional is pulled towards the prior
        # mean m0 = init_mean even without a drift/intercept.
        miss <- M[idx]

        # Missing states (leading pad, structural zero, zero-total row): the
        # full conditional is exactly N(cm, cvar) -- no likelihood term -- so
        # draw them directly (Gibbs).
        if (any(miss)) {
          im <- idx[miss]
          Z[im, k] <- stats::rnorm(length(im), cm[miss], sqrt(cvar[miss]))
        }

        # Observed states: adaptive random-walk Metropolis targeting a 0.44
        # acceptance rate (the univariate single-site optimum). The proposal
        # is symmetric, so the ratio is prior x likelihood only.
        if (!all(miss)) {
          io <- idx[!miss]
          cmo <- cm[!miss]
          sdo <- sqrt(cvar[!miss])
          j <- io - 1L                              # observation indices
          z_cur  <- Z[io, k]
          z_prop <- stats::rnorm(length(io), z_cur, sqrt(tau_z[io, k]))
          ll_diff <- obs_loglik(family, j, k, z_prop, Z, Y, Ntot, offset) -
                     obs_loglik(family, j, k, z_cur,  Z, Y, Ntot, offset)
          pr_diff <- stats::dnorm(z_prop, cmo, sdo, log = TRUE) -
                     stats::dnorm(z_cur,  cmo, sdo, log = TRUE)
          alpha <- pmin(1, exp(pr_diff + ll_diff))
          # A non-finite ratio (e.g. -Inf - -Inf from an overflowing linear
          # predictor) must count as a rejection; letting NaN reach the
          # adaptation below would freeze tau_z (and the site) permanently.
          alpha[!is.finite(alpha)] <- 0
          acc <- stats::runif(length(io)) < alpha
          Z[io[acc], k] <- z_prop[acc]
          tau_z[io, k] <- exp(log(tau_z[io, k]) + irep^(-0.6) * (alpha - 0.44))
        }
      }
    }

    # --- 2. update innovation variances and (mu, rho), column by column ------
    for (k in seq_len(J)) {
      zk <- Z[, k]
      dz <- zk[-1] - rho[k] * zk[-N] - mu[k]   # e_t = z_t - rho z_{t-1} - mu
      inno[[k]] <- update_innovations(innovations, inno[[k]], dz, prior, irep)

      # --- 2b. update (mu, rho): exact conjugate Gibbs draw --------------------
      # The state equation is linear in (mu, rho) and the fixed z[1] prior does
      # not involve them, so the conjugate Gaussian law is the full conditional
      # (see the file header): draw rho from its marginal truncated to (-1, 1),
      # then mu | rho exactly. Since the truncation involves rho only, the pair
      # is an exact sample from the joint truncated full conditional.
      #   ar1 (+ mu) : joint truncated-conjugate Gaussian Gibbs draw of (rho, mu)
      #   rw  + mu   : exact conjugate Gaussian draw of mu       (rho = 1 fixed)
      #   rw  only   : nothing                                   (rho = 1, mu = 0)
      w <- 1 / inno[[k]]$sig2

      if (is_ar1) {                               # AR(1) always includes intercept
        av <- zk[-N]; yv <- zk[-1]
        Lam0 <- diag(c(1 / prior$mu_sd^2, 1 / prior$ar_rho_sd^2))
        m0p  <- c(prior$mu_mean, prior$ar_rho_mean)
        XtWX <- matrix(c(sum(w), sum(w * av), sum(w * av), sum(w * av^2)), 2L, 2L)
        XtWy <- c(sum(w * yv), sum(w * av * yv))
        Vpost <- solve(Lam0 + XtWX)
        mpost <- as.vector(Vpost %*% (Lam0 %*% m0p + XtWy))
        rho[k] <- rtnorm(mpost[2L], sqrt(Vpost[2L, 2L]), -1, 1)
        mu_cm  <- mpost[1L] + Vpost[1L, 2L] / Vpost[2L, 2L] * (rho[k] - mpost[2L])
        mu_cv  <- Vpost[1L, 1L] - Vpost[1L, 2L]^2 / Vpost[2L, 2L]
        mu[k]  <- stats::rnorm(1, mu_cm, sqrt(max(mu_cv, 0)))
      } else if (include_mu) {                    # RW drift: rho = 1 fixed, exact
        prec_post <- 1 / prior$mu_sd^2 + sum(w)
        mean_post <- (prior$mu_mean / prior$mu_sd^2 + sum(w * (zk[-1] - zk[-N]))) / prec_post
        mu[k] <- stats::rnorm(1, mean_post, sqrt(1 / prec_post))
      }

      # rebuild the GMRF precision bands and drift linear term (O(N)) -----------
      bands <- gmrf_bands(w, rho[k], mu[k], v0, m0, N)
      Pdiag[, k] <- bands$diag; Poff[, k] <- bands$off; drift_b[, k] <- bands$b
    }

    # --- 3. zero inflation: gate indicators and gate-open probability ---------
    if (zero_inflated) {
      z_y <- Z[obs_rows, 1L]                       # latent states at observations
      v <- ifelse(y == 0, 0, 1)
      idx0 <- which(y == 0)
      if (length(idx0)) {
        # P(Y = 0 | gate open) under the observation model
        if (is_pois) {
          p0 <- stats::dpois(0, exp(offset[idx0, 1L] + z_y[idx0]))
        } else {
          p0 <- stats::dbinom(0, size = Ntot[idx0],
                              prob = stats::plogis(offset[idx0, 1L] + z_y[idx0]))
        }
        # P(gate open | y = 0) = pi * P0 / (pi * P0 + (1 - pi))
        post_open <- pi_open * p0 / (pi_open * p0 + (1 - pi_open))
        v[idx0] <- as.numeric(post_open > stats::runif(length(idx0)))
      }
      M <- c(TRUE, v == 0)                         # structural zeros are missing
      pi_open <- stats::rbeta(1, prior$zi_open_a + sum(v),
                              prior$zi_open_b + n - sum(v))
    }

    # --- 4. store -------------------------------------------------------------
    if (irep > nburn && ((irep - nburn) %% thin == 0)) {
      pos <- (irep - nburn) %/% thin
      z_store[pos, , ] <- Z[obs_rows, ]
      z0_store[pos, ]  <- Z[1L, ]
      for (k in seq_len(J)) sig2_store[pos, , k] <- inno[[k]]$sig2
      Eta <- offset + Z[obs_rows, , drop = FALSE]    # linear predictor incl. offset
      if (is_multi) {
        # shares P (n x K, user column order), expected counts N_t P, and a
        # multinomial replicate at each t
        P <- alr_to_prob(Eta, baseline)
        prob_store[pos, , ]   <- P
        fitted_store[pos, , ] <- Ntot * P
        yrep_store[pos, , ]   <- rmultinom_vec(Ntot, P)
      } else {
        eta <- Eta[, 1L]
        # gate-open (conditional) mean and replicate from the observation model
        if (is_pois) {
          m_open <- exp(eta)
          y_open <- stats::rpois(n, m_open)
        } else {
          m_open <- Ntot * stats::plogis(eta)
          y_open <- stats::rbinom(n, size = Ntot, prob = stats::plogis(eta))
        }
        fitted_open_store[pos, ] <- m_open
        yrep_open_store[pos, ]   <- y_open
        if (zero_inflated) {
          # unconditional versions include the zero-inflation gate: the mean is
          # scaled by pi_open, and the replicate applies a fresh Bernoulli(pi_open)
          # gate (marginalising the imputed gate states v), mirroring the gate
          # applied to the forecast draws.
          fitted_store[pos, , 1L] <- pi_open * m_open
          open_rep <- stats::runif(n) < pi_open
          yrep_store[pos, , 1L] <- ifelse(open_rep, y_open, 0)
          v_store[pos, ] <- v
          pi_open_store[pos] <- pi_open
        } else {
          fitted_store[pos, , 1L] <- m_open           # identical without a gate
          yrep_store[pos, , 1L]   <- y_open
        }
      }
      for (k in seq_len(J)) {
        ik <- inno[[k]]
        rho_store[pos, k] <- rho[k]
        mu_store[pos, k]  <- mu[k]
        if (!is.null(scale_store)) scale_store[pos, k] <- ik$scale
        if (!is.null(nu_store))    nu_store[pos, k]    <- ik$nu
        # representative increment variance for the innovation-SD summary:
        # the estimated constant variance (gaussian), the marginal variance
        # implied by the scale and mixing distribution (t, mixture), or the
        # average per-increment variance over the in-sample transitions (sv).
        innov_var_store[pos, k] <- if (innovations == "t") {
          ik$scale * (if (ik$nu > 2) ik$nu / (ik$nu - 2) else 3)
        } else if (innovations == "mixture") {
          ik$scale * sum(ik$eta * ik$sig2_h)
        } else if (innovations == "sv") {
          mean(ik$sig2[-1L])
        } else {                       # gaussian
          ik$scale
        }
        if (innovations == "mixture") {
          # components are exchangeable; report them ordered by variance
          cv <- ik$scale * ik$sig2_h
          o  <- order(cv)
          mixv_store[pos, , k] <- cv[o]
          mixw_store[pos, , k] <- ik$eta[o]
        }
        if (innovations == "sv") {
          svmu_store[pos, k]  <- ik$sv_startpara$mu
          svphi_store[pos, k] <- ik$sv_startpara$phi
          svsig_store[pos, k] <- ik$sv_startpara$sigma
        }
      }
    }
    tick(irep)
  }

  # output: multinomial keeps arrays (with category dimnames); the univariate
  # families are dropped back to the matrix/vector layout ---------------------
  if (is_multi) {
    lat_names <- cat_names[nb_idx]
    dimnames(z_store)    <- list(NULL, NULL, lat_names)
    dimnames(sig2_store) <- list(NULL, NULL, lat_names)
    dimnames(fitted_store) <- dimnames(yrep_store) <- dimnames(prob_store) <-
      list(NULL, NULL, cat_names)
    colnames(z0_store) <- colnames(rho_store) <- colnames(mu_store) <-
      colnames(innov_var_store) <- lat_names
    for (nm in c("scale_store", "nu_store", "svmu_store", "svphi_store", "svsig_store")) {
      m <- get(nm)
      if (!is.null(m)) { colnames(m) <- lat_names; assign(nm, m) }
    }
    if (!is.null(mixw_store)) {
      dimnames(mixw_store) <- dimnames(mixv_store) <- list(NULL, NULL, lat_names)
    }
  } else {
    # drop ONLY the trailing category axis (keeps draws x time matrices even
    # when nsave = 1)
    drop3 <- function(a) matrix(a, dim(a)[1L], dim(a)[2L])
    col1  <- function(m) if (is.null(m)) NULL else m[, 1L]
    z_store      <- drop3(z_store)
    sig2_store   <- drop3(sig2_store)
    fitted_store <- drop3(fitted_store)
    yrep_store   <- drop3(yrep_store)
    z0_store <- col1(z0_store)
    rho_store <- col1(rho_store); mu_store <- col1(mu_store)
    innov_var_store <- col1(innov_var_store)
    scale_store <- col1(scale_store); nu_store <- col1(nu_store)
    svmu_store <- col1(svmu_store); svphi_store <- col1(svphi_store)
    svsig_store <- col1(svsig_store)
    if (!is.null(mixw_store)) {
      mixw_store <- drop3(mixw_store); mixv_store <- drop3(mixv_store)
    }
    offset <- offset[, 1L]
  }

  list(
    draws = list(
      z           = z_store,            # observation-aligned latent states
      z0          = z0_store,           # initial (pad) state
      sig2        = sig2_store,         # variance of the increment into z_t
      fitted      = fitted_store,       # unconditional mean (includes ZI gate)
      fitted_open = fitted_open_store,  # mean conditional on gate open (NULL: multinomial)
      yrep        = yrep_store,         # unconditional posterior predictive
      yrep_open   = yrep_open_store,    # replicate conditional on gate open (NULL: multinomial)
      fitted_prob = prob_store,         # category shares (multinomial only)
      gate        = v_store,            # NULL unless zeros = "inflated"
      pi_open     = pi_open_store,      # NULL unless zeros = "inflated"
      rho         = rho_store,
      mu          = mu_store,
      innov_var   = innov_var_store,
      scale       = scale_store,        # NULL for "sv"
      nu          = nu_store,           # NULL unless "t"
      mix_weight  = mixw_store,         # NULL unless "mixture"
      mix_var     = mixv_store,         # NULL unless "mixture"
      sv_mu       = svmu_store,         # NULL unless "sv"
      sv_phi      = svphi_store,
      sv_sigma    = svsig_store
    ),
    offset = offset
  )
}
