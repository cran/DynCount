# Innovation structures -------------------------------------------------------
#
# These functions update the per-increment innovation variance vector `sig2`
# (length n) given the current latent increments `dz`. Each updater takes and
# returns a small mutable `state` list so that the auxiliary quantities
# (Student-t mixing weights, degrees of freedom, mixture allocations, SV
# parameters) persist across MCMC iterations.

#' Initialise innovation-structure state
#' @keywords internal
#' @noRd
init_innovation_state <- function(innovations, n_incr, prior) {
  state <- list(
    sig2  = rep(1, n_incr),       # per-increment variance (length n)
    omega = rep(1, n_incr),       # Student-t / mixture mixing weights
    scale = 1,                    # overall variance/scale (gaussian, t, mixture)
    nu    = prior$df_min + prior$df_mean_excess,  # current t degrees of freedom
    tau_nu = 0.01,                # RW-MH scale for log(nu)
    sig2_h = rep(1, prior$mix_components),         # mixture component variances
    eta    = rep(1 / prior$mix_components, prior$mix_components)  # mixture weights
  )
  if (innovations == "sv") {
    if (!requireNamespace("stochvol", quietly = TRUE)) {
      stop("`innovations = \"sv\"` requires the 'stochvol' package.", call. = FALSE)
    }
    state$sv_prior <- prior$sv_prior %||% stochvol::specify_priors()
    state$sv_startpara <- list(mu = 0, phi = 0.9, sigma = 0.1,
                               nu = Inf, rho = 0, beta = NA, latent0 = 0)
    state$sv_startlatent <- rep(0, n_incr)
  }
  state
}

#' Update the innovation variance, dispatching on the structure
#'
#' @param innovations one of "gaussian", "t", "mixture", "sv"
#' @param state innovation state list (see [init_innovation_state()])
#' @param dz numeric vector of latent increments (length n)
#' @param prior a [dynamic_prior()] object
#' @param iter current MCMC iteration (for the t-df adaptive MH step)
#' @return updated `state` (with refreshed `state$sig2` of length n)
#' @keywords internal
#' @noRd
update_innovations <- function(innovations, state, dz, prior, iter) {
  a0 <- prior$var_shape
  b0 <- prior$var_rate
  n  <- length(dz)

  if (innovations == "gaussian") {
    # Constant-variance Gaussian increments: a single variance shared across t.
    v <- 1 / stats::rgamma(1, a0 + 0.5 * n, b0 + 0.5 * sum(dz^2))
    state$scale <- v
    state$sig2  <- rep(v, n)

  } else if (innovations == "t") {
    # Student-t scale mixture. `scale` is the overall variance; `omega` are the
    # per-increment Gamma mixing weights; `nu` (degrees of freedom) is updated
    # by an adaptive random-walk Metropolis step on the log scale.
    scale <- 1 / stats::rgamma(1, a0 + 0.5 * n, b0 + 0.5 * sum(state$omega * dz^2))
    state$omega <- stats::rgamma(n, 0.5 + state$nu / 2,
                                 state$nu / 2 + 0.5 * (dz^2 / scale))
    state$scale <- scale
    state$sig2  <- scale / state$omega

    nu_prop <- exp(stats::rnorm(1, log(state$nu), sqrt(state$tau_nu)))
    ll_diff <- sum(stats::dgamma(state$omega, shape = nu_prop / 2,
                                 rate = nu_prop / 2, log = TRUE)) -
               sum(stats::dgamma(state$omega, shape = state$nu / 2,
                                 rate = state$nu / 2, log = TRUE))
    # prior: nu - df_min ~ Exp(rate = 1 / df_mean_excess)
    rate <- 1 / prior$df_mean_excess
    pr_diff <- stats::dexp(nu_prop - prior$df_min, rate = rate, log = TRUE) -
               stats::dexp(state$nu - prior$df_min, rate = rate, log = TRUE)
    jacobian <- log(nu_prop / state$nu)        # log-scale RW proposal Jacobian
    alpha <- min(1, exp(jacobian + pr_diff + ll_diff))
    # treat a non-finite ratio as a rejection so NaN never reaches the
    # adaptation of tau_nu (which would freeze the nu update permanently)
    if (!is.finite(alpha)) alpha <- 0
    if (stats::runif(1) < alpha) state$nu <- nu_prop
    # target a 0.44 acceptance rate (univariate single-site RW-MH optimum)
    state$tau_nu <- exp(log(state$tau_nu) + iter^(-0.6) * (alpha - 0.44))

  } else if (innovations == "mixture") {
    # Finite scale mixture of normals (K components) using the Gumbel-max trick
    # for the categorical allocation, with Dirichlet weights.
    H <- prior$mix_components
    scale <- 1 / stats::rgamma(1, a0 + 0.5 * n, b0 + 0.5 * sum(state$omega * dz^2))
    pp <- sapply(seq_len(H), function(h)
      stats::dnorm(dz, 0, sqrt(scale * state$sig2_h[h]), log = TRUE))
    pp <- t(log(c(state$eta)) + t(pp))
    vmix <- max.col(pp + matrix(rgumbel0(n * H), n, H))
    counts <- tabulate(vmix, H)
    state$eta <- rdirichlet1(prior$mix_concentration + counts)
    state$sig2_h <- sapply(seq_len(H), function(h)
      1 / stats::rgamma(1, a0 + 0.5 * counts[h],
                        b0 + 0.5 * sum((dz[vmix == h])^2 / scale)))
    state$scale <- scale
    state$omega <- 1 / state$sig2_h[vmix]
    state$sig2 <- scale / state$omega

  } else if (innovations == "sv") {
    # Stochastic volatility: log-variance follows an AR(1), fit with stochvol.
    res <- stochvol::svsample_fast_cpp(
      dz, startpara = state$sv_startpara,
      startlatent = log(state$sig2), priorspec = state$sv_prior)
    state$sv_startpara[c("mu", "phi", "sigma")] <-
      as.list(res$para[, c("mu", "phi", "sigma")])
    state$sv_startlatent <- drop(res$latent)
    state$sig2 <- as.numeric(exp(state$sv_startlatent))

  } else {
    stop("Unknown innovation structure: ", innovations, call. = FALSE)
  }

  state
}

# Names of the four supported innovation structures.
.innovation_choices <- c("gaussian", "t", "mixture", "sv")

#' @keywords internal
#' @noRd
normalise_innovations <- function(innovations) {
  # If the caller passed the full default vector, take the first element.
  if (length(innovations) > 1L) innovations <- innovations[1L]
  match.arg(innovations, .innovation_choices)
}
