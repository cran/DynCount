# User-facing fitting function ------------------------------------------------

#' Fit a Bayesian dynamic count / binomial time-series model
#'
#' Fits a GMRF state-space model in which a latent trajectory \eqn{z_t} evolves
#' as either a first-order random walk (`latent_dynamics = "rw"`, the default)
#' or a stationary AR(1) process (`latent_dynamics = "ar1"`), and the
#' observations are linked to it through a Poisson (log link) or binomial
#' (logit link) observation model. See [DynCount-package] for an overview.
#'
#' @details
#' \strong{Zero handling.} Under `zeros = "inflated"` a latent gate decides,
#' for each observed zero, whether it is structural (gate closed) or an
#' ordinary sampling zero produced by the Poisson/binomial process (gate
#' open); the gate-open probability is a single time-constant parameter. See
#' [structural_zero_prob()].
#'
#' \strong{Forecasting.} Forecasts are produced \emph{during} sampling: set
#' `horizon = H` (a positive integer) and `H` extra latent states are appended
#' to the state vector and sampled as missing values inside the MCMC. A response is
#' drawn from the observation model at each sampled forecast state. The stored
#' forecast draws are retrieved with [forecast.dynamic_fit()]. With the default
#' `horizon = 0` no forecast is produced.
#'
#' @param y Numeric vector of non-negative integer observations. For the
#'   Poisson family these are counts; for the binomial family these are the
#'   numbers of successes.
#' @param family Observation model, `"poisson"` (default) or `"binomial"`.
#' @param trials For `family = "binomial"`, the number of trials. Either a
#'   single number (recycled) or a vector the same length as `y`. Each `y` must
#'   not exceed its number of trials. Ignored for the Poisson family.
#' @param innovations Distribution of the latent increments, one of
#'   `"gaussian"` (default), `"t"`, `"mixture"`, `"sv"`. See
#'   [DynCount-package] for details. `"sv"` requires the \pkg{stochvol} package.
#' @param latent_dynamics Latent state evolution: `"rw"` (default; a
#'   first-order random walk, \eqn{\rho = 1} fixed) or `"ar1"` (a stationary
#'   AR(1) with \eqn{\rho} sampled on \eqn{(-1, 1)}). `"ar1"` always includes
#'   an intercept (`include_mu` is forced to `TRUE`).
#' @param include_mu Logical; include a scalar \eqn{\mu} in the state equation
#'   \eqn{z_t = \mu + \rho z_{t-1} + \varepsilon_t}. It is a \strong{drift} under
#'   `latent_dynamics = "rw"` (\eqn{\rho = 1}) and an \strong{intercept} under
#'   `"ar1"`. With `FALSE` (default) \eqn{\mu = 0} and is not sampled; with
#'   `TRUE` it has the Gaussian prior in [dynamic_prior()]. When
#'   `latent_dynamics = "ar1"` this is forced to `TRUE` regardless of the value
#'   supplied: the intercept gives the process a non-zero stationary mean
#'   \eqn{\mu / (1 - \rho)}, without which the zero-mean stationary assumption is
#'   rarely appropriate for a log-rate/logit series.
#' @param zeros Zero handling for both families: `"none"` (default),
#'   `"inflated"` (time-constant zero inflation; see
#'   [structural_zero_prob()]), or `"missing"` (observed zeros treated as
#'   missing data). `zero_inflation = TRUE` is shorthand for
#'   `zeros = "inflated"`. See Details.
#' @param zero_inflation Logical convenience flag. If `TRUE`, sets
#'   `zeros = "inflated"`. Ignored if `zeros` is supplied explicitly.
#' @param prior A [dynamic_prior()] object giving the prior hyperparameters.
#' @param nsave Number of posterior draws to keep. Default `4000`.
#' @param nburn Number of burn-in iterations. Default `1000`.
#' @param thin Thinning interval: one draw is kept every `thin` iterations, so
#'   the sampler runs `nburn + nsave * thin` iterations in total and retains
#'   `nsave` draws. Default `1`.
#' @param horizon Forecast horizon `H` (a non-negative integer). When
#'   `H >= 1`, forecasts are produced inside the sampler and stored for
#'   retrieval with [forecast.dynamic_fit()]; when `H = 0` (default) no forecast
#'   is produced.
#' @param forecast_trials For the binomial family, the number of trials in the
#'   forecast period (length 1, recycled, or length `horizon`). If omitted it
#'   defaults to the last observed trials count. Ignored when `horizon = 0`.
#' @param forecast_offset Known per-period offset over the forecast horizon
#'   (length 1, recycled, or length `horizon`). Defaults to `0`. Ignored when
#'   `horizon = 0`.
#' @param offset Optional known per-observation offset on the linear-predictor
#'   scale (length 1 or `n`). For the Poisson family this is a log-exposure
#'   term, so the mean is \eqn{\exp(\mathrm{offset}_t + z_t)}; for the binomial
#'   family it shifts the logit. Default `NULL` (no offset).
#' @param verbose Logical; print a progress bar. Default `FALSE`.
#' @param seed Optional integer random seed for reproducibility.
#'
#' @return An object of class `"dynamic_fit"`: a list with three elements.
#'   \describe{
#'     \item{`draws`}{Posterior draws only -- a list of matrices/vectors with
#'       one row (or element) per stored draw:
#'       \describe{
#'         \item{`z`, `sig2`}{Latent states and increment variances.}
#'         \item{`fitted`, `yrep`}{\strong{Unconditional} fitted means and
#'           posterior predictive replicates; under `zeros = "inflated"` they
#'           include the zero-inflation gate.}
#'         \item{`fitted_open`, `yrep_open`}{Their
#'           \strong{conditional-on-gate-open} counterparts: the latent-implied
#'           mean, and a replicate drawn straight from the observation model.}
#'         \item{`gate`, `pi_open`}{Zero-inflation gate indicators and
#'           gate-open probability. For fits without zero inflation these are
#'           stored as placeholders: `gate` is the constant indicator
#'           `y > 0` and `pi_open` is `1`.}
#'         \item{`rho`, `mu`, `nu`, `innov_var`}{AR(1) coefficient,
#'           drift/intercept, Student-t degrees of freedom, and a
#'           representative innovation variance whose definition depends on
#'           `innovations`: the estimated constant increment variance
#'           \eqn{\sigma^2} for `"gaussian"`; the marginal increment variance
#'           for `"t"` (\eqn{\sigma^2 \nu / (\nu - 2)}) and `"mixture"`
#'           (\eqn{\sigma^2 \sum_h \eta_h \sigma_h^2}); and the average of
#'           the per-increment variances \eqn{\exp(h_t)} over the series for
#'           `"sv"`. For `"t"`, if a draw has \eqn{\nu \le 2} (possible only
#'           when `df_min <= 2`), the undefined variance factor is replaced
#'           by the convention `3`.}
#'         \item{`forecast_z`, `forecast_y`}{Latent and response forecasts
#'           when `horizon >= 1` (`NULL` otherwise); `forecast_y` is
#'           unconditional, i.e. includes the gate.}
#'       }
#'       Use `yrep`, not `yrep_open`, for posterior predictive checks; without
#'       zero inflation the conditional and unconditional pairs are identical.
#'       See [predict.dynamic_fit()].}
#'     \item{`data`}{The observed inputs: `y`, `trials`, the series length `n`,
#'       and the resolved per-observation `offset`.}
#'     \item{`spec`}{The model and MCMC specification: `family`, `innovations`,
#'       `latent_dynamics`, `include_mu`, `zeros`, `prior`, `nsave`, `nburn`,
#'       `thin`, `horizon`, and the fixed forecast-period inputs
#'       `forecast_offset` / `forecast_trials` (`NULL` when not applicable).}
#'   }
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2)
#' fit <- fit_dynamic_model(sim$y, family = "poisson", nsave = 300, nburn = 200)
#' summary(fit)
#'
#' @seealso [dynamic_prior()], [forecast.dynamic_fit()], [plot_fitted()],
#'   [structural_zero_prob()]
#' @export
fit_dynamic_model <- function(y,
                              family = c("poisson", "binomial"),
                              trials = NULL,
                              innovations = c("gaussian", "t", "mixture", "sv"),
                              latent_dynamics = c("rw", "ar1"),
                              include_mu = FALSE,
                              zeros = c("none", "inflated", "missing"),
                              zero_inflation = FALSE,
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

  family      <- match.arg(family)
  latent_dynamics <- match.arg(latent_dynamics)
  stopifnot(is.logical(include_mu), length(include_mu) == 1L)
  # AR(1) always carries an intercept: without it the latent process is assumed
  # to have a stationary mean of zero, which is rarely appropriate for a
  # log-rate/logit series and forces rho -> 1. Enable it silently.
  if (identical(latent_dynamics, "ar1")) include_mu <- TRUE
  innovations <- normalise_innovations(innovations)

  # resolve zero handling -----------------------------------------------------
  if (missing(zeros) && isTRUE(zero_inflation)) {
    zeros <- "inflated"
  } else {
    zeros <- match.arg(zeros)
  }
  if (!inherits(prior, "dynamic_prior")) {
    stop("`prior` must be created with dynamic_prior().", call. = FALSE)
  }

  # input validation ----------------------------------------------------------
  if (!is.numeric(y) || anyNA(y) || any(y < 0) || any(y != round(y))) {
    stop("`y` must be a vector of non-negative integers.", call. = FALSE)
  }
  y <- as.numeric(y)
  n <- length(y)
  if (n < 3) stop("`y` must have at least 3 observations.", call. = FALSE)

  if (family == "binomial") {
    if (is.null(trials)) stop("`trials` is required for `family = \"binomial\"`.",
                              call. = FALSE)
    if (!is.numeric(trials) || anyNA(trials)) {
      stop("`trials` must be a numeric vector without missing values.",
           call. = FALSE)
    }
    if (length(trials) == 1L) trials <- rep(trials, n)
    if (length(trials) != n) {
      stop("`trials` must be length 1 or length(y).", call. = FALSE)
    }
    if (any(trials <= 0) || any(trials != round(trials))) {
      stop("`trials` must be positive integers.", call. = FALSE)
    }
    if (any(y > trials)) stop("Each `y` must be <= its number of `trials`.",
                              call. = FALSE)
    trials <- as.numeric(trials)
  }

  nsave <- check_count(nsave, "nsave")
  nburn <- check_count0(nburn, "nburn")
  thin  <- check_count(thin, "thin")

  # horizon: a single non-negative integer (0 = no forecast) ------------------
  if (length(horizon) != 1L || !is.finite(horizon) || horizon < 0 ||
      horizon != round(horizon)) {
    stop("`horizon` must be a single non-negative integer.", call. = FALSE)
  }
  horizon <- as.integer(horizon)

  # validate offsets: known, fixed inputs must be finite numerics --------------
  if (!is.null(offset) &&
      (!is.numeric(offset) || anyNA(offset) || any(!is.finite(offset)))) {
    stop("`offset` must be a finite numeric vector.", call. = FALSE)
  }
  if (!is.null(forecast_offset) &&
      (!is.numeric(forecast_offset) || anyNA(forecast_offset) ||
       any(!is.finite(forecast_offset)))) {
    stop("`forecast_offset` must be a finite numeric vector.", call. = FALSE)
  }

  # validate forecast_trials if given -----------------------------------------
  if (family == "binomial" && !is.null(forecast_trials)) {
    if (!is.numeric(forecast_trials) || any(forecast_trials <= 0) ||
        any(forecast_trials != round(forecast_trials))) {
      stop("`forecast_trials` must be positive integers.", call. = FALSE)
    }
    forecast_trials <- as.numeric(forecast_trials)
  }

  out <- dynamic_sampler(
    y = y, family = family, trials = trials, innovations = innovations,
    zeros = zeros, latent_dynamics = latent_dynamics, include_mu = include_mu,
    prior = prior, nsave = nsave, nburn = nburn, thin = thin,
    horizon = horizon, forecast_trials = forecast_trials,
    forecast_offset = forecast_offset,
    offset = offset, verbose = verbose, seed = seed
  )

  structure(
    list(
      draws = out$draws,
      data  = list(y = y, trials = trials, n = n, offset = out$offset),
      spec  = list(family = family, innovations = innovations,
                   latent_dynamics = latent_dynamics, include_mu = include_mu,
                   zeros = zeros, prior = prior, nsave = nsave, nburn = nburn,
                   thin = thin, horizon = horizon,
                   forecast_offset = out$forecast_offset,
                   forecast_trials = out$forecast_trials)
    ),
    class = "dynamic_fit"
  )
}

#' @export
print.dynamic_fit <- function(x, ...) {
  sp <- x$spec
  cat("<dynamic_fit>\n")
  cat(sprintf("  family      : %s\n", sp$family))
  cat(sprintf("  dynamics    : %s%s\n", sp$latent_dynamics,
              if (identical(sp$latent_dynamics, "ar1"))
                sprintf("  (rho post. mean = %.3f)", mean(x$draws$rho)) else "  (rho = 1)"))
  if (isTRUE(sp$include_mu)) {
    cat(sprintf("  %-12s: post. mean = %.3f\n",
                if (identical(sp$latent_dynamics, "ar1")) "intercept mu" else "drift mu",
                mean(x$draws$mu)))
  }
  cat(sprintf("  innovations : %s\n", sp$innovations))
  cat(sprintf("  zeros       : %s\n", sp$zeros))
  cat(sprintf("  observations: %d  (zeros: %d)\n",
              x$data$n, sum(x$data$y == 0)))
  cat(sprintf("  draws kept  : %d\n", nrow(x$draws$z)))
  if (sp$horizon >= 1L) {
    cat(sprintf("  forecast    : horizon %d (stored)\n", sp$horizon))
  }
  invisible(x)
}
