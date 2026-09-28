# User-facing fitting function ------------------------------------------------

#' Fit a Bayesian dynamic count / binomial / multinomial time-series model
#'
#' Fits a GMRF state-space model in which a latent trajectory \eqn{z_t} evolves
#' as either a first-order random walk (`latent_dynamics = "rw"`, the default)
#' or a stationary AR(1) process (`latent_dynamics = "ar1"`), and the
#' observations are linked to it through a Poisson (log link), binomial
#' (logit link) or multinomial (additive-log-ratio link) observation model.
#' See [DynCount-package] for an overview.
#'
#' @details
#' \strong{Estimation.} The model is estimated by Metropolis-within-Gibbs
#' MCMC. The latent states are updated one site at a time by adaptive
#' random-walk Metropolis steps that use their Gaussian Markov random field
#' (GMRF) full conditionals. States without an observation (the initial
#' state, structural zeros, zero-total rows) are drawn exactly from their
#' Gaussian full conditionals. The
#' innovation parameters, the drift/intercept \eqn{\mu} and the AR(1)
#' coefficient \eqn{\rho} are updated by Gibbs steps (with a Metropolis step
#' for the Student-t degrees of freedom).
#'
#' \strong{Multinomial family.} With `family = "multinomial"`, `y` is an
#' \eqn{n \times K} matrix of category counts (rows = time, columns =
#' categories, \eqn{K \ge 2}); the row totals \eqn{N_t} are treated as known
#' trials. One category \eqn{b} is the \strong{baseline}, and the remaining
#' \eqn{K - 1} categories each get their own latent additive-log-ratio (ALR)
#' series \eqn{z_{t,k} = \log(p_{t,k} / p_{t,b})}, so that
#' \deqn{y_t \sim \mathrm{Multinomial}(N_t, p_t), \qquad
#'       p_{t,k} = \frac{\exp(o_{t,k} + z_{t,k})}{1 + \sum_{j \ne b}
#'       \exp(o_{t,j} + z_{t,j})}, \qquad p_{t,b} = \frac{1}{1 + \sum_{j \ne b}
#'       \exp(o_{t,j} + z_{t,j})},}
#' with known offsets \eqn{o_{t,k}} (zero by default). The chosen latent
#' dynamics, innovation structure and drift/intercept setting apply to every
#' ALR series, but \emph{no parameters are shared across categories}. Each
#' series has its own innovation variance (and, where relevant, degrees of
#' freedom, mixture components or volatility path), its own \eqn{\rho} and
#' \eqn{\mu}, and its own copy of `prior`. The series are coupled only
#' through the multinomial likelihood, and each series is updated with the
#' other categories held at their current values. A row with
#' \eqn{N_t = 0} carries no information about the shares and is handled like
#' a missing observation. Zero inflation is not available for this family
#' (`zeros` must be `"none"`). Running time grows linearly in \eqn{K - 1}.
#'
#' \strong{Zero handling.} Under `zeros = "inflated"` a latent gate decides,
#' for each observed zero, whether it is structural (gate closed) or an
#' ordinary sampling zero produced by the Poisson/binomial process (gate
#' open); the gate-open probability is a single parameter that is constant
#' over time. See [structural_zero_prob()]. Under `zeros = "missing"` the
#' observed zeros are treated as missing values.
#'
#' \strong{Forecasting.} Forecasts are obtained by forward simulation from the
#' posterior draws, so they can be computed after fitting for any horizon with
#' [forecast.dynamic_fit()]. Setting `horizon = H` additionally simulates an
#' `H`-step forecast right after sampling and stores it in the fit (under the
#' same `seed`).
#'
#' \strong{Reproducibility.} With a `seed`, the sampler and the fit-time
#' forecast run with that seed, and the previous state of the global random
#' number generator is restored afterwards.
#'
#' @param y For the Poisson and binomial families a numeric vector of
#'   non-negative integer observations (counts, or numbers of successes). For
#'   the multinomial family an \eqn{n \times K} matrix (or data frame) of
#'   non-negative integer category counts with one row per time point; its
#'   column names, if any, are used as category labels.
#' @param family Observation model, `"poisson"` (default), `"binomial"` or
#'   `"multinomial"`.
#' @param trials For `family = "binomial"`, the number of trials. Either a
#'   single number (recycled) or a vector the same length as `y`. Each `y` must
#'   not exceed its number of trials. For the multinomial family the row
#'   totals of `y` are used and this argument must be left `NULL`; for the
#'   Poisson family it is ignored with a warning.
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
#' @param zeros Zero handling for the Poisson and binomial families: `"none"`
#'   (default), `"inflated"` (zero inflation with a time-constant gate-open
#'   probability; see [structural_zero_prob()]), or `"missing"` (observed zeros
#'   treated as missing data). Must be `"none"` for the multinomial family. See
#'   Details.
#' @param zero_inflation A single `TRUE` or `FALSE`. `TRUE` is shorthand for
#'   `zeros = "inflated"`; combining it with a different explicit `zeros` is an
#'   error. Unlike the `zero_inflation` argument of
#'   [simulate_dynamic_poisson()] and [simulate_dynamic_binomial()], it is not
#'   a probability.
#' @param prior A [dynamic_prior()] object giving the prior hyperparameters.
#'   For the multinomial family the same prior is applied independently to
#'   every non-baseline category.
#' @param nsave Number of posterior draws to keep. Default `4000`.
#' @param nburn Number of burn-in iterations. Default `1000`.
#' @param thin Thinning interval: one draw is kept every `thin` iterations, so
#'   the sampler runs `nburn + nsave * thin` iterations in total and retains
#'   `nsave` draws. Default `1`.
#' @param horizon Forecast horizon `H` (a non-negative integer). When
#'   `H >= 1`, an `H`-step forecast is simulated after sampling and stored for
#'   retrieval with [forecast.dynamic_fit()]; when `H = 0` (default) nothing is
#'   stored, and forecasts can still be computed later with
#'   `forecast(fit, horizon = H)`.
#' @param forecast_trials For the binomial and multinomial families, the number
#'   of trials (binomial) or the total count per period (multinomial) over the
#'   forecast horizon (length 1, recycled, or length `horizon`). If omitted it
#'   defaults to the last observed number of trials (binomial) or the last
#'   non-zero row total (multinomial). Used only when `horizon >= 1`.
#' @param forecast_offset Known offset over the forecast horizon. For the
#'   Poisson and binomial families a vector of length 1 (recycled) or
#'   `horizon`; for the multinomial family a scalar, a length \eqn{K - 1}
#'   vector (one constant per non-baseline category) or an
#'   \eqn{H \times (K - 1)} matrix on the ALR scale (columns aligned with the
#'   non-baseline categories in their original order). Defaults to `0`, with a
#'   warning if `offset` is non-zero. Used only when `horizon >= 1`.
#' @param offset Optional known per-observation offset on the linear-predictor
#'   scale. For the Poisson family (length 1 or `n`) this is a log-exposure
#'   term, so the mean is \eqn{\exp(\mathrm{offset}_t + z_t)}; for the binomial
#'   family it shifts the logit. For the multinomial family it is a scalar, a
#'   length \eqn{K - 1} vector (one constant per non-baseline category) or an
#'   \eqn{n \times (K - 1)} matrix of per-category shifts on the ALR scale
#'   (columns aligned with the non-baseline categories in their original
#'   order; the baseline has no offset). Default `NULL` (no offset).
#' @param baseline Multinomial family only: the baseline category, given as a
#'   column index or a column name of `y`, or `"largest"` (default) to use the
#'   category with the largest total count over the series (ties broken by
#'   column order). A category with many observations makes the ALR
#'   transform numerically well behaved. Note that the model is \emph{not}
#'   invariant to this choice: the latent dynamics are placed on the
#'   log-ratios relative to the baseline, so if the baseline's own share
#'   moves a lot every ALR series inherits that movement. A large category
#'   with a stable share is the natural choice. For the other families a
#'   non-default value is ignored with a warning.
#' @param verbose Logical; print a progress bar. Default `FALSE`.
#' @param seed Optional random seed for reproducibility. The previous state of
#'   the global random number generator is restored after fitting.
#'
#' @return An object of class `"dynamic_fit"`: a list with three elements.
#'   \describe{
#'     \item{`draws`}{Posterior draws only. For the Poisson and binomial
#'       families a list of matrices/vectors with one row (or element) per
#'       stored draw:
#'       \describe{
#'         \item{`z`}{Latent states aligned with the observations
#'           (`draws x n`).}
#'         \item{`z0`}{The initial latent state, one period before the first
#'           observation, which carries the `init_mean` / `init_var` prior.}
#'         \item{`sig2`}{Variance of the increment leading into each
#'           \eqn{z_t} (`draws x n`); the first column is the increment from
#'           `z0`.}
#'         \item{`fitted`, `yrep`}{\strong{Unconditional} fitted means and
#'           posterior predictive replicates; under `zeros = "inflated"` they
#'           include the zero-inflation gate.}
#'         \item{`fitted_open`, `yrep_open`}{Their
#'           \strong{conditional-on-gate-open} counterparts: the latent-implied
#'           mean, and a replicate drawn straight from the observation model.}
#'         \item{`gate`, `pi_open`}{Zero-inflation gate indicators and
#'           gate-open probability; `NULL` unless `zeros = "inflated"`.}
#'         \item{`rho`, `mu`}{AR(1) coefficient (`1` under the random walk)
#'           and drift/intercept (`0` unless `include_mu = TRUE`).}
#'         \item{`innov_var`}{A representative innovation variance whose
#'           definition depends on `innovations`: the estimated constant
#'           increment variance \eqn{\sigma^2} for `"gaussian"`; the marginal
#'           increment variance for `"t"` (\eqn{\sigma^2 \nu / (\nu - 2)}) and
#'           `"mixture"` (\eqn{\sigma^2 \sum_h \eta_h \sigma_h^2}); and the
#'           average of the per-increment variances \eqn{\exp(h_t)} over the
#'           in-sample transitions for `"sv"`. For `"t"`, if a draw has
#'           \eqn{\nu \le 2} (possible only when `df_min <= 2`), the undefined
#'           variance factor is replaced by the convention `3`.}
#'         \item{`scale`}{The overall innovation scale \eqn{\sigma^2}
#'           (`NULL` for `"sv"`).}
#'         \item{`nu`}{Student-t degrees of freedom (`NULL` unless
#'           `innovations = "t"`).}
#'         \item{`mix_weight`, `mix_var`}{Mixture weights and component
#'           variances (\eqn{\sigma^2 \sigma_h^2}), `draws x mix_components`
#'           (`NULL` unless `innovations = "mixture"`). The components are
#'           exchangeable and not identified individually. Within each draw
#'           they are stored in order of increasing variance, which is a
#'           labelling convention rather than an identification, so
#'           per-component summaries are meaningful only if the components
#'           are clearly separated. Label-invariant quantities, such as
#'           `innov_var` and the forecasts, are unaffected.}
#'         \item{`sv_mu`, `sv_phi`, `sv_sigma`}{Level, persistence and
#'           volatility of the AR(1) log-variance process (`NULL` unless
#'           `innovations = "sv"`).}
#'         \item{`forecast_z`, `forecast_y`}{Latent and response forecasts
#'           when `horizon >= 1` (`NULL` otherwise); `forecast_y` is
#'           unconditional, i.e. includes the gate.}
#'       }
#'       Use `yrep`, not `yrep_open`, for posterior predictive checks; without
#'       zero inflation the conditional and unconditional pairs are identical.
#'       See [predict.dynamic_fit()].
#'       \cr\cr
#'       For the multinomial family the same components carry an extra
#'       trailing \emph{category} dimension, named by the category labels:
#'       latent quantities (`z`, `sig2`, `forecast_z`, `mix_weight`,
#'       `mix_var`) are arrays with one slice per non-baseline category, and
#'       `z0`, `rho`, `mu`, `innov_var`, `scale`, `nu` and the `sv_*`
#'       parameters are `draws x (K - 1)` matrices; response-scale quantities
#'       (`fitted` -- the expected counts \eqn{N_t p_{t,k}} --, `yrep`,
#'       `forecast_y`, and the additional `fitted_prob` / `forecast_prob`
#'       holding the category shares) are `draws x time x K` arrays including
#'       the baseline, in the original column order. `fitted_open`,
#'       `yrep_open`, `gate` and `pi_open` are `NULL`.}
#'     \item{`data`}{The observed inputs: `y`, `trials` (the row totals for the
#'       multinomial family), the series length `n`, and the resolved
#'       per-observation `offset`. For the multinomial family also `K`, the
#'       `categories` (column labels), and the `baseline` index and
#'       `baseline_name`.}
#'     \item{`spec`}{The model and MCMC specification: `family`, `innovations`,
#'       `latent_dynamics`, `include_mu`, `zeros`, `prior`, `nsave`, `nburn`,
#'       `thin`, `horizon`, and the forecast-period inputs `forecast_offset` /
#'       `forecast_trials` used for the stored forecast (`NULL` when not
#'       applicable).}
#'   }
#'
#' @examples
#' sim <- simulate_dynamic_poisson(n = 60, sigma = 0.2, log_rate0 = 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, family = "poisson", nsave = 300, nburn = 200,
#'                          seed = 1)
#' summary(fit)
#' forecast(fit, horizon = 5)
#'
#' # multinomial choice counts: three categories, baseline chosen automatically
#' simm <- simulate_dynamic_multinomial(n = 40, sigma = 0.15, trials = 100,
#'                                      alr0 = c(-1, -0.5), seed = 2)
#' fitm <- fit_dynamic_model(simm$y, family = "multinomial",
#'                           nsave = 200, nburn = 100, seed = 2)
#' fitm
#' summary(fitm)$params
#'
#' @seealso [dynamic_prior()], [forecast.dynamic_fit()], [predict.dynamic_fit()],
#'   [plot_fitted()], [structural_zero_prob()], [simulate_dynamic_multinomial()]
#' @export
fit_dynamic_model <- function(y,
                              family = c("poisson", "binomial", "multinomial"),
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
                              baseline = "largest",
                              verbose = FALSE,
                              seed = NULL) {

  family      <- match.arg(family)
  latent_dynamics <- match.arg(latent_dynamics)
  stopifnot(is.logical(include_mu), length(include_mu) == 1L, !is.na(include_mu))
  # AR(1) always carries an intercept: without it the latent process is assumed
  # to have a stationary mean of zero, which is rarely appropriate for a
  # log-rate/logit series and forces rho -> 1. Enable it silently.
  if (identical(latent_dynamics, "ar1")) include_mu <- TRUE
  innovations <- normalise_innovations(innovations)
  is_multi <- family == "multinomial"

  # resolve zero handling -----------------------------------------------------
  if (!is.logical(zero_inflation) || length(zero_inflation) != 1L ||
      is.na(zero_inflation)) {
    stop("`zero_inflation` must be TRUE or FALSE. In fit_dynamic_model() it ",
         "only switches zero inflation on (the gate-open probability is ",
         "estimated); it is not a probability as in the simulate_*() functions.",
         call. = FALSE)
  }
  if (missing(zeros)) {
    zeros <- if (zero_inflation) "inflated" else "none"
  } else {
    zeros <- match.arg(zeros)
    if (zero_inflation && zeros != "inflated") {
      stop("`zero_inflation = TRUE` conflicts with `zeros = \"", zeros, "\"`.",
           call. = FALSE)
    }
  }
  if (is_multi && zeros != "none") {
    stop("Zero inflation / missing zeros are not available for ",
         "`family = \"multinomial\"`; use `zeros = \"none\"`.", call. = FALSE)
  }
  if (!inherits(prior, "dynamic_prior")) {
    stop("`prior` must be created with dynamic_prior().", call. = FALSE)
  }

  # input validation ----------------------------------------------------------
  multi <- NULL
  if (is_multi) {
    if (is.data.frame(y)) y <- as.matrix(y)
    if (!is.matrix(y) || !is.numeric(y) || ncol(y) < 2L) {
      stop("For `family = \"multinomial\"`, `y` must be a numeric matrix with ",
           "at least 2 columns (categories).", call. = FALSE)
    }
    if (anyNA(y) || any(y < 0) || any(y != round(y))) {
      stop("`y` must contain non-negative integer counts.", call. = FALSE)
    }
    if (!is.null(trials)) {
      stop("`trials` must be NULL for `family = \"multinomial\"`: the row ",
           "totals of `y` are the trials.", call. = FALSE)
    }
    storage.mode(y) <- "double"
    n <- nrow(y); K <- ncol(y)
    if (n < 3) stop("`y` must have at least 3 observations.", call. = FALSE)
    if (is.null(colnames(y))) colnames(y) <- paste0("cat", seq_len(K))
    if (anyDuplicated(colnames(y))) stop("Column names of `y` must be unique.",
                                         call. = FALSE)
    if (all(rowSums(y) == 0)) stop("All row totals of `y` are zero.", call. = FALSE)
    baseline <- resolve_baseline(baseline, y)
    trials <- rowSums(y)
    multi <- list(K = K, categories = colnames(y), baseline = baseline,
                  baseline_name = colnames(y)[baseline])
  } else {
    if (!is.numeric(y) || anyNA(y) || any(y < 0) || any(y != round(y))) {
      stop("`y` must be a vector of non-negative integers.", call. = FALSE)
    }
    y <- as.numeric(y)
    n <- length(y)
    if (n < 3) stop("`y` must have at least 3 observations.", call. = FALSE)
    if (!identical(baseline, "largest")) {
      warning("`baseline` is only used for the multinomial family and is ignored.",
              call. = FALSE)
    }
  }

  if (family == "poisson" && !is.null(trials)) {
    warning("`trials` is ignored for the Poisson family.", call. = FALSE)
    trials <- NULL
  }
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

  # horizon: a single non-negative integer (0 = no stored forecast) ------------
  if (!is.numeric(horizon) || length(horizon) != 1L || !is.finite(horizon) ||
      horizon < 0 || horizon != round(horizon)) {
    stop("`horizon` must be a single non-negative integer.", call. = FALSE)
  }
  horizon <- as.integer(horizon)

  # validate offsets: known, fixed inputs must be finite numerics --------------
  if (!is.null(offset) &&
      (!is.numeric(offset) || anyNA(offset) || any(!is.finite(offset)))) {
    stop("`offset` must be a finite numeric vector or matrix.", call. = FALSE)
  }
  if (is_multi) {
    # multinomial offsets: scalar, a length-(K - 1) vector (one constant per
    # non-baseline category, recycled over time), or a full n x (K - 1) matrix
    offset <- expand_alr_offset(offset, n, K - 1L, "offset", "n")
  }
  has_offset <- !is.null(offset) && any(offset != 0)

  # forecast-period inputs (validated before sampling) -------------------------
  fc_inputs <- NULL
  if (horizon >= 1L) {
    fc_inputs <- resolve_forecast_inputs(family, horizon,
                                         if (is_multi) K - 1L else 1L, trials,
                                         has_offset = has_offset,
                                         forecast_offset = forecast_offset,
                                         forecast_trials = forecast_trials)
  } else if (!is.null(forecast_offset) || !is.null(forecast_trials)) {
    warning("`forecast_offset` and `forecast_trials` are ignored when ",
            "`horizon = 0`; pass them to forecast() together with a horizon.",
            call. = FALSE)
  }

  with_seed(seed, {
    out <- dynamic_sampler(
      y = y, family = family, trials = trials,
      baseline = if (is_multi) baseline else NULL,
      innovations = innovations,
      zeros = zeros, latent_dynamics = latent_dynamics, include_mu = include_mu,
      prior = prior, nsave = nsave, nburn = nburn, thin = thin,
      offset = offset, verbose = verbose
    )

    fit <- structure(
      list(
        draws = c(out$draws,
                  list(forecast_z = NULL, forecast_y = NULL, forecast_prob = NULL)),
        data  = c(list(y = y, trials = trials, n = n, offset = out$offset), multi),
        spec  = list(family = family, innovations = innovations,
                     latent_dynamics = latent_dynamics, include_mu = include_mu,
                     zeros = zeros, prior = prior, nsave = nsave, nburn = nburn,
                     thin = thin, horizon = horizon,
                     forecast_offset = NULL, forecast_trials = NULL)
      ),
      class = "dynamic_fit"
    )

    if (horizon >= 1L) {
      sims <- simulate_forecast(fit, horizon, fc_inputs$offset, fc_inputs$trials)
      fit$draws[c("forecast_z", "forecast_y", "forecast_prob")] <-
        list(sims$z, sims$y, sims$prob)
      fit$spec[c("forecast_offset", "forecast_trials")] <-
        list(if (is_multi) fc_inputs$offset else fc_inputs$offset[, 1L],
             fc_inputs$trials)
    }
    fit
  })
}

#' Expand a multinomial offset to a `rows x J` matrix (or leave it `NULL`)
#' @keywords internal
#' @noRd
expand_alr_offset <- function(x, rows, J, name, rows_name) {
  if (is.null(x) || length(x) == 1L) return(x)
  if (is.null(dim(x)) && length(x) == J) x <- matrix(x, rows, J, byrow = TRUE)
  x <- as.matrix(x)
  if (!identical(dim(x), c(rows, J))) {
    stop(sprintf("`%s` must be a scalar, a length K - 1 vector or a %s x (K - 1) matrix for the multinomial family.",
                 name, rows_name), call. = FALSE)
  }
  x
}

#' Resolve the multinomial baseline category to a column index
#' @keywords internal
#' @noRd
resolve_baseline <- function(baseline, y) {
  K <- ncol(y)
  if (is.null(baseline) || identical(baseline, "largest")) {
    return(which.max(colSums(y)))                # first maximum on ties
  }
  if (length(baseline) != 1L || anyNA(baseline)) {
    stop("`baseline` must be a single column index, column name or \"largest\".",
         call. = FALSE)
  }
  if (is.character(baseline)) {
    idx <- match(baseline, colnames(y))
    if (is.na(idx)) stop(sprintf("`baseline` \"%s\" is not a column of `y`.", baseline),
                         call. = FALSE)
    return(idx)
  }
  if (!is.numeric(baseline) || baseline != round(baseline) ||
      baseline < 1 || baseline > K) {
    stop("`baseline` must be a column index in 1..K, a column name or \"largest\".",
         call. = FALSE)
  }
  as.integer(baseline)
}

#' @export
print.dynamic_fit <- function(x, ...) {
  sp <- x$spec
  is_multi <- identical(sp$family, "multinomial")
  cat("<dynamic_fit>\n")
  cat(sprintf("  family      : %s\n", sp$family))
  if (is_multi) {
    cat(sprintf("  categories  : %d (%s); baseline = %s\n", x$data$K,
                paste(x$data$categories, collapse = ", "), x$data$baseline_name))
  }
  fmt3 <- function(v) paste(sprintf("%.3f", v), collapse = ", ")
  cat(sprintf("  dynamics    : %s%s\n", sp$latent_dynamics,
              if (identical(sp$latent_dynamics, "ar1"))
                sprintf("  (rho post. mean = %s)", fmt3(post_mean_cols(x$draws$rho)))
              else "  (rho = 1)"))
  if (isTRUE(sp$include_mu)) {
    cat(sprintf("  %-12s: post. mean = %s\n",
                if (identical(sp$latent_dynamics, "ar1")) "intercept mu" else "drift mu",
                fmt3(post_mean_cols(x$draws$mu))))
  }
  cat(sprintf("  innovations : %s\n", sp$innovations))
  cat(sprintf("  zeros       : %s\n", sp$zeros))
  if (is_multi) {
    cat(sprintf("  observations: %d  (zero-total rows: %d)\n",
                x$data$n, sum(x$data$trials == 0)))
  } else {
    cat(sprintf("  observations: %d  (zeros: %d)\n",
                x$data$n, sum(x$data$y == 0)))
  }
  cat(sprintf("  draws kept  : %d\n", sp$nsave))
  if (sp$horizon >= 1L) {
    cat(sprintf("  forecast    : horizon %d (stored)\n", sp$horizon))
  }
  invisible(x)
}

#' Posterior means of a draws vector or of each column of a draws matrix
#' @keywords internal
#' @noRd
post_mean_cols <- function(d) {
  if (is.null(dim(d))) mean(d) else colMeans(d)
}
