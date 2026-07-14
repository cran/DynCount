# Simulation ------------------------------------------------------------------

#' Simulate a Poisson dynamic series
#'
#' Generates a latent log-rate process (random walk or AR(1)) and Poisson
#' counts, optionally with zero inflation.
#'
#' @param n Number of observations.
#' @param sigma Standard deviation of the latent increments (Gaussian).
#' @param log_rate0 Initial log-rate \eqn{z_1}. Default `1`.
#' @param zero_inflation Probability that the gate is **closed** (i.e. the
#'   probability of a structural zero) at each time point. `0` (default) gives
#'   an ordinary Poisson series.
#' @param rho AR(1) coefficient of the latent process
#'   \eqn{z_t = \mu + \rho z_{t-1} + \varepsilon_t}. Default `1` (a random walk).
#' @param mu Drift (random walk) / intercept (AR(1)) of the latent process.
#'   Default `0`.
#' @param offset Known log-exposure offset (length 1 or `n`); the Poisson mean
#'   is \eqn{\exp(\mathrm{offset}_t + z_t)}. Default `0`.
#' @param seed Optional random seed.
#'
#' @return A list with components `y` (observed counts), `log_rate` (the latent
#'   log-rate path \eqn{z_t}), `rate` (the mean `exp(offset + log_rate)`),
#'   `offset`, and `structural` (logical, `TRUE` where a structural zero was
#'   forced).
#'
#' @examples
#' sim <- simulate_dynamic_poisson(n = 50, sigma = 0.2, log_rate0 = 2,
#'                            zero_inflation = 0.2, seed = 1)
#' table(sim$y == 0, sim$structural)
#' @export
simulate_dynamic_poisson <- function(n, sigma, log_rate0 = 1,
                                zero_inflation = 0, rho = 1, mu = 0,
                                offset = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- check_count(n, "n")
  stopifnot(sigma >= 0, zero_inflation >= 0, zero_inflation < 1)
  if (length(offset) == 1L) offset <- rep(offset, n)
  if (length(offset) != n) stop("`offset` must be length 1 or n.", call. = FALSE)
  z <- numeric(n); z[1] <- log_rate0
  eps <- stats::rnorm(n - 1, 0, sigma)
  for (t in seq_len(n - 1)) z[t + 1] <- mu + rho * z[t] + eps[t]
  rate <- exp(offset + z)                     # mean incl. offset
  y <- stats::rpois(n, rate)
  structural <- stats::runif(n) < zero_inflation
  y[structural] <- 0
  list(y = y, log_rate = z, rate = rate, offset = offset,
       structural = structural)
}

#' Simulate a binomial dynamic series
#'
#' Generates a latent logit process (random walk or AR(1)) and binomial counts, optionally with
#' structural (zero-inflation) zeros.
#'
#' @param n Number of observations.
#' @param sigma Standard deviation of the latent increments (Gaussian).
#' @param trials Number of trials: a single number (recycled) or a length-`n`
#'   vector.
#' @param logit0 Initial logit \eqn{z_1}. Default `0`.
#' @param zero_inflation Structural-zero probability. With probability
#'   `zero_inflation` an observation is forced to a structural zero (gate
#'   closed) regardless of the binomial draw. Default `0` (no inflation).
#' @param rho AR(1) coefficient of the latent process
#'   \eqn{z_t = \mu + \rho z_{t-1} + \varepsilon_t}. Default `1` (a random walk).
#' @param mu Drift (random walk) / intercept (AR(1)) of the latent process.
#'   Default `0`.
#' @param offset Known offset on the logit scale (length 1 or `n`); the success
#'   probability is \eqn{\mathrm{logit}^{-1}(\mathrm{offset}_t + z_t)}.
#'   Default `0`.
#' @param seed Optional random seed.
#'
#' @return A list with components `y` (successes), `trials`, `logit` (latent
#'   path \eqn{z_t}), `prob` (`plogis(offset + logit)`), `offset`, and
#'   `structural` (logical, `TRUE` for structural zeros).
#'
#' @examples
#' sim <- simulate_dynamic_binomial(n = 50, sigma = 0.15, trials = 40, seed = 1)
#' head(sim$y)
#' # with structural zeros:
#' zi <- simulate_dynamic_binomial(50, 0.15, trials = 40, zero_inflation = 0.2, seed = 1)
#' mean(zi$structural)
#' @export
simulate_dynamic_binomial <- function(n, sigma, trials, logit0 = 0,
                                 zero_inflation = 0, rho = 1, mu = 0,
                                 offset = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- check_count(n, "n")
  stopifnot(sigma >= 0, zero_inflation >= 0, zero_inflation < 1)
  if (length(trials) == 1L) trials <- rep(trials, n)
  if (length(trials) != n) stop("`trials` must be length 1 or n.", call. = FALSE)
  if (length(offset) == 1L) offset <- rep(offset, n)
  if (length(offset) != n) stop("`offset` must be length 1 or n.", call. = FALSE)
  z <- numeric(n); z[1] <- logit0
  eps <- stats::rnorm(n - 1, 0, sigma)
  for (t in seq_len(n - 1)) z[t + 1] <- mu + rho * z[t] + eps[t]
  prob <- stats::plogis(offset + z)
  y <- stats::rbinom(n, size = trials, prob = prob)
  structural <- stats::runif(n) < zero_inflation
  y[structural] <- 0
  list(y = y, trials = trials, logit = z, prob = prob, offset = offset,
       structural = structural)
}
