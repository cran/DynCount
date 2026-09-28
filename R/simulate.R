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
#' @param seed Optional random seed. The previous state of the global random
#'   number generator is restored afterwards.
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
  n <- check_count(n, "n")
  stopifnot(sigma >= 0, zero_inflation >= 0, zero_inflation < 1)
  if (length(offset) == 1L) offset <- rep(offset, n)
  if (length(offset) != n) stop("`offset` must be length 1 or n.", call. = FALSE)
  with_seed(seed, {
    z <- numeric(n); z[1] <- log_rate0
    eps <- stats::rnorm(n - 1, 0, sigma)
    for (t in seq_len(n - 1)) z[t + 1] <- mu + rho * z[t] + eps[t]
    rate <- exp(offset + z)                     # mean incl. offset
    y <- stats::rpois(n, rate)
    structural <- stats::runif(n) < zero_inflation
    y[structural] <- 0
    list(y = y, log_rate = z, rate = rate, offset = offset,
         structural = structural)
  })
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
#' @param seed Optional random seed. The previous state of the global random
#'   number generator is restored afterwards.
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
  n <- check_count(n, "n")
  stopifnot(sigma >= 0, zero_inflation >= 0, zero_inflation < 1)
  if (length(trials) == 1L) trials <- rep(trials, n)
  if (length(trials) != n) stop("`trials` must be length 1 or n.", call. = FALSE)
  if (!is.numeric(trials) || anyNA(trials) || any(trials <= 0) ||
      any(trials != round(trials))) {
    stop("`trials` must be positive integers.", call. = FALSE)
  }
  if (length(offset) == 1L) offset <- rep(offset, n)
  if (length(offset) != n) stop("`offset` must be length 1 or n.", call. = FALSE)
  with_seed(seed, {
    z <- numeric(n); z[1] <- logit0
    eps <- stats::rnorm(n - 1, 0, sigma)
    for (t in seq_len(n - 1)) z[t + 1] <- mu + rho * z[t] + eps[t]
    prob <- stats::plogis(offset + z)
    y <- stats::rbinom(n, size = trials, prob = prob)
    structural <- stats::runif(n) < zero_inflation
    y[structural] <- 0
    list(y = y, trials = trials, logit = z, prob = prob, offset = offset,
         structural = structural)
  })
}

#' Simulate a multinomial dynamic series
#'
#' Generates \eqn{K - 1} independent latent additive-log-ratio (ALR) processes
#' (random walk or AR(1)), one per non-baseline category, and multinomial
#' choice counts with known totals. The ALR series
#' \eqn{z_{t,k} = \log(p_{t,k} / p_{t,b})} share no parameters: `sigma`,
#' `rho` and `mu` may each be a single value (recycled) or a vector of length
#' \eqn{K - 1} giving one value per non-baseline category.
#'
#' @param n Number of observations.
#' @param sigma Standard deviation(s) of the latent increments (Gaussian);
#'   length 1 or `K - 1`.
#' @param trials Total count per period: a single number (recycled) or a
#'   length-`n` vector.
#' @param alr0 Numeric vector of length \eqn{K - 1}: the initial ALR values
#'   \eqn{z_{1,k}} of the non-baseline categories, in column order. Its
#'   length determines the number of categories \eqn{K}. Default `c(0, 0)`
#'   (three equally likely categories).
#' @param baseline The baseline category: a column index in `1..K` of the
#'   returned count matrix, or one of the `categories`. Default `K` (last
#'   column). Pass the returned `baseline` to [fit_dynamic_model()] to fit the
#'   model on the same ALR scale as the simulation.
#' @param rho AR(1) coefficient(s) of the latent processes
#'   \eqn{z_{t,k} = \mu_k + \rho_k z_{t-1,k} + \varepsilon_{t,k}}; length 1 or
#'   `K - 1`. Default `1` (random walks).
#' @param mu Drift (random walk) / intercept (AR(1)) of the latent processes;
#'   length 1 or `K - 1`. Default `0`.
#' @param offset Known offset on the ALR scale: a scalar, a length
#'   \eqn{K - 1} vector (one constant per non-baseline category) or an
#'   \eqn{n \times (K - 1)} matrix (columns aligned with the non-baseline
#'   categories). Default `0`.
#' @param categories Optional character vector of length \eqn{K} with the
#'   category labels (column names of the returned counts). Default
#'   `"cat1"`, ..., `"catK"`.
#' @param seed Optional random seed. The previous state of the global random
#'   number generator is restored afterwards.
#'
#' @return A list with components `y` (an `n x K` matrix of counts with the
#'   baseline in column `baseline`), `trials`, `alr` (the `n x (K - 1)` latent
#'   ALR paths \eqn{z_{t,k}}), `prob` (the `n x K` matrix of category
#'   probabilities), `offset`, `baseline` (the column index) and `categories`.
#'
#' @examples
#' sim <- simulate_dynamic_multinomial(n = 50, sigma = 0.15, trials = 200,
#'                                     alr0 = c(-0.5, 0.5), seed = 1)
#' head(sim$y)
#' colSums(sim$y)
#' # fit on the simulated ALR scale by passing the simulated baseline
#' fit <- fit_dynamic_model(sim$y, family = "multinomial", baseline = sim$baseline,
#'                          nsave = 200, nburn = 100, seed = 1)
#' summary(fit)$params
#' @export
simulate_dynamic_multinomial <- function(n, sigma, trials, alr0 = c(0, 0),
                                         baseline = length(alr0) + 1L,
                                         rho = 1, mu = 0, offset = 0,
                                         categories = NULL, seed = NULL) {
  n <- check_count(n, "n")
  J <- length(alr0)
  if (J < 1L || anyNA(alr0)) stop("`alr0` must be a numeric vector of length K - 1 >= 1.",
                                  call. = FALSE)
  K <- J + 1L
  rec <- function(x, name) {
    if (length(x) == 1L) x <- rep(x, J)
    if (length(x) != J) stop(sprintf("`%s` must have length 1 or K - 1.", name),
                             call. = FALSE)
    x
  }
  sigma <- rec(sigma, "sigma"); rho <- rec(rho, "rho"); mu <- rec(mu, "mu")
  stopifnot(all(sigma >= 0))
  if (is.null(categories)) categories <- paste0("cat", seq_len(K))
  if (length(categories) != K) stop("`categories` must have length K.", call. = FALSE)
  if (is.character(baseline) && length(baseline) == 1L) {
    b <- match(baseline, categories)
    if (is.na(b)) stop(sprintf("`baseline` \"%s\" is not one of `categories`.", baseline),
                       call. = FALSE)
    baseline <- b
  }
  if (!is.numeric(baseline) || length(baseline) != 1L || baseline < 1 ||
      baseline > K || baseline != round(baseline)) {
    stop("`baseline` must be a column index in 1..K or one of `categories`.",
         call. = FALSE)
  }
  baseline <- as.integer(baseline)
  if (length(trials) == 1L) trials <- rep(trials, n)
  if (length(trials) != n) stop("`trials` must be length 1 or n.", call. = FALSE)
  if (length(offset) == 1L) {
    offset <- matrix(offset, n, J)
  } else if (is.null(dim(offset)) && length(offset) == J) {
    offset <- matrix(offset, n, J, byrow = TRUE)
  }
  offset <- as.matrix(offset)
  if (!identical(dim(offset), c(n, J)))
    stop("`offset` must be a scalar, a length K - 1 vector or an n x (K - 1) matrix.",
         call. = FALSE)

  with_seed(seed, {
    z <- matrix(NA_real_, n, J)
    z[1L, ] <- alr0
    for (k in seq_len(J)) {
      eps <- stats::rnorm(n - 1, 0, sigma[k])
      for (t in seq_len(n - 1)) z[t + 1L, k] <- mu[k] + rho[k] * z[t, k] + eps[t]
    }
    prob <- alr_to_prob(offset + z, baseline)
    y <- t(rmultinom_rows(trials, prob))
    storage.mode(y) <- "double"
    colnames(y) <- colnames(prob) <- categories
    colnames(z) <- categories[-baseline]
    list(y = y, trials = trials, alr = z, prob = prob, offset = offset,
         baseline = baseline, categories = categories)
  })
}
