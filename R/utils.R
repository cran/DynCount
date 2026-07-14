# Internal helpers ------------------------------------------------------------

# Null-coalescing helper (base R gained `%||%` in 4.4; kept here for R (>= 3.5)).
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Generalised first-difference (GMRF structure) matrix
#'
#' Returns the \eqn{(N-1) \times N} operator \eqn{D_\rho} whose rows encode the
#' latent state evolution \eqn{z_t - \rho\, z_{t-1}}. The GMRF precision is then
#' \eqn{P = D_\rho^\top \mathrm{diag}(1/\sigma^2) D_\rho}. With \eqn{\rho = 1}
#' this is the ordinary first-difference matrix `diff(diag(N))` (the random-walk
#' case); with \eqn{\rho \ne 1} it gives the AR(1) GMRF.
#'
#' @param N integer; length of the (padded) latent state vector.
#' @param rho numeric; AR(1) coefficient (`1` for the random walk).
#' @return an `(N-1) x N` matrix.
#' @keywords internal
#' @noRd
make_gmrf_diff <- function(N, rho = 1) {
  D <- matrix(0, N - 1L, N)
  idx <- seq_len(N - 1L)
  D[cbind(idx, idx)]      <- -rho
  D[cbind(idx, idx + 1L)] <- 1
  D
}

#' Tridiagonal bands of the GMRF precision and the drift linear term
#'
#' For the AR(1)/RW GMRF precision \eqn{P = D_\rho^\top \mathrm{diag}(w) D_\rho}
#' (with \eqn{D_\rho} the generalised first-difference operator) plus a proper
#' \eqn{N(m_0, v_0)} prior on the first state, this returns the diagonal and the
#' single off-diagonal of the (symmetric tridiagonal) \eqn{P}, together with the
#' linear term \eqn{b = \mu\, D_\rho^\top w} (with the initial-state prior term
#' added to \eqn{b_1}). Computed directly in \eqn{O(N)} without forming the
#' dense matrix. `off[i]` is \eqn{P_{i,i+1} = P_{i+1,i}}.
#'
#' @param w numeric vector of increment precisions `1 / sig2` (length `N - 1`).
#' @param rho AR(1) coefficient (`1` for the random walk).
#' @param mu drift/intercept (`0` if disabled).
#' @param v0,m0 variance and mean of the proper prior on the first state.
#' @param N integer length of the (padded) latent state vector.
#' @return a list with `diag` (length `N`), `off` (length `N - 1`) and `b`
#'   (length `N`).
#' @keywords internal
#' @noRd
gmrf_bands <- function(w, rho, mu, v0, m0, N) {
  d <- numeric(N)
  d[1] <- rho^2 * w[1]
  if (N >= 3L) d[2:(N - 1L)] <- rho^2 * w[2:(N - 1L)] + w[1:(N - 2L)]
  d[N] <- w[N - 1L]
  d[1] <- d[1] + 1 / v0                          # proper initial-state prior

  off <- -rho * w                                 # P[i, i+1], length N - 1

  b <- numeric(N)
  if (mu != 0) {                                  # b = mu * (D_rho' w)
    b[1] <- -rho * w[1]
    if (N >= 3L) b[2:(N - 1L)] <- -rho * w[2:(N - 1L)] + w[1:(N - 2L)]
    b[N] <- w[N - 1L]
    b <- mu * b
  }
  b[1] <- b[1] + m0 / v0
  list(diag = d, off = off, b = b)
}

#' Numerically stable log-sum-exp
#' @param x numeric vector
#' @return `log(sum(exp(x)))` computed stably
#' @keywords internal
#' @noRd
log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

#' Numerically stable log of the mean of exp(x)
#' @param x numeric vector
#' @return `log(mean(exp(x)))` computed stably
#' @keywords internal
#' @noRd
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

#' Stable log(exp(a) + exp(b))
#' @keywords internal
#' @noRd
logspace_add <- function(a, b) {
  m <- pmax(a, b)
  ifelse(is.finite(m), m + log(exp(a - m) + exp(b - m)), m)
}

#' Draw from a Gumbel(0, 1) distribution (base-R replacement for
#' `extraDistr::rgumbel`). Used for the Gumbel-max categorical sampler.
#' @keywords internal
#' @noRd
rgumbel0 <- function(n) {
  -log(-log(stats::runif(n)))
}

#' Draw a single Dirichlet vector (base-R replacement for
#' `MCMCpack::rdirichlet`).
#' @param alpha numeric vector of concentration parameters
#' @keywords internal
#' @noRd
rdirichlet1 <- function(alpha) {
  g <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
  g / sum(g)
}

#' Posterior summary (mean, sd, quantiles) of a draws matrix, column-wise
#' @param draws matrix with one MCMC draw per row
#' @param probs quantile probabilities
#' @keywords internal
#' @noRd
summarise_draws <- function(draws, probs = c(0.025, 0.5, 0.975)) {
  if (is.null(dim(draws))) draws <- matrix(draws, ncol = 1L)
  qs <- t(apply(draws, 2L, stats::quantile, probs = probs, names = FALSE))
  out <- data.frame(
    mean = colMeans(draws),
    sd   = apply(draws, 2L, stats::sd)
  )
  for (j in seq_along(probs)) {
    out[[paste0("q", probs[j] * 100)]] <- qs[, j]
  }
  out
}

#' Extract the (lower, median, upper) quantile columns of a [summarise_draws()]
#' table by name, given the half-alpha `a` used to build it (probs
#' `c(a, 0.5, 1 - a)`). Indexing by name avoids depending on column position.
#' @keywords internal
#' @noRd
qband <- function(s, a) {
  list(lo = s[[paste0("q", a * 100)]],
       md = s[["q50"]],
       hi = s[[paste0("q", (1 - a) * 100)]])
}

#' Light progress reporter that does not require the `progress` package.
#' @keywords internal
#' @noRd
make_progress <- function(total, verbose) {
  if (!isTRUE(verbose)) {
    return(function(i) invisible(NULL))
  }
  width <- 40L
  last <- -1L
  function(i) {
    pct <- floor(100 * i / total)
    if (pct != last) {
      done <- floor(width * i / total)
      bar  <- paste0(strrep("=", done), strrep(" ", width - done))
      message(sprintf("\r sampling... [%s] %3d%%", bar, pct),
              appendLF = (i == total))
      last <<- pct
    }
    invisible(NULL)
  }
}

#' Validate that a value is a single positive integer-ish number
#' @keywords internal
#' @noRd
check_count <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x <= 0 || x != round(x)) {
    stop(sprintf("`%s` must be a single positive integer.", name), call. = FALSE)
  }
  as.integer(x)
}

#' Validate that a value is a single non-negative integer-ish number
#'
#' Like [check_count()] but permits `0` (used for `nburn`, where a zero
#' burn-in is a legitimate request).
#' @keywords internal
#' @noRd
check_count0 <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x < 0 || x != round(x)) {
    stop(sprintf("`%s` must be a single non-negative integer.", name),
         call. = FALSE)
  }
  as.integer(x)
}
