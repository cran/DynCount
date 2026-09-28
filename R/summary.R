# Summary ---------------------------------------------------------------------

#' Summarise a fitted dynamic model
#'
#' @param object A `"dynamic_fit"` object.
#' @param probs Quantile probabilities. Default `c(0.025, 0.5, 0.975)`.
#' @param ... Unused.
#'
#' @return An object of class `"summary.dynamic_fit"` containing posterior
#'   summaries of the global parameters (`params`, one row per parameter) and
#'   of the fitted values (`fitted`). The rows of `params` and the draws they
#'   summarise are
#'   \describe{
#'     \item{`innov_sd`}{`sqrt(draws$innov_var)`: the estimated increment SD
#'       for `"gaussian"` innovations, the marginal SD for `"t"` and
#'       `"mixture"`, and the root of the series-average variance for `"sv"`.}
#'     \item{`t_df`}{`draws$nu`, the Student-t degrees of freedom.}
#'     \item{`sv_mu`, `sv_phi`, `sv_sigma`}{`draws$sv_mu`, `draws$sv_phi` and
#'       `draws$sv_sigma`, the parameters of the AR(1) log-variance process.}
#'     \item{`ar1_rho`}{`draws$rho`, the AR(1) coefficient.}
#'     \item{`drift_mu` / `intercept_mu`}{`draws$mu`, the random-walk drift or
#'       the AR(1) intercept.}
#'     \item{`gate_open_prob`}{`draws$pi_open`, the probability that the
#'       zero-inflation gate is open (one minus the structural-zero
#'       probability).}
#'   }
#'   Rows are included only where relevant. For `"mixture"` innovations only
#'   `innov_sd` is reported. The mixture components are exchangeable and not
#'   identified individually, so per-component summaries are not given. The
#'   draws of the weights and component variances are in `draws$mix_weight`
#'   and `draws$mix_var`. For the multinomial family every
#'   latent parameter is reported once per non-baseline category, with the
#'   category label in square brackets (e.g. `innov_sd[B]`), and the fitted
#'   summary is a long data frame with a `category` column covering all `K`
#'   categories.
#' @examples
#' sim <- simulate_dynamic_poisson(50, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 200, nburn = 100, seed = 1)
#' summary(fit)
#' @export
summary.dynamic_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  d <- object$draws
  sp <- object$spec
  is_multi <- identical(sp$family, "multinomial")

  # one row per latent series (1 for poisson/binomial, K - 1 for multinomial)
  lat_names <- if (is_multi) object$data$categories[-object$data$baseline] else NULL
  as_mat <- function(v) if (is.null(dim(v))) matrix(v, ncol = 1L) else v
  label <- function(base) if (is_multi) sprintf("%s[%s]", base, lat_names) else base

  # innovation standard deviation, per draw: the estimated SD for "gaussian",
  # the marginal SD implied by the scale and mixing distribution for "t" and
  # "mixture", and the root of the series-average variance for "sv".
  params <- summarise_draws(sqrt(as_mat(d$innov_var)), probs)
  rownames(params) <- label("innov_sd")

  add <- function(params, draws, base) {
    s <- summarise_draws(as_mat(draws), probs)
    rownames(s) <- label(base)
    rbind(params, s)
  }
  if (sp$innovations == "t") params <- add(params, d$nu, "t_df")
  if (sp$innovations == "sv" && !is.null(d$sv_mu)) {
    params <- add(params, d$sv_mu, "sv_mu")
    params <- add(params, d$sv_phi, "sv_phi")
    params <- add(params, d$sv_sigma, "sv_sigma")
  }
  if (identical(sp$latent_dynamics, "ar1")) params <- add(params, d$rho, "ar1_rho")
  if (isTRUE(sp$include_mu)) {
    params <- add(params, d$mu,
                  if (identical(sp$latent_dynamics, "ar1")) "intercept_mu" else "drift_mu")
  }
  if (sp$zeros == "inflated") {
    params <- rbind(params, summarise_draws(d$pi_open, probs))
    rownames(params)[nrow(params)] <- "gate_open_prob"
  }

  fitted_summary <- if (is_multi) {
    summarise_draws_long(d$fitted, probs, "time")
  } else {
    summarise_draws(d$fitted, probs)
  }

  out <- list(
    spec = sp,
    n = object$data$n,
    n_zero = if (is_multi) sum(object$data$trials == 0) else sum(object$data$y == 0),
    K = object$data$K,
    baseline_name = object$data$baseline_name,
    params = params,
    fitted = fitted_summary
  )
  class(out) <- "summary.dynamic_fit"
  out
}

#' @export
print.summary.dynamic_fit <- function(x, ...) {
  is_multi <- identical(x$spec$family, "multinomial")
  cat("Dynamic count model summary\n")
  cat(sprintf("  family = %s | dynamics = %s | innovations = %s | zeros = %s\n",
              x$spec$family, x$spec$latent_dynamics %||% "rw",
              x$spec$innovations, x$spec$zeros))
  if (is_multi) {
    cat(sprintf("  %d categories, baseline = %s\n", x$K, x$baseline_name))
    cat(sprintf("  %d observations (%d zero-total rows), %d posterior draws\n\n",
                x$n, x$n_zero, x$spec$nsave))
  } else {
    cat(sprintf("  %d observations (%d zeros), %d posterior draws\n\n",
                x$n, x$n_zero, x$spec$nsave))
  }
  cat("Global parameters (posterior summaries):\n")
  print(round(x$params, 4))
  if (is_multi) {
    rng <- tapply(x$fitted$mean, x$fitted$category, range)
    cat("\nFitted values: range of posterior means by category\n")
    for (nm in names(rng)) cat(sprintf("  %-12s [%.2f, %.2f]\n", nm, rng[[nm]][1], rng[[nm]][2]))
  } else {
    cat(sprintf("\nFitted values: range of posterior means [%.2f, %.2f]\n",
                min(x$fitted$mean), max(x$fitted$mean)))
  }
  invisible(x)
}
