# Summary ---------------------------------------------------------------------

#' Summarise a fitted dynamic model
#'
#' @param object A `"dynamic_fit"` object.
#' @param probs Quantile probabilities. Default `c(0.025, 0.5, 0.975)`.
#' @param ... Unused.
#'
#' @return An object of class `"summary.dynamic_fit"` containing posterior
#'   summaries of the global parameters (the innovation standard deviation
#'   `innov_sd` -- the estimated SD for `"gaussian"` innovations, the marginal
#'   SD for `"t"` and `"mixture"`, and the series-average SD for `"sv"` --
#'   and where relevant the t degrees of freedom, the AR(1) coefficient, the
#'   drift/intercept, and the zero-inflation gate-open probability) and the
#'   fitted-value summary.
#' @examples
#' sim <- simulate_dynamic_poisson(50, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 200, nburn = 100)
#' summary(fit)
#' @export
summary.dynamic_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  d <- object$draws
  sp <- object$spec

  # innovation standard deviation, per draw: the estimated SD for "gaussian",
  # the marginal SD implied by the scale and mixing distribution for "t" and
  # "mixture", and the root of the series-average variance for "sv".
  inc_sd <- sqrt(d$innov_var)
  params <- summarise_draws(inc_sd, probs)
  rownames(params) <- "innov_sd"

  if (sp$innovations == "t") {
    params <- rbind(params, summarise_draws(d$nu, probs))
    rownames(params)[nrow(params)] <- "t_df"
  }
  if (identical(sp$latent_dynamics, "ar1")) {
    params <- rbind(params, summarise_draws(d$rho, probs))
    rownames(params)[nrow(params)] <- "ar1_rho"
  }
  if (isTRUE(sp$include_mu)) {
    params <- rbind(params, summarise_draws(d$mu, probs))
    rownames(params)[nrow(params)] <-
      if (identical(sp$latent_dynamics, "ar1")) "intercept_mu" else "drift_mu"
  }
  if (sp$zeros == "inflated") {
    params <- rbind(params, summarise_draws(d$pi_open, probs))
    rownames(params)[nrow(params)] <- "gate_open_prob"
  }

  fitted_summary <- summarise_draws(d$fitted, probs)

  out <- list(
    spec = sp,
    n = object$data$n,
    n_zero = sum(object$data$y == 0),
    params = params,
    fitted = fitted_summary
  )
  class(out) <- "summary.dynamic_fit"
  out
}

#' @export
print.summary.dynamic_fit <- function(x, ...) {
  cat("Dynamic count model summary\n")
  cat(sprintf("  family = %s | dynamics = %s | innovations = %s | zeros = %s\n",
              x$spec$family, x$spec$latent_dynamics %||% "rw",
              x$spec$innovations, x$spec$zeros))
  cat(sprintf("  %d observations (%d zeros), %d posterior draws\n\n",
              x$n, x$n_zero, x$spec$nsave))
  cat("Global parameters (posterior summaries):\n")
  print(round(x$params, 4))
  cat(sprintf("\nFitted values: range of posterior means [%.2f, %.2f]\n",
              min(x$fitted$mean), max(x$fitted$mean)))
  invisible(x)
}
