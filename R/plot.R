# Plotting --------------------------------------------------------------------
#
# All plots use base graphics to keep dependencies minimal. Each function draws
# a posterior median with a shaded credible band where appropriate.

#' @keywords internal
#' @noRd
band <- function(x, lo, hi, col) {
  graphics::polygon(c(x, rev(x)), c(lo, rev(hi)), col = col, border = NA)
}

#' Plot the fitted latent trajectory
#'
#' Shows the posterior median of the latent state (log-rate for the Poisson
#' family, logit for the binomial family) with a credible band.
#'
#' @param object A `"dynamic_fit"` object.
#' @param level Credible level for the band. Default `0.95`.
#' @param ... Passed to [graphics::plot()].
#' @return Invisibly, the summary data frame used for plotting.
#' @examples
#' sim <- simulate_dynamic_poisson(60, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 200)
#' plot_latent(fit)
#' @export
plot_latent <- function(object, level = 0.95, ...) {
  a <- (1 - level) / 2
  n <- object$data$n
  z <- object$draws$z[, 2:(n + 1L), drop = FALSE]        # observation states only
  s <- summarise_draws(z, c(a, 0.5, 1 - a))
  q <- qband(s, a)
  tt <- seq_len(nrow(s))
  ylab <- if (object$spec$family == "poisson") "latent log-rate" else "latent logit"
  graphics::plot(tt, q$md, type = "n",
                 ylim = range(q$lo, q$hi), xlab = "time", ylab = ylab,
                 main = "Latent trajectory", ...)
  band(tt, q$lo, q$hi, grDevices::adjustcolor("steelblue", 0.3))
  graphics::lines(tt, q$md, col = "steelblue", lwd = 2)
  invisible(s)
}

#' Plot observed versus fitted values
#'
#' Observed series with the posterior median fitted mean and a credible band on
#' the response scale.
#'
#' @inheritParams plot_latent
#' @param level Credible level. Default `0.95`.
#' @return Invisibly, the fitted summary.
#' @examples
#' sim <- simulate_dynamic_poisson(60, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 200)
#' plot_fitted(fit)
#' @export
plot_fitted <- function(object, level = 0.95, ...) {
  a <- (1 - level) / 2
  s <- summarise_draws(object$draws$fitted, c(a, 0.5, 1 - a))
  q <- qband(s, a)
  y <- object$data$y
  tt <- seq_along(y)
  graphics::plot(tt, y, type = "n", ylim = range(0, y, q$hi),
                 xlab = "time", ylab = "count",
                 main = "Observed vs fitted", ...)
  band(tt, q$lo, q$hi, grDevices::adjustcolor("tomato", 0.25))
  graphics::lines(tt, q$md, col = "tomato", lwd = 2)
  graphics::points(tt, y, pch = 16, cex = 0.6, col = "grey20")
  graphics::legend("topleft", c("observed", "fitted (median)", "band"),
                   pch = c(16, NA, 15), lty = c(NA, 1, NA),
                   col = c("grey20", "tomato", grDevices::adjustcolor("tomato", 0.25)),
                   bty = "n", cex = 0.8)
  invisible(s)
}

#' Plot observed history, fitted values and forecast with uncertainty
#'
#' Shows the observed series together with the in-sample fitted mean (median and
#' credible band) and the forecast (median and credible band) over the forecast
#' horizon. The forecast draws are those produced during sampling, so the model
#' must have been fitted with `horizon >= 1` (see [forecast.dynamic_fit()]).
#'
#' @param object A `"dynamic_fit"` object fitted with `horizon >= 1`.
#' @param level Credible level. Default `0.95`.
#' @param ... Passed to [graphics::plot()].
#' @return Invisibly, the `"dynamic_forecast"` object.
#' @examples
#' sim <- simulate_dynamic_poisson(60, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 200, horizon = 12)
#' plot_forecast(fit)
#' @export
plot_forecast <- function(object, level = 0.95, ...) {
  a  <- (1 - level) / 2
  pr <- c(a, 0.5, 1 - a)
  fc <- forecast(object, probs = pr)
  fit_s <- summarise_draws(object$draws$fitted, pr)   # in-sample fitted band
  qf <- qband(fit_s, a)                               # fitted lo/md/hi
  qc <- qband(fc$summary, a)                          # forecast lo/md/hi
  y <- object$data$y
  n <- length(y)
  h <- fc$horizon
  fh <- (n + 1):(n + h)
  tt <- seq_len(n)

  ylim <- range(0, y, qf$lo, qf$hi, qc$lo, qc$hi)
  graphics::plot(tt, y, type = "n", xlim = c(1, n + h), ylim = ylim,
                 xlab = "time", ylab = "count",
                 main = "Observed, fitted and forecast", ...)
  # in-sample fitted band + median
  band(tt, qf$lo, qf$hi, grDevices::adjustcolor("tomato", 0.20))
  graphics::lines(tt, qf$md, col = "tomato", lwd = 2)
  # observed
  graphics::points(tt, y, pch = 16, cex = 0.6, col = "grey20")
  # forecast band + median
  band(fh, qc$lo, qc$hi, grDevices::adjustcolor("purple", 0.20))
  graphics::lines(fh, qc$md, col = "purple", lwd = 2)
  graphics::abline(v = n + 0.5, lty = 3, col = "grey60")
  graphics::legend("topleft",
                   c("observed", "fitted (median)", "forecast (median)"),
                   pch = c(16, NA, NA), lty = c(NA, 1, 1),
                   col = c("grey20", "tomato", "purple"), bty = "n", cex = 0.8)
  invisible(fc)
}

#' Plot zero-inflation diagnostics
#'
#' Bar chart of the posterior probability that each observed zero is a
#' structural zero. Requires a model fitted with `zeros = "inflated"`.
#'
#' @param object A `"dynamic_fit"` object.
#' @param ... Passed to [graphics::barplot()].
#' @return Invisibly, the [structural_zero_prob()] data frame.
#' @examples
#' sim <- simulate_dynamic_poisson(80, 0.2, 2, zero_inflation = 0.3, seed = 1)
#' fit <- fit_dynamic_model(sim$y, zero_inflation = TRUE, nsave = 300, nburn = 200)
#' plot_zero_inflation(fit)
#' @export
plot_zero_inflation <- function(object, ...) {
  szp <- structural_zero_prob(object, zeros_only = TRUE)
  if (nrow(szp) == 0) {
    message("No observed zeros to display.")
    return(invisible(szp))
  }
  graphics::barplot(szp$p_structural, names.arg = szp$time, las = 2,
                    col = grDevices::adjustcolor("darkorange", 0.7),
                    border = NA, ylim = c(0, 1),
                    ylab = "P(structural zero)", xlab = "time index",
                    main = "Probability each observed zero is structural", ...)
  graphics::abline(h = 0.5, lty = 3, col = "grey50")
  invisible(szp)
}

#' Plot method for fitted dynamic models
#'
#' A convenience wrapper that dispatches to the dedicated plotting functions.
#'
#' @param x A `"dynamic_fit"` object.
#' @param which One of `"fitted"` (default), `"latent"`, `"forecast"`,
#'   `"zeros"`.
#' @param ... Passed to the underlying plotting function.
#' @return Invisibly, the result of the underlying plotting function.
#' @export
plot.dynamic_fit <- function(x, which = c("fitted", "latent",
                                          "forecast", "zeros"), ...) {
  which <- match.arg(which)
  switch(which,
         fitted   = plot_fitted(x, ...),
         latent   = plot_latent(x, ...),
         forecast = plot_forecast(x, ...),
         zeros    = plot_zero_inflation(x, ...))
}
