# Plotting --------------------------------------------------------------------
#
# All plots use base graphics to keep dependencies minimal. Each function draws
# a posterior median with a shaded credible band where appropriate. For the
# multinomial family the functions draw one panel per category (all
# categories by default, or the ones selected with `category`), arranging the
# panels in a grid and restoring the graphics parameters afterwards.

#' @keywords internal
#' @noRd
band <- function(x, lo, hi, col) {
  graphics::polygon(c(x, rev(x)), c(lo, rev(hi)), col = col, border = NA)
}

#' Resolve the `category` argument of the plot functions to column indices
#' among `choices` (a character vector of category labels).
#' @keywords internal
#' @noRd
resolve_categories <- function(category, choices) {
  if (is.null(category)) return(seq_along(choices))
  if (is.character(category)) {
    idx <- match(category, choices)
    if (anyNA(idx)) stop("Unknown category: ",
                         paste(category[is.na(idx)], collapse = ", "), call. = FALSE)
    return(idx)
  }
  if (!is.numeric(category) || any(category < 1) || any(category > length(choices)))
    stop("`category` must be category names or indices.", call. = FALSE)
  as.integer(category)
}

#' Call a base graphics function with defaults that `...` may override
#'
#' `defaults` and `dots` are named lists; entries in `dots` replace entries of
#' the same name in `defaults`, so users can pass e.g. `main`, `xlab` or `ylim`
#' without triggering "matched by multiple actual arguments" errors.
#' @keywords internal
#' @noRd
call_with_defaults <- function(fun, args, defaults, dots) {
  do.call(fun, c(args, utils::modifyList(defaults, dots)))
}

#' Lay out `m` panels in a grid and return a function restoring `par()`.
#' Does nothing (and returns a no-op) for a single panel.
#' @keywords internal
#' @noRd
panel_grid <- function(m) {
  if (m <= 1L) return(function() invisible(NULL))
  nc <- ceiling(sqrt(m)); nr <- ceiling(m / nc)
  op <- graphics::par(mfrow = c(nr, nc))
  function() graphics::par(op)
}

#' Plot the fitted latent trajectory
#'
#' Shows the posterior median of the latent state (log-rate for the Poisson
#' family, logit for the binomial family, one additive-log-ratio series per
#' non-baseline category for the multinomial family) with a credible band.
#'
#' @param object A `"dynamic_fit"` object.
#' @param level Credible level for the band. Default `0.95`.
#' @param category Multinomial family only: which categories to draw (names
#'   or indices among the non-baseline categories for [plot_latent()]; among
#'   all `K` categories for [plot_fitted()] and [plot_forecast()]). `NULL`
#'   (default) draws one panel per category.
#' @param ... Passed to [graphics::plot()]; may override the default `main`,
#'   `xlab`, `ylab`, `ylim` and other plot arguments.
#' @return Invisibly, the summary data frame used for plotting (in long format
#'   with a `category` column for the multinomial family).
#' @examples
#' sim <- simulate_dynamic_poisson(60, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 200, seed = 1)
#' plot_latent(fit)
#' @export
plot_latent <- function(object, level = 0.95, category = NULL, ...) {
  dots <- list(...)
  a <- (1 - level) / 2
  pr <- c(a, 0.5, 1 - a)
  fam <- object$spec$family
  if (identical(fam, "multinomial")) {
    lat_names <- object$data$categories[-object$data$baseline]
    ks <- resolve_categories(category, lat_names)
    z <- object$draws$z[, , ks, drop = FALSE]
    s <- summarise_draws_long(z, pr, "time")
    restore <- panel_grid(length(ks)); on.exit(restore())
    for (k in seq_along(ks)) {
      sk <- s[s$category == lat_names[ks[k]], , drop = FALSE]
      q <- qband(sk, a)
      call_with_defaults(graphics::plot, list(sk$time, q$md),
        list(type = "n", ylim = range(q$lo, q$hi), xlab = "time",
             ylab = sprintf("latent log-ratio vs %s", object$data$baseline_name),
             main = sprintf("Latent trajectory: %s", lat_names[ks[k]])), dots)
      band(sk$time, q$lo, q$hi, grDevices::adjustcolor("steelblue", 0.3))
      graphics::lines(sk$time, q$md, col = "steelblue", lwd = 2)
    }
    return(invisible(s))
  }
  s <- summarise_draws(object$draws$z, pr)
  q <- qband(s, a)
  tt <- seq_len(nrow(s))
  ylab <- if (fam == "poisson") "latent log-rate" else "latent logit"
  call_with_defaults(graphics::plot, list(tt, q$md),
    list(type = "n", ylim = range(q$lo, q$hi), xlab = "time", ylab = ylab,
         main = "Latent trajectory"), dots)
  band(tt, q$lo, q$hi, grDevices::adjustcolor("steelblue", 0.3))
  graphics::lines(tt, q$md, col = "steelblue", lwd = 2)
  invisible(s)
}

#' Plot observed versus fitted values
#'
#' Observed series with the posterior median fitted mean and a credible band on
#' the response scale. For the multinomial family, one panel per category
#' showing the observed and expected counts \eqn{N_t p_{t,k}}.
#'
#' @inheritParams plot_latent
#' @param level Credible level. Default `0.95`.
#' @return Invisibly, the fitted summary.
#' @examples
#' sim <- simulate_dynamic_poisson(60, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 200, seed = 1)
#' plot_fitted(fit)
#' @export
plot_fitted <- function(object, level = 0.95, category = NULL, ...) {
  dots <- list(...)
  a <- (1 - level) / 2
  pr <- c(a, 0.5, 1 - a)
  ylab <- if (identical(object$spec$family, "binomial")) "successes" else "count"
  draw_one <- function(tt, y, q, main) {
    call_with_defaults(graphics::plot, list(tt, y),
      list(type = "n", ylim = range(0, y, q$hi), xlab = "time", ylab = ylab,
           main = main), dots)
    band(tt, q$lo, q$hi, grDevices::adjustcolor("tomato", 0.25))
    graphics::lines(tt, q$md, col = "tomato", lwd = 2)
    graphics::points(tt, y, pch = 16, cex = 0.6, col = "grey20")
    graphics::legend("topleft", c("observed", "fitted (median)", "band"),
                     pch = c(16, NA, 15), lty = c(NA, 1, NA),
                     col = c("grey20", "tomato", grDevices::adjustcolor("tomato", 0.25)),
                     bty = "n", cex = 0.8)
  }
  if (identical(object$spec$family, "multinomial")) {
    cats <- object$data$categories
    ks <- resolve_categories(category, cats)
    s <- summarise_draws_long(object$draws$fitted[, , ks, drop = FALSE], pr, "time")
    restore <- panel_grid(length(ks)); on.exit(restore())
    for (k in ks) {
      sk <- s[s$category == cats[k], , drop = FALSE]
      draw_one(sk$time, object$data$y[, k], qband(sk, a),
               sprintf("Observed vs fitted: %s", cats[k]))
    }
    return(invisible(s))
  }
  s <- summarise_draws(object$draws$fitted, pr)
  y <- object$data$y
  draw_one(seq_along(y), y, qband(s, a), "Observed vs fitted")
  invisible(s)
}

#' Plot observed history, fitted values and forecast with uncertainty
#'
#' Shows the observed series together with the in-sample fitted mean (median and
#' credible band) and the forecast (median and credible band) over the forecast
#' horizon. The forecast is the one stored in the fit, or a new one computed
#' with [forecast.dynamic_fit()] when `horizon` (or another forecast input) is
#' supplied. For the multinomial family, one panel per category.
#'
#' @param object A `"dynamic_fit"` object.
#' @param level Credible level. Default `0.95`.
#' @inheritParams plot_latent
#' @param horizon,forecast_offset,forecast_trials,seed Passed to
#'   [forecast.dynamic_fit()]. With the defaults, the forecast stored in the
#'   fit is plotted.
#' @param ... Passed to [graphics::plot()]; may override the default `main`,
#'   `xlab`, `ylab`, `ylim` and other plot arguments.
#' @return Invisibly, the `"dynamic_forecast"` object.
#' @examples
#' sim <- simulate_dynamic_poisson(60, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 200, seed = 1)
#' plot_forecast(fit, horizon = 12)
#' @export
plot_forecast <- function(object, level = 0.95, category = NULL, horizon = NULL,
                          forecast_offset = NULL, forecast_trials = NULL,
                          seed = NULL, ...) {
  dots <- list(...)
  a  <- (1 - level) / 2
  pr <- c(a, 0.5, 1 - a)
  fc <- forecast(object, horizon = horizon, forecast_offset = forecast_offset,
                 forecast_trials = forecast_trials, probs = pr, seed = seed)
  n <- object$data$n
  h <- fc$horizon
  fh <- (n + 1):(n + h)
  tt <- seq_len(n)
  ylab <- if (identical(object$spec$family, "binomial")) "successes" else "count"
  draw_one <- function(y, qf, qc, main) {
    ylim <- range(0, y, qf$lo, qf$hi, qc$lo, qc$hi)
    call_with_defaults(graphics::plot, list(tt, y),
      list(type = "n", xlim = c(1, n + h), ylim = ylim, xlab = "time",
           ylab = ylab, main = main), dots)
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
  }
  if (identical(object$spec$family, "multinomial")) {
    cats <- object$data$categories
    ks <- resolve_categories(category, cats)
    fit_s <- summarise_draws_long(object$draws$fitted, pr, "time")
    restore <- panel_grid(length(ks)); on.exit(restore())
    for (k in ks) {
      qf <- qband(fit_s[fit_s$category == cats[k], , drop = FALSE], a)
      qc <- qband(fc$summary[fc$summary$category == cats[k], , drop = FALSE], a)
      draw_one(object$data$y[, k], qf, qc,
               sprintf("Observed, fitted and forecast: %s", cats[k]))
    }
    return(invisible(fc))
  }
  fit_s <- summarise_draws(object$draws$fitted, pr)   # in-sample fitted band
  draw_one(object$data$y, qband(fit_s, a), qband(fc$summary, a),
           "Observed, fitted and forecast")
  invisible(fc)
}

#' Plot zero-inflation diagnostics
#'
#' Bar chart of the posterior probability that each observed zero is a
#' structural zero (see [structural_zero_prob()]), one bar per observed zero
#' labelled by its time index, with a dotted reference line at 0.5. Requires
#' a model fitted with `zeros = "inflated"` (Poisson or binomial family).
#'
#' @param object A `"dynamic_fit"` object.
#' @param ... Passed to [graphics::barplot()]; may override the defaults
#'   (e.g. `main`, `col`).
#' @return Invisibly, the [structural_zero_prob()] data frame.
#' @examples
#' sim <- simulate_dynamic_poisson(80, 0.2, 2, zero_inflation = 0.3, seed = 1)
#' fit <- fit_dynamic_model(sim$y, zero_inflation = TRUE, nsave = 300, nburn = 200,
#'                          seed = 1)
#' plot_zero_inflation(fit)
#' @export
plot_zero_inflation <- function(object, ...) {
  szp <- structural_zero_prob(object, zeros_only = TRUE)
  if (nrow(szp) == 0) {
    message("No observed zeros to display.")
    return(invisible(szp))
  }
  call_with_defaults(graphics::barplot, list(szp$p_structural),
    list(names.arg = szp$time, las = 2,
         col = grDevices::adjustcolor("darkorange", 0.7),
         border = NA, ylim = c(0, 1),
         ylab = "P(structural zero)", xlab = "time index",
         main = "Probability each observed zero is structural"), list(...))
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
#' @param ... Passed to the underlying plotting function (e.g. `category`
#'   for multinomial fits, or `horizon` for `which = "forecast"`).
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
