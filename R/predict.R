# Prediction and forecasting --------------------------------------------------

#' In-sample fitted values and posterior predictive replicates
#'
#' @param object A `"dynamic_fit"` object.
#' @param type `"mean"` (default) returns the posterior of the mean of `y`;
#'   `"response"` returns posterior predictive replicates of `y`.
#' @param conditional Logical. Leave `FALSE` (the default) for anything
#'   compared against observed data -- posterior predictive checks,
#'   calibration -- and set `TRUE` only to inspect the latent intensity
#'   process. With `FALSE`, means and replicates come from the full
#'   zero-inflated model (the gate is included, so structural zeros are
#'   reproduced); with `TRUE`, they are conditional on the gate being open,
#'   i.e. drawn straight from the Poisson/binomial observation model, and show
#'   systematically too few zeros under zero inflation. Without zero inflation
#'   the two versions are identical. See [DynCount-package] for the gate
#'   notation.
#' @param probs Quantile probabilities for the summary. Default
#'   `c(0.025, 0.5, 0.975)`.
#' @param ... Unused.
#'
#' @return A list with `summary` (a data frame, one row per observation) and
#'   `draws` (the underlying draws matrix, draws x time).
#' @examples
#' sim <- simulate_dynamic_poisson(40, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 200, nburn = 100)
#' head(predict(fit)$summary)
#' @export
predict.dynamic_fit <- function(object, type = c("mean", "response"),
                                conditional = FALSE,
                                probs = c(0.025, 0.5, 0.975), ...) {
  type <- match.arg(type)
  stopifnot(is.logical(conditional), length(conditional) == 1L)
  draws <- if (type == "mean") {
    if (conditional) object$draws$fitted_open else object$draws$fitted
  } else {
    if (conditional) object$draws$yrep_open else object$draws$yrep
  }
  s <- summarise_draws(draws, probs)
  s <- cbind(time = seq_len(nrow(s)), observed = object$data$y, s)
  list(summary = s, draws = draws)
}

#' Generic forecast function
#'
#' Generic for extracting forecasts from a fitted model. See
#' [forecast.dynamic_fit()] for the method for `"dynamic_fit"` objects.
#'
#' @param object A fitted model object.
#' @param ... Passed to methods.
#' @return The value returned by the method dispatched on `object`; for
#'   `"dynamic_fit"` objects a `"dynamic_forecast"` object (see
#'   [forecast.dynamic_fit()]).
#' @export
forecast <- function(object, ...) UseMethod("forecast")

#' Forecast a fitted dynamic model
#'
#' Returns the forecast draws that were produced during MCMC sampling.
#' Forecasts are generated inside the sampler when the model is fitted with
#' `horizon = H >= 1`. This method extracts and summarises those stored
#' draws.
#'
#' @param object A `"dynamic_fit"` object that was fitted with `horizon >= 1`.
#' @param probs Quantile probabilities for the summary. Default
#'   `c(0.025, 0.5, 0.975)`.
#' @param ... Unused.
#'
#' @details
#' \strong{Zero inflation.} For fits with `zeros = "inflated"` the response
#' forecast is \strong{unconditional}: each draw includes the zero-inflation
#' gate, so structural zeros are reproduced and the draws are directly
#' comparable to future observations (like the in-sample `yrep`, unlike
#' `yrep_open`). The latent forecast summary in `$latent` is gate-free by
#' construction.
#'
#' @return An object of class `"dynamic_forecast"`: a list with `summary` (one
#'   row per horizon \eqn{1, \dots, H}, on the response scale), `latent`
#'   (summary of the forecast latent path), `draws` (the predictive draws,
#'   draws x `H`), `final` (the single-row summary of the final `H`-step-ahead
#'   forecast), `final_draws` (the predictive draws at horizon `H`), and the
#'   forecast horizon `horizon`.
#'
#' @examples
#' sim <- simulate_dynamic_poisson(60, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 200, horizon = 8)
#' fc <- forecast(fit)
#' fc$summary
#' @export
forecast.dynamic_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  d  <- object$draws
  sp <- object$spec
  if (is.null(d$forecast_y) || sp$horizon < 1L) {
    stop("No forecast was produced. Refit with `horizon >= 1` to generate ",
         "forecasts during sampling.", call. = FALSE)
  }
  h      <- sp$horizon
  ydraw  <- d$forecast_y                  # draws x H response forecast
  latent <- d$forecast_z                  # draws x H latent forecast
  fam    <- sp$family

  resp <- summarise_draws(ydraw, probs)
  resp <- cbind(horizon = seq_len(h), resp)
  lat  <- summarise_draws(latent, probs)
  lat  <- cbind(horizon = seq_len(h), lat)

  structure(
    list(summary = resp, latent = lat, draws = ydraw,
         final = resp[h, , drop = FALSE],          # the final H-step-ahead row
         final_draws = ydraw[, h],                 # response draws at horizon H
         horizon = h, family = fam, n_obs = object$data$n),
    class = "dynamic_forecast"
  )
}

#' @export
print.dynamic_forecast <- function(x, ...) {
  h <- x$horizon
  cat(sprintf("<dynamic_forecast> %s, horizon %d\n", x$family, h))
  print(utils::head(x$summary, 10), row.names = FALSE)
  if (h > 10L) {
    # the final row is not visible above, so show it separately
    cat(sprintf("final %d-step-ahead forecast:\n", h))
    print(x$final, row.names = FALSE)
  }
  invisible(x)
}
