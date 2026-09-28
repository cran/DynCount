# Prediction and forecasting --------------------------------------------------

#' In-sample fitted values and posterior predictive replicates
#'
#' @param object A `"dynamic_fit"` object.
#' @param type `"mean"` (default) returns the posterior of the mean of `y`;
#'   `"response"` returns posterior predictive replicates of `y`. For the
#'   multinomial family `"prob"` returns the posterior of the category
#'   shares \eqn{p_{t,k}} (the mean is then the expected count
#'   \eqn{N_t p_{t,k}}).
#' @param conditional Logical. Leave `FALSE` (the default) for anything
#'   compared against observed data -- posterior predictive checks,
#'   calibration -- and set `TRUE` only to inspect the latent intensity
#'   process. With `FALSE`, means and replicates come from the full
#'   zero-inflated model (the gate is included, so structural zeros are
#'   reproduced); with `TRUE`, they are conditional on the gate being open,
#'   i.e. drawn straight from the Poisson/binomial observation model, and show
#'   systematically too few zeros under zero inflation. Without zero inflation
#'   the two versions are identical, and the argument is ignored for the
#'   multinomial family (which has no gate). See [DynCount-package] for the
#'   gate notation.
#' @param probs Quantile probabilities for the summary. Default
#'   `c(0.025, 0.5, 0.975)`.
#' @param ... Unused.
#'
#' @return A list with `summary` (a data frame, one row per observation) and
#'   `draws` (the underlying draws matrix, draws x time). For the multinomial
#'   family `summary` is in long format with one row per observation and
#'   category (columns `time`, `category`, `observed`, ...) and `draws` is a
#'   `draws x time x K` array.
#' @examples
#' sim <- simulate_dynamic_poisson(40, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 200, nburn = 100, seed = 1)
#' head(predict(fit)$summary)                     # posterior of the mean of y
#' head(predict(fit, type = "response")$summary)  # posterior predictive
#' @export
predict.dynamic_fit <- function(object, type = c("mean", "response", "prob"),
                                conditional = FALSE,
                                probs = c(0.025, 0.5, 0.975), ...) {
  type <- match.arg(type)
  stopifnot(is.logical(conditional), length(conditional) == 1L)
  is_multi <- identical(object$spec$family, "multinomial")
  if (type == "prob" && !is_multi) {
    stop("`type = \"prob\"` is only available for the multinomial family.",
         call. = FALSE)
  }
  d <- object$draws
  draws <- switch(type,
    mean     = if (conditional) d$fitted_open %||% d$fitted else d$fitted,
    response = if (conditional) d$yrep_open %||% d$yrep else d$yrep,
    prob     = d$fitted_prob)
  if (is_multi) {
    s <- summarise_draws_long(draws, probs, "time")
    obs <- object$data$y
    if (type == "prob") obs <- obs / pmax(object$data$trials, 1)
    s$observed <- obs[cbind(s$time, match(s$category, object$data$categories))]
    s <- s[, c("time", "category", "observed",
               setdiff(names(s), c("time", "category", "observed")))]
  } else {
    s <- summarise_draws(draws, probs)
    s <- cbind(time = seq_len(nrow(s)), observed = object$data$y, s)
  }
  list(summary = s, draws = draws)
}
