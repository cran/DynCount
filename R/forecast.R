# Forecasting -----------------------------------------------------------------
#
# Forecasts are produced after sampling by forward simulation. Future latent
# states carry no likelihood, so given one posterior draw of the parameters
# and of the last in-sample state z_n, their joint conditional is simply the
# transition law of the latent process. Simulating
#   z_{n+h} = mu + rho z_{n+h-1} + eps_{n+h},   h = 1, ..., H,
# with eps drawn from the fitted innovation structure, and then drawing the
# response from the observation model, gives exact draws from the posterior
# predictive distribution. This is the same distribution as the one obtained
# by appending the H future states to the latent vector and sampling them as
# missing values inside the MCMC, but it does not require the horizon to be
# fixed before fitting.

#' @importFrom generics forecast
#' @export
generics::forecast

#' Forecast a fitted dynamic model
#'
#' Returns posterior predictive forecasts for the `H` periods after the end of
#' the series. The forecasts are obtained by forward simulation from the
#' stored posterior draws, so any horizon can be requested after fitting. If
#' the model was fitted with `horizon >= 1`, the forecast draws stored in the
#' fit are returned unless `horizon`, `forecast_offset` or `forecast_trials`
#' is supplied.
#'
#' @details
#' For every stored posterior draw, the latent path is propagated from the
#' last in-sample state \eqn{z_n} with
#' \eqn{z_{n+h} = \mu + \rho z_{n+h-1} + \varepsilon_{n+h}}, drawing the
#' increments from the fitted innovation structure (Gaussian with the drawn
#' variance, Student-t with the drawn scale and degrees of freedom, the drawn
#' scale mixture, or the stochastic-volatility process continued from its last
#' in-sample log-variance). A response is then drawn from the observation
#' model at each simulated state. The intervals therefore reflect parameter,
#' state and innovation uncertainty.
#'
#' \strong{Zero inflation.} For fits with `zeros = "inflated"` the response
#' forecast is \strong{unconditional}: each draw includes the zero-inflation
#' gate, so structural zeros are reproduced and the draws are directly
#' comparable to future observations (like the in-sample `yrep`, unlike
#' `yrep_open`). The latent forecast summary in `$latent` is gate-free by
#' construction.
#'
#' \strong{Multinomial family.} Response forecasts are multinomial draws with
#' the totals given by `forecast_trials`. The summaries are in long format
#' with a `category` column. `summary` and `final` cover all `K` categories on
#' the count scale, `prob` gives the forecast category shares, and `latent`
#' the `K - 1` ALR series. `draws` is a `draws x H x K` array and
#' `final_draws` a `draws x K` matrix.
#'
#' @param object A `"dynamic_fit"` object.
#' @param horizon Forecast horizon `H` (a positive integer). If `NULL`
#'   (default), the forecast stored in the fit is returned; the fit must then
#'   have been created with `horizon >= 1`.
#' @param forecast_offset Known offset over the forecast horizon, in the same
#'   format as in [fit_dynamic_model()]. Defaults to `0`; a warning is issued
#'   if the model was fitted with a non-zero offset but no `forecast_offset`
#'   is supplied.
#' @param forecast_trials Binomial and multinomial families: the number of
#'   trials (binomial) or the total count per period (multinomial) over the
#'   forecast horizon, of length 1 (recycled) or `horizon`. Defaults to the
#'   last observed number of trials (binomial) or the last non-zero row total
#'   (multinomial).
#' @param probs Quantile probabilities for the summary. Default
#'   `c(0.025, 0.5, 0.975)`.
#' @param seed Optional random seed for the forward simulation. The previous
#'   state of the global random number generator is restored afterwards.
#' @param ... Unused.
#'
#' @return An object of class `"dynamic_forecast"`: a list with `summary` (one
#'   row per horizon \eqn{1, \dots, H}, on the response scale), `latent`
#'   (summary of the forecast latent path), `draws` (the predictive draws,
#'   draws x `H`), `latent_draws` (the latent forecast draws), `final` (the
#'   single-row summary of the final `H`-step-ahead forecast), `final_draws`
#'   (the predictive draws at horizon `H`), the forecast `horizon`, and the
#'   `forecast_offset` / `forecast_trials` that were used. For the
#'   multinomial family additionally `prob` (share summary) and `prob_draws`
#'   (`draws x H x K` share draws); see Details for the layout.
#'
#' @examples
#' sim <- simulate_dynamic_poisson(60, 0.2, 2, seed = 1)
#' fit <- fit_dynamic_model(sim$y, nsave = 300, nburn = 200, seed = 1)
#' fc <- forecast(fit, horizon = 8)
#' fc$summary
#' fc$final
#' @seealso [fit_dynamic_model()], [plot_forecast()]
#' @export
forecast.dynamic_fit <- function(object, horizon = NULL, forecast_offset = NULL,
                                 forecast_trials = NULL,
                                 probs = c(0.025, 0.5, 0.975), seed = NULL, ...) {
  d  <- object$draws
  sp <- object$spec
  fam <- sp$family
  is_multi <- identical(fam, "multinomial")
  stored <- !is.null(d$forecast_y) && isTRUE(sp$horizon >= 1L)
  new_sim <- !is.null(horizon) || !is.null(forecast_offset) ||
    !is.null(forecast_trials) || !stored

  if (new_sim) {
    H <- horizon %||% (if (stored) sp$horizon else NULL)
    if (is.null(H)) {
      stop("No forecast is stored in this fit. Supply a horizon, e.g. ",
           "`forecast(fit, horizon = 8)`.", call. = FALSE)
    }
    H <- check_count(H, "horizon")
    if (is.null(d$z0) || (sp$innovations != "sv" && is.null(d$scale))) {
      stop("This fit does not contain the draws needed for forward simulation; ",
           "refit it with the installed version of DynCount.", call. = FALSE)
    }
    J <- if (is_multi) object$data$K - 1L else 1L
    fi <- resolve_forecast_inputs(fam, H, J, object$data$trials,
                                  has_offset = any(object$data$offset != 0),
                                  forecast_offset = forecast_offset,
                                  forecast_trials = forecast_trials)
    sims <- with_seed(seed, simulate_forecast(object, H, fi$offset, fi$trials))
    ydraw  <- sims$y
    latent <- sims$z
    pdraw  <- sims$prob
    h <- H
    used_offset <- if (is_multi) fi$offset else fi$offset[, 1L]
    used_trials <- fi$trials
  } else {
    ydraw  <- d$forecast_y
    latent <- d$forecast_z
    pdraw  <- d$forecast_prob
    h <- sp$horizon
    used_offset <- sp$forecast_offset
    used_trials <- sp$forecast_trials
  }

  if (is_multi) {
    resp <- summarise_draws_long(ydraw, probs, "horizon")
    lat  <- summarise_draws_long(latent, probs, "horizon")
    prob <- summarise_draws_long(pdraw, probs, "horizon")
    final <- resp[resp$horizon == h, , drop = FALSE]
    rownames(final) <- NULL
    final_draws <- matrix(ydraw[, h, ], dim(ydraw)[1L], dim(ydraw)[3L],
                          dimnames = list(NULL, dimnames(ydraw)[[3L]]))
  } else {
    resp <- summarise_draws(ydraw, probs)
    resp <- cbind(horizon = seq_len(h), resp)
    lat  <- summarise_draws(latent, probs)
    lat  <- cbind(horizon = seq_len(h), lat)
    prob <- NULL
    final <- resp[h, , drop = FALSE]          # the final H-step-ahead row
    rownames(final) <- NULL
    final_draws <- ydraw[, h]                 # response draws at horizon H
  }

  structure(
    list(summary = resp, latent = lat, prob = prob, draws = ydraw,
         latent_draws = latent,
         prob_draws = if (is_multi) pdraw else NULL,
         final = final, final_draws = final_draws,
         horizon = h, forecast_offset = used_offset, forecast_trials = used_trials,
         family = fam, n_obs = object$data$n,
         categories = object$data$categories, baseline = object$data$baseline),
    class = "dynamic_forecast"
  )
}

#' @export
print.dynamic_forecast <- function(x, ...) {
  h <- x$horizon
  cat(sprintf("<dynamic_forecast> %s, horizon %d\n", x$family, h))
  if (identical(x$family, "multinomial")) {
    # show the first horizons of every category, then the final step
    hh <- min(h, 5L)
    print(x$summary[x$summary$horizon <= hh, , drop = FALSE], row.names = FALSE)
    if (h > hh) {
      cat(sprintf("final %d-step-ahead forecast:\n", h))
      print(x$final, row.names = FALSE)
    }
    return(invisible(x))
  }
  print(utils::head(x$summary, 10), row.names = FALSE)
  if (h > 10L) {
    # the final row is not visible above, so show it separately
    cat(sprintf("final %d-step-ahead forecast:\n", h))
    print(x$final, row.names = FALSE)
  }
  invisible(x)
}

#' Validate and expand the forecast-period inputs
#'
#' Returns `offset`, an `H x J` matrix, and `trials`, a length-`H` vector
#' (`NULL` for the Poisson family).
#' @keywords internal
#' @noRd
resolve_forecast_inputs <- function(family, H, J, trials, has_offset,
                                    forecast_offset = NULL,
                                    forecast_trials = NULL) {
  is_multi <- identical(family, "multinomial")
  if (!is.null(forecast_offset) &&
      (!is.numeric(forecast_offset) || anyNA(forecast_offset) ||
       any(!is.finite(forecast_offset)))) {
    stop("`forecast_offset` must be a finite numeric vector or matrix.", call. = FALSE)
  }
  if (is.null(forecast_offset)) {
    if (isTRUE(has_offset)) {
      warning("The model was fitted with a non-zero `offset` but no ",
              "`forecast_offset` was supplied; the forecasts assume an offset of 0.",
              call. = FALSE)
    }
    forecast_offset <- 0
  }
  if (is_multi) {
    fo <- expand_alr_offset(forecast_offset, H, J, "forecast_offset", "H")
    if (length(fo) == 1L) fo <- matrix(fo, H, J)
  } else {
    fo <- as.numeric(forecast_offset)
    if (length(fo) == 1L) fo <- rep(fo, H)
    if (length(fo) != H) {
      stop("`forecast_offset` must have length 1 or `horizon`.", call. = FALSE)
    }
    fo <- matrix(fo, H, 1L)
  }

  ft <- NULL
  if (identical(family, "poisson")) {
    if (!is.null(forecast_trials)) {
      warning("`forecast_trials` is ignored for the Poisson family.", call. = FALSE)
    }
  } else if (!is.null(forecast_trials)) {
    if (!is.numeric(forecast_trials) || anyNA(forecast_trials) ||
        any(forecast_trials < 0) || any(forecast_trials != round(forecast_trials)) ||
        (identical(family, "binomial") && any(forecast_trials <= 0))) {
      stop("`forecast_trials` must be non-negative integers (positive for the ",
           "binomial family).", call. = FALSE)
    }
    ft <- as.numeric(forecast_trials)
    if (length(ft) == 1L) ft <- rep(ft, H)
    if (length(ft) != H) {
      stop("`forecast_trials` must have length 1 or `horizon`.", call. = FALSE)
    }
  } else {
    # last observed number of trials / last non-zero row total
    last <- utils::tail(which(trials > 0), 1L)
    ft <- rep(as.numeric(trials[last]), H)
  }
  list(offset = fo, trials = ft)
}

#' Draw a mixture component per row of a weight matrix
#' @keywords internal
#' @noRd
draw_component <- function(w) {
  S <- nrow(w); Hm <- ncol(w)
  if (Hm == 1L) return(rep(1L, S))
  cw <- t(apply(w, 1L, cumsum))
  1L + as.integer(rowSums(stats::runif(S) > cw[, -Hm, drop = FALSE]))
}

#' Forward-simulate latent and response forecasts from the stored draws
#'
#' @param object a `"dynamic_fit"`.
#' @param H forecast horizon.
#' @param forecast_offset `H x J` matrix of offsets.
#' @param forecast_trials length-`H` trials / totals (`NULL` for Poisson).
#' @return list with `z` (latent draws), `y` (response draws) and `prob`
#'   (multinomial shares, otherwise `NULL`), in the layouts documented in
#'   [fit_dynamic_model()].
#' @keywords internal
#' @noRd
simulate_forecast <- function(object, H, forecast_offset, forecast_trials) {
  d  <- object$draws
  sp <- object$spec
  is_multi <- identical(sp$family, "multinomial")
  n  <- object$data$n
  as3 <- function(a) if (length(dim(a)) == 3L) a else array(a, c(dim(a), 1L))
  asm <- function(v) if (is.null(v) || !is.null(dim(v))) v else matrix(v, ncol = 1L)
  z     <- as3(d$z)
  S <- dim(z)[1L]; J <- dim(z)[3L]
  rho   <- asm(d$rho);   mu    <- asm(d$mu)
  scale <- asm(d$scale); nu    <- asm(d$nu)
  svmu  <- asm(d$sv_mu); svphi <- asm(d$sv_phi); svsig <- asm(d$sv_sigma)
  mixw  <- if (!is.null(d$mix_weight)) as3(d$mix_weight)
  mixv  <- if (!is.null(d$mix_var)) as3(d$mix_var)
  inn   <- sp$innovations

  Zf <- array(NA_real_, c(S, H, J))
  for (k in seq_len(J)) {
    zc <- z[, n, k]                                  # last in-sample state
    if (inn == "sv") lv <- log(as3(d$sig2)[, n, k])  # last in-sample log-variance
    for (h in seq_len(H)) {
      v <- switch(inn,
        gaussian = scale[, k],
        t = scale[, k] / stats::rgamma(S, shape = nu[, k] / 2, rate = nu[, k] / 2),
        mixture = {
          vv <- matrix(mixv[, , k], S)
          vv[cbind(seq_len(S), draw_component(matrix(mixw[, , k], S)))]
        },
        sv = {
          lv <- svmu[, k] + svphi[, k] * (lv - svmu[, k]) + svsig[, k] * stats::rnorm(S)
          exp(lv)
        })
      zc <- mu[, k] + rho[, k] * zc + stats::rnorm(S, 0, sqrt(v))
      Zf[, h, k] <- zc
    }
  }

  # linear predictor over the horizon (forecast_offset is H x J)
  Eta_f <- Zf + rep(forecast_offset, each = S)
  if (is_multi) {
    K <- object$data$K
    b <- object$data$baseline
    cats <- object$data$categories
    Pf <- Yf <- array(NA_real_, c(S, H, K))
    for (h in seq_len(H)) {
      Ph <- alr_to_prob(matrix(Eta_f[, h, ], S, J), b)
      Pf[, h, ] <- Ph
      Yf[, h, ] <- rmultinom_vec(rep(forecast_trials[h], S), Ph)
    }
    dimnames(Zf) <- list(NULL, NULL, cats[-b])
    dimnames(Pf) <- dimnames(Yf) <- list(NULL, NULL, cats)
    return(list(z = Zf, y = Yf, prob = Pf))
  }

  eta <- matrix(Eta_f[, , 1L], S, H)
  Yf <- if (identical(sp$family, "poisson")) {
    stats::rpois(S * H, exp(eta))
  } else {
    stats::rbinom(S * H, size = rep(forecast_trials, each = S),
                  prob = stats::plogis(eta))
  }
  Yf <- matrix(Yf, S, H)
  if (identical(sp$zeros, "inflated")) {
    # gate closed -> structural zero; pi_open is recycled down the columns
    closed <- matrix(stats::runif(S * H), S, H) > d$pi_open
    Yf[closed] <- 0
  }
  list(z = matrix(Zf[, , 1L], S, H), y = Yf, prob = NULL)
}
