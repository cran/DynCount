#' Specify priors for a dynamic count / binomial model
#'
#' Builds the prior hyperparameters used by [fit_dynamic_model()]. Called with
#' no arguments it returns weakly informative defaults.
#'
#' @details
#' The model places a GMRF latent process on the series \eqn{z_t} (a log-rate
#' for the Poisson family, a logit for the binomial family): either a first-order
#' random walk (`latent_dynamics = "rw"`) or an AR(1) process
#' (`latent_dynamics = "ar1"`). The increments
#' \eqn{\varepsilon_t = z_t - \mu - \rho\, z_{t-1}} (with \eqn{\rho = 1} for the
#' random walk and \eqn{\mu = 0} unless a drift/intercept is included) are given
#' one of four innovation distributions, each governed by some of the priors
#' below.
#'
#' \strong{Innovation variance (all structures).} The baseline increment
#' variance \eqn{\sigma^2} (or, for the `"t"`/`"mixture"` structures, the
#' overall scale) has an inverse-gamma prior
#' \deqn{\sigma^2 \sim \mathrm{InvGamma}(\code{var\_shape}, \code{var\_rate}).}
#' Smaller `var_rate` gives a smoother (more strongly shrunk) latent path.
#'
#' \strong{Student-t degrees of freedom (`innovations = "t"`).} The degrees of
#' freedom are modelled as \eqn{\nu = \code{df\_min} + E}, where
#' \eqn{E \sim \mathrm{Exp}(\mathrm{rate} = 1/\code{df\_mean\_excess})}. Thus
#' \eqn{\nu \ge \code{df\_min}} and its prior mean is
#' \eqn{\code{df\_min} + \code{df\_mean\_excess}}. Large \eqn{\nu} approaches
#' the Gaussian case.
#'
#' \strong{Scale mixture (`innovations = "mixture"`).} The mixture uses
#' `mix_components` variance components, each with an
#' \eqn{\mathrm{InvGamma}(\code{var\_shape}, \code{var\_rate})} prior, and
#' symmetric Dirichlet weights with concentration `mix_concentration`.
#'
#' \strong{Stochastic volatility (`innovations = "sv"`).} The log-variance
#' AR(1) priors are delegated to \pkg{stochvol}. Pass a prior specification
#' created with [stochvol::specify_priors()] via `sv_prior`, or leave it `NULL`
#' to use the \pkg{stochvol} defaults.
#'
#' \strong{Zero inflation (`zeros = "inflated"`).} The probability that the
#' observation "gate" is open (i.e. that a zero is an ordinary sampling zero
#' rather than a structural zero) has a
#' \eqn{\mathrm{Beta}(\code{zi\_open\_a}, \code{zi\_open\_b})} prior. The
#' default `Beta(1, 1)` is uniform.
#'
#' \strong{AR(1) coefficient (`latent_dynamics = "ar1"`).} When the latent
#' state follows \eqn{z_t = \mu + \rho\, z_{t-1} + \varepsilon_t}, the
#' coefficient \eqn{\rho} is given a Gaussian prior
#' \eqn{\rho \sim \mathrm{N}(\code{ar\_rho\_mean}, \code{ar\_rho\_sd}^2)},
#' truncated to the stationary region \eqn{\rho \in (-1, 1)}, and is sampled
#' jointly with \eqn{\mu} by an exact conjugate Gibbs draw (see
#' [DynCount-package]). AR(1) always carries an intercept
#' (`include_mu = TRUE`). Under `latent_dynamics = "rw"` the coefficient is
#' fixed at \eqn{\rho = 1} and this prior is unused.
#'
#' \strong{Drift / intercept (`include_mu = TRUE`).} A scalar \eqn{\mu} in the
#' state equation \eqn{z_t = \mu + \rho\, z_{t-1} + \varepsilon_t}: a
#' \emph{drift} under the random walk (\eqn{\rho = 1}) and an \emph{intercept}
#' under AR(1) (where it is always enabled). It is given a Gaussian prior
#' \eqn{\mu \sim \mathrm{N}(\code{mu\_mean}, \code{mu\_sd}^2)}. With
#' `include_mu = FALSE` (random walk only), \eqn{\mu = 0} and is not sampled.
#'
#' \strong{Initial state.} A proper, fixed
#' \eqn{\mathrm{N}(\code{init\_mean}, \code{init\_var})} prior (default
#' \eqn{N(0, 100)}) anchors the otherwise-improper GMRF on the first latent
#' state, under both `"rw"` and `"ar1"` dynamics; keeping it free of
#' \eqn{(\mu, \rho)} is what makes their update conjugate (see
#' [DynCount-package]). With the diffuse default the initial state is
#' effectively determined by the data.
#'
#' @param var_shape Shape of the inverse-gamma prior on the innovation
#'   variance/scale. Default `2.5`.
#' @param var_rate Rate of the inverse-gamma prior on the innovation
#'   variance/scale. Default `0.5`.
#' @param df_min Lower bound for the Student-t degrees of freedom. Default `3`.
#' @param df_mean_excess Prior mean of \eqn{\nu - \code{df\_min}}. Default `6`.
#' @param mix_components Number of components in the scale-mixture innovation
#'   structure. Default `2`.
#' @param mix_concentration Symmetric Dirichlet concentration for the mixture
#'   weights. Default `1`.
#' @param sv_prior Optional \pkg{stochvol} prior specification for the
#'   stochastic-volatility innovation structure. Default `NULL`.
#' @param zi_open_a,zi_open_b Beta prior parameters for the gate-open
#'   probability in zero-inflated models. Default `1` and `1`.
#' @param ar_rho_mean,ar_rho_sd Mean and standard deviation of the Gaussian
#'   prior on the AR(1) coefficient \eqn{\rho} (used only when
#'   `latent_dynamics = "ar1"`). The conjugate Gibbs update truncates
#'   \eqn{\rho} to the stationary region \eqn{\rho \in (-1, 1)}. Defaults `0`
#'   and `1`.
#' @param mu_mean,mu_sd Mean and standard deviation of the Gaussian prior on the
#'   drift/intercept \eqn{\mu} (used only when `include_mu = TRUE`). Defaults
#'   `0` and `1`.
#' @param init_mean,init_var Mean and \strong{variance} of the Gaussian prior on
#'   the initial latent state, used under both random-walk and AR(1) dynamics.
#'   Defaults `0` and `100` (a diffuse \eqn{N(0, 100)}).
#'
#' @return An object of class `"dynamic_prior"`: a named list of hyperparameters.
#'
#' @examples
#' # Defaults
#' dynamic_prior()
#'
#' # Smoother latent path and heavier-tailed t innovations
#' dynamic_prior(var_rate = 0.1, df_min = 2, df_mean_excess = 3)
#'
#' @seealso [fit_dynamic_model()]
#' @export
dynamic_prior <- function(var_shape = 2.5,
                          var_rate = 0.5,
                          df_min = 3,
                          df_mean_excess = 6,
                          mix_components = 2,
                          mix_concentration = 1,
                          sv_prior = NULL,
                          zi_open_a = 1,
                          zi_open_b = 1,
                          ar_rho_mean = 0,
                          ar_rho_sd = 1,
                          mu_mean = 0,
                          mu_sd = 1,
                          init_mean = 0,
                          init_var = 100) {
  stopifnot(
    var_shape > 0, var_rate > 0,
    df_min > 0, df_mean_excess > 0,
    mix_components >= 1, mix_components == round(mix_components),
    mix_concentration > 0,
    zi_open_a > 0, zi_open_b > 0,
    is.finite(ar_rho_mean), ar_rho_sd > 0,
    is.finite(mu_mean), mu_sd > 0,
    is.finite(init_mean), init_var > 0
  )
  structure(
    list(
      var_shape = var_shape,
      var_rate = var_rate,
      df_min = df_min,
      df_mean_excess = df_mean_excess,
      mix_components = as.integer(mix_components),
      mix_concentration = mix_concentration,
      sv_prior = sv_prior,
      zi_open_a = zi_open_a,
      zi_open_b = zi_open_b,
      ar_rho_mean = ar_rho_mean,
      ar_rho_sd = ar_rho_sd,
      mu_mean = mu_mean,
      mu_sd = mu_sd,
      init_mean = init_mean,
      init_var = init_var
    ),
    class = "dynamic_prior"
  )
}

#' @export
print.dynamic_prior <- function(x, ...) {
  cat("<dynamic_prior>\n")
  cat(sprintf("  innovation variance ~ InvGamma(shape = %g, rate = %g)\n",
              x$var_shape, x$var_rate))
  cat(sprintf("  t degrees of freedom = %g + Exp(mean = %g)\n",
              x$df_min, x$df_mean_excess))
  cat(sprintf("  mixture: %d components, Dirichlet concentration = %g\n",
              x$mix_components, x$mix_concentration))
  cat(sprintf("  zero-inflation gate-open prob ~ Beta(%g, %g)\n",
              x$zi_open_a, x$zi_open_b))
  cat(sprintf("  AR(1) rho ~ N(mean = %g, sd = %g) truncated to (-1, 1)  [ar1 only]\n",
              x$ar_rho_mean, x$ar_rho_sd))
  cat(sprintf("  drift/intercept mu ~ N(mean = %g, sd = %g)  [include_mu only]\n",
              x$mu_mean, x$mu_sd))
  cat(sprintf("  initial state ~ N(%g, %g)  [rw and ar1]\n",
              x$init_mean, x$init_var))
  cat(sprintf("  sv_prior: %s\n", if (is.null(x$sv_prior)) "stochvol defaults" else "user supplied"))
  invisible(x)
}
