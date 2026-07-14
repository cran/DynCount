#' DynCount: Bayesian Dynamic Models for Poisson and Binomial Time Series
#'
#' The \pkg{DynCount} package fits state-space models for non-Gaussian
#' time series. A latent trajectory \eqn{z_t} follows flexible dynamics --
#' a first-order random walk or a stationary AR(1) process -- and the
#' observations are linked to it through either a Poisson (log link) or a
#' binomial (logit link) observation model. An optional known `offset` may be
#' added to the observation model's linear predictor (a log-exposure for the
#' Poisson mean, a logit shift for the binomial probability); it is a fixed,
#' user-supplied input and is not part of the latent process \eqn{z_t}.
#'
#' The package implements and extends the methodology of Zens and Bijak
#' (2026, \doi{10.1214/26-AOAS2171}). It supports several innovation structures (Gaussian,
#' Student-t, finite scale mixture, stochastic volatility) and time-constant zero
#' inflation.
#'
#' @section Main entry points:
#' \describe{
#'   \item{[fit_dynamic_model()]}{Fit a Poisson or binomial dynamic model
#'         (random walk or stationary AR(1)), optionally forecasting during
#'         sampling via `horizon`.}
#'   \item{[forecast()] / [predict()]}{Extract the forecast draws produced
#'         during sampling, or the in-sample fitted values.}
#'   \item{[summary()]}{Posterior summaries of a fitted model.}
#'   \item{[simulate_dynamic_poisson()], [simulate_dynamic_binomial()]}{Simulate
#'         data.}
#'   \item{Plotting}{[plot_latent()], [plot_fitted()], [plot_forecast()],
#'         [plot_zero_inflation()].}
#' }
#'
#' @section Latent dynamics:
#' The `latent_dynamics` argument of [fit_dynamic_model()] selects the GMRF
#' state evolution \eqn{z_t = \mu + \rho z_{t-1} + \varepsilon_t}:
#' \describe{
#'   \item{`"rw"`}{A first-order random walk, i.e. \eqn{\rho = 1} fixed (the
#'         default).}
#'   \item{`"ar1"`}{A stationary AR(1) process, always with an intercept
#'         (`include_mu = TRUE`).}
#' }
#' Both share the precision \eqn{P = D_\rho^\top \mathrm{diag}(1/\sigma^2)
#' D_\rho}, where \eqn{D_\rho} is the generalised first-difference operator; it
#' is symmetric tridiagonal, so its bands are formed and the latent full
#' conditionals evaluated in \eqn{O(N)}. A scalar \eqn{\mu} (set
#' `include_mu = TRUE`) adds a \emph{drift} (RW) or \emph{intercept} (AR(1)) to
#' the state equation; it enters as a linear term in the latent full
#' conditionals, leaving the sparse precision unchanged.
#'
#' The otherwise-improper GMRF is anchored by a proper, fixed
#' \eqn{N(\code{init\_mean}, \code{init\_var})} prior on the first latent
#' state (default \eqn{N(0, 100)}), under both dynamics. Because this prior
#' does not involve \eqn{(\mu, \rho)}, the joint full conditional of the
#' drift/intercept and the AR(1) coefficient is conjugate Gaussian, and both
#' are sampled by exact Gibbs draws -- \eqn{\rho} from its marginal truncated
#' to the stationary region \eqn{(-1, 1)}, then \eqn{\mu} given \eqn{\rho} --
#' with no Metropolis-Hastings correction.
#'
#' @section The four innovation structures:
#' The `innovations` argument controls the distribution of the latent
#' increments \eqn{\varepsilon_t = z_t - \mu - \rho z_{t-1}} (with
#' \eqn{\rho = 1} for the random walk and \eqn{\mu = 0} unless a drift/intercept
#' is included):
#' \describe{
#'   \item{`"gaussian"`}{\eqn{\varepsilon_t \sim N(0, \sigma^2)} with a single,
#'         constant variance.}
#'   \item{`"t"`}{A Student-t scale mixture: \eqn{\varepsilon_t \sim t_\nu(0,
#'         \sigma^2)}, robust to occasional large jumps. The degrees of
#'         freedom \eqn{\nu} are estimated from the data.}
#'   \item{`"mixture"`}{A finite scale mixture of normals with `K` components,
#'         resulting in a flexible increment distribution.}
#'   \item{`"sv"`}{Stochastic volatility: \eqn{\log\sigma^2_t} follows an AR(1)
#'         process (delegated to the \pkg{stochvol} package).}
#' }
#'
#' @section Zero inflation:
#' Set `zeros = "inflated"` (or `zero_inflation = TRUE`) to fit time-constant
#' zero inflation with either family: the observed count is
#' \eqn{y_t = v_t \tilde y_t}, where the gate
#' \eqn{v_t \sim \mathrm{Bernoulli}(\pi_{\mathrm{open}})} decides whether an
#' observed zero is \emph{structural} (gate closed) or a \emph{sampling} zero
#' produced by the Poisson/binomial process (gate open). The gate-open
#' probability \eqn{\pi_{\mathrm{open}}} is a single parameter that does not
#' vary over time or with covariates. See [structural_zero_prob()].
#'
#' Every fit stores two flavours of in-sample fitted values and replicates:
#' \emph{unconditional} draws (`fitted`, `yrep`), which include the gate and
#' are the quantities to compare with observed data (posterior predictive
#' checks), and \emph{conditional-on-gate-open} draws (`fitted_open`,
#' `yrep_open`), which come straight from the observation model and describe
#' the latent intensity process. Without zero inflation the pairs coincide
#' exactly. Response forecasts (`forecast_y`) are always unconditional -- the
#' gate is applied. See [predict.dynamic_fit()].
#'
#' @section Forecasting:
#' Forecasts are produced \emph{during} sampling: fit with `horizon = H` and `H`
#' latent states are appended to the state vector and sampled as missing values
#' inside the MCMC (no likelihood contribution, drawn from their GMRF full
#' conditionals). A response is drawn from the observation model at each sampled 
#' forecast state. Retrieve the stored draws with [forecast.dynamic_fit()].
#'
#' @references
#' Zens, G. and Bijak, J. (2026). Dynamic Count Models with Flexible Innovation
#' Processes for Irregular Maritime Migration. \emph{The Annals of Applied
#' Statistics}. \doi{10.1214/26-AOAS2171}.
#'
#' @aliases DynCount DynCount-package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats dbinom dnorm dpois plogis qlogis quantile rbeta
#'   rgamma rnorm rpois runif sd rbinom dgamma dexp
#' @importFrom graphics abline lines points polygon legend barplot
#' @importFrom grDevices adjustcolor
#' @importFrom utils head
## usethis namespace: end
NULL
