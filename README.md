# DynCount

**Bayesian dynamic models for Poisson and binomial time series.**

`DynCount` fits state-space models to non-Gaussian time series. A latent
trajectory `z[t]` follows flexible dynamics — a first-order **random walk** or
a stationary **AR(1)** process — and the observations are linked to it through
either a **Poisson** (log link) or a **binomial** (logit link) observation
model. Estimation is by Metropolis-within-Gibbs MCMC using the Gaussian Markov
random field full conditionals.

The package implements and extends the methodology of Zens and Bijak (2026),
*Dynamic Count Models with Flexible Innovation Processes for Irregular
Maritime Migration*, *The Annals of Applied Statistics*,
[doi:10.1214/26-AOAS2171](https://doi.org/10.1214/26-AOAS2171); see
`citation("DynCount")`.

## Installation

```r
install.packages("DynCount")
```

## Features

* **Latent dynamics**: first-order random walk (`"rw"`, default) or
  stationary AR(1) (`"ar1"`, with an intercept; `rho` sampled on `(-1, 1)`).
* Four **innovation structures** for the latent increments: `"gaussian"`,
  `"t"` (degrees of freedom sampled), `"mixture"` (finite scale mixture
  of normals), `"sv"` (stochastic volatility, via **stochvol**).
* **Forecasting during MCMC.** Set `horizon = H` when fitting and forecasts are
  produced inside the sampler, propagating parameter, state and innovation
  uncertainty. `forecast()` then extracts the stored draws.
* User-modifiable priors via `dynamic_prior()`; simulation, summaries, and
  plotting tools.
* **Optional drift / intercept** `mu` (`include_mu = TRUE`): a drift under the
  random walk and an intercept under AR(1). Disabled by default (`mu = 0`)
  under the random walk.
* **Poisson offset** (`offset =`): a known log-exposure term, so the mean is
  `exp(offset_t + z_t)`.
* Poisson (log link) and binomial (logit link, known trials) observation
  models, sharing one interface.
* **Time-constant zero inflation** for both families, separating structural
  from sampling zeros. Fitted values, replicates and forecasts are stored
  both unconditionally and conditionally on the gate being open; see
  `?predict.dynamic_fit`.


## Quick start

```r
library(DynCount)

## simulate a Poisson random walk and recover the latent rate
sim <- simulate_dynamic_poisson(n = 80, sigma = 0.18, log_rate0 = 2.5, seed = 1)
fit <- fit_dynamic_model(sim$y, family = "poisson", seed = 1)  # latent_dynamics = "rw"

summary(fit)
plot_fitted(fit)

## forecast 20 steps ahead: request the horizon when fitting
fit_fc <- fit_dynamic_model(sim$y, family = "poisson", horizon = 20, seed = 1)
fc <- forecast(fit_fc)
fc$summary      # full path, one row per horizon
fc$final        # the single 20-step-ahead forecast
plot_forecast(fit_fc)
```

### AR(1) dynamics

AR(1) always carries an intercept: `include_mu` is enabled automatically,
giving the process a non-zero stationary mean `mu / (1 - rho)`.

```r
# stationary AR(1) log-rate with mean 4  (mu = 4 * (1 - rho))
sim_ar <- simulate_dynamic_poisson(200, sigma = 0.2, log_rate0 = 4,
                                   rho = 0.9, mu = 0.4, seed = 3)
# include_mu is switched on automatically for AR(1)
fit_ar <- fit_dynamic_model(sim_ar$y, latent_dynamics = "ar1", seed = 3)
summary(fit_ar)            # posteriors of ar1_rho (in (-1, 1)) and intercept_mu
```

### Drift / intercept and offset

```r
## random-walk drift
sim_d <- simulate_dynamic_poisson(200, sigma = 0.12, log_rate0 = 1, mu = 0.03, seed = 4)
fit_d <- fit_dynamic_model(sim_d$y, include_mu = TRUE, seed = 4)   # rho = 1, mu sampled

## Poisson offset (log-exposure); supply forecast_offset for the future
expo <- log(runif(200, 50, 200))
fit_o <- fit_dynamic_model(sim_ar$y, offset = expo, horizon = 8,
                           forecast_offset = log(120), seed = 3)
forecast(fit_o)$final
```

### Zero inflation

```r
data(uk_weekly)
fit_zip <- fit_dynamic_model(uk_weekly$count, zero_inflation = TRUE, seed = 1)
structural_zero_prob(fit_zip)
plot_zero_inflation(fit_zip)
```

### Binomial model

```r
sim_b <- simulate_dynamic_binomial(n = 80, sigma = 0.12, trials = 50, seed = 1)
fit_b <- fit_dynamic_model(sim_b$y, family = "binomial", trials = sim_b$trials,
                           horizon = 8, forecast_trials = 50)
forecast(fit_b)$final
```

## Example data

* `uk_weekly` — weekly English Channel crossings.
* `med_weekly` — a longer weekly Mediterranean-crossings series with larger
  counts and few zeros.

Both are weekly aggregates of detected irregular maritime crossings as used in
Zens and Bijak (2026).

## Documentation

```r
vignette("DynCount-intro", package = "DynCount")
```

