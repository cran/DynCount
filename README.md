# DynCount

**Bayesian dynamic models for count time series.**

`DynCount` fits Bayesian state-space models to count time series. A latent
trajectory `z[t]` follows flexible dynamics — a first-order **random walk** or
a stationary **AR(1)** process — and the observations are linked to it through
a **Poisson** (log link), a **binomial** (logit link) or a **multinomial**
(additive-log-ratio link, one latent series per non-baseline category)
observation model. Estimation is by Metropolis-within-Gibbs MCMC using the
Gaussian Markov random field full conditionals, and forecasts for any horizon
are obtained by forward simulation from the posterior draws.

The package implements and extends the methodology of Zens and Bijak (2026),
*Dynamic Count Models with Flexible Innovation Processes for Irregular
Maritime Migration*, *The Annals of Applied Statistics*,
[doi:10.1214/26-AOAS2171](https://doi.org/10.1214/26-AOAS2171); see
`citation("DynCount")`.

## Installation

```r
install.packages("DynCount")
```

## Quick start

```r
library(DynCount)

## simulate a Poisson random walk and recover the latent rate
sim <- simulate_dynamic_poisson(n = 80, sigma = 0.18, log_rate0 = 2.5, seed = 1)
fit <- fit_dynamic_model(sim$y, family = "poisson", seed = 1)  # latent_dynamics = "rw"

summary(fit)
plot_fitted(fit)

## forecast 20 steps ahead from the posterior draws
fc <- forecast(fit, horizon = 20, seed = 1)
fc$summary      # full path, one row per horizon
fc$final        # the single 20-step-ahead forecast
plot_forecast(fit, horizon = 20, seed = 1)
```

### AR(1) dynamics

AR(1) always carries an intercept. `include_mu` is enabled automatically,
which gives the process a non-zero stationary mean `mu / (1 - rho)`.

```r
# stationary AR(1) log-rate with mean 4  (mu = 4 * (1 - rho))
sim_ar <- simulate_dynamic_poisson(200, sigma = 0.2, log_rate0 = 4,
                                   rho = 0.9, mu = 0.4, seed = 3)
# include_mu is switched on automatically for AR(1)
fit_ar <- fit_dynamic_model(sim_ar$y, latent_dynamics = "ar1", seed = 3)
summary(fit_ar)            # posteriors of ar1_rho (in (-1, 1)) and intercept_mu
```

### Drift and offset

```r
## random-walk drift
sim_d <- simulate_dynamic_poisson(200, sigma = 0.12, log_rate0 = 1, mu = 0.03, seed = 4)
fit_d <- fit_dynamic_model(sim_d$y, include_mu = TRUE, seed = 4)   # rho = 1, mu sampled

## Poisson offset (log-exposure); supply forecast_offset for the future
expo  <- log(runif(200, 50, 200))
sim_o <- simulate_dynamic_poisson(200, sigma = 0.12, log_rate0 = -3.5,
                                  offset = expo, seed = 5)
fit_o <- fit_dynamic_model(sim_o$y, offset = expo, seed = 5)
forecast(fit_o, horizon = 8, forecast_offset = log(120))$final
```

### Heavy-tailed and time-varying innovations

```r
fit_t   <- fit_dynamic_model(med_weekly$count, innovations = "t", seed = 1)
fit_mix <- fit_dynamic_model(med_weekly$count, innovations = "mixture", seed = 1)
fit_sv  <- fit_dynamic_model(med_weekly$count, innovations = "sv", seed = 1)  # needs stochvol
summary(fit_t)             # t_df: degrees of freedom of the increments
```

### Zero inflation

```r
fit_zip <- fit_dynamic_model(uk_weekly$count, zero_inflation = TRUE, seed = 1)
summary(fit_zip)                 # gate_open_prob = P(gate open)
structural_zero_prob(fit_zip)    # P(structural) for each observed zero
plot_zero_inflation(fit_zip)
```

### Binomial model

```r
sim_b <- simulate_dynamic_binomial(n = 80, sigma = 0.12, trials = 50, seed = 1)
fit_b <- fit_dynamic_model(sim_b$y, family = "binomial", trials = sim_b$trials, seed = 1)
forecast(fit_b, horizon = 8, forecast_trials = 50)$final
```

### Multinomial choice counts

```r
# three categories, counts by week; C is the baseline of the simulation
sim_m <- simulate_dynamic_multinomial(n = 100, sigma = c(0.15, 0.1), trials = 300,
                                      alr0 = c(-1, 0.5), categories = c("A", "B", "C"),
                                      seed = 1)
fit_m <- fit_dynamic_model(sim_m$y, family = "multinomial",
                           baseline = sim_m$baseline, seed = 1)
fit_m                              # shows the categories and the baseline
summary(fit_m)                     # one innov_sd[...] row per non-baseline category
predict(fit_m, type = "prob")$summary   # posterior category shares, long format
forecast(fit_m, horizon = 8, forecast_trials = 300)$final
plot_fitted(fit_m)                 # one panel per category
plot_latent(fit_m, category = "A") # latent log-ratio of A vs the baseline
```

Posterior draws of multinomial fits carry a trailing category dimension, e.g.
`fit_m$draws$fitted_prob[, , "A"]` are the draws of the share of category A.
Without `baseline`, the fit uses the category with the largest total count.
The model is not invariant to this choice, because the dynamics are placed on
the log-ratios relative to the baseline, so pick a large category whose share
is stable.

## Example data

* `uk_weekly` — weekly English Channel crossings, 2018-W01 to 2025-W11.
* `med_weekly` — a longer weekly Mediterranean-crossings series with larger
  counts and few zeros, 2015-W40 to 2025-W11.

Both are weekly aggregates of detected irregular maritime crossings as used in
Zens and Bijak (2026). Both data frames have the columns `week` (ISO week label), `count` and `date`
(the Monday of the week).

## Documentation

```r
vignette("DynCount-intro", package = "DynCount")
```

