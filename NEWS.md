# DynCount 0.2.0

## New features

* New observation family `family = "multinomial"` for choice counts: `y` is an
  `n x K` matrix, the row totals are the known trials, and each non-baseline
  category has its own latent additive-log-ratio (ALR) series
  `log(p_k / p_baseline)` following the chosen latent dynamics and innovation
  structure, with no parameters shared across categories. The baseline is
  selected with the new `baseline` argument (default `"largest"`). Per-category
  ALR offsets, forecasting (`forecast_trials` = future totals, by default the
  last non-zero row total), `summary()`, `predict()` (including
  `type = "prob"`), `forecast()` and the plot functions (one panel per
  category, `category` argument) all support the new family; posterior draws
  carry a trailing category dimension. With `K = 2` the model coincides with
  the binomial family.
* New simulator `simulate_dynamic_multinomial()`. Its `baseline` can be given
  as a column index or a category name, and its `offset` as a scalar, a
  per-category vector or a matrix.
* Forecasting is now post hoc: `forecast(fit, horizon = H)` forward-simulates
  the latent path from every stored posterior draw (with increments from the
  fitted innovation structure) and draws responses from the observation
  model. This targets the same posterior predictive distribution as the
  previous in-sampler forecasts, but the horizon, the forecast offset and the
  forecast trials can now be chosen after fitting. `forecast()` and
  `plot_forecast()` gain `horizon`, `forecast_offset`, `forecast_trials` and
  `seed` arguments. Fitting with `horizon = H` still stores an `H`-step
  forecast in the fit.
* The draws of the stochastic-volatility parameters (`sv_mu`, `sv_phi`,
  `sv_sigma`) are now stored and reported by `summary()`. The draws of the
  mixture weights and component variances (`mix_weight`, `mix_var`) and of
  the overall innovation scale (`scale`) are stored as well.
* The example data sets gain a `date` column (the Monday of each ISO week).

## Changes to defaults

* The default prior on the innovation variance is now the vague
  `InvGamma(0.01, 0.01)` instead of `InvGamma(2.5, 0.5)`, which put almost no
  mass below `sigma = 0.1` and could force rough latent paths on smooth series.
  Pass `dynamic_prior(var_shape = 2.5, var_rate = 0.5)` to recover the old
  behaviour. See `?dynamic_prior` for the trade-offs.
* The relative variances of the `"mixture"` components now have their own
  hyperparameters `mix_var_shape` / `mix_var_rate` (defaults `2.5` / `0.5`,
  i.e. unchanged) rather than reusing `var_shape` / `var_rate`.

## Breaking changes

* `draws$z` and `draws$sig2` now have one column per observation. The initial
  latent state is stored separately as `draws$z0`, and forecast states are
  only in `draws$forecast_z`.
* `draws$gate` and `draws$pi_open` are `NULL` for fits without zero
  inflation, and `draws$nu` is `NULL` unless `innovations = "t"`.
* `forecast()` is now the generic from the generics package (re-exported), so
  `forecast(fit)` also works when the forecast or fable packages are
  attached.
* `zero_inflation` in `fit_dynamic_model()` must be a single `TRUE`/`FALSE`;
  other values (e.g. a probability, as in the simulators) are an error, as is
  combining `zero_inflation = TRUE` with a different explicit `zeros`.

## Sampler

* The sampler works on a matrix of latent series; the Poisson and binomial
  families are the single-column case.
* The latent states are updated in a checkerboard order. As the GMRF
  precision is tridiagonal, the states with odd and with even time index are
  conditionally independent given each other, so each half is updated in one
  vectorised step with the same single-site adaptive Metropolis moves. This
  makes fitting roughly 7 to 17 times faster.
* The Student-t degrees of freedom are updated with the mixing weights
  integrated out (a partially collapsed Gibbs step), which improves their
  effective sample size several-fold.

## Minor improvements and fixes

* A `seed` argument no longer changes the global random number stream. The
  previous RNG state is restored after fitting, forecasting and simulating.
* Arguments passed through `...` to the plot functions (e.g. `main`, `xlab`,
  `ylim`) now override the defaults instead of causing an error.
* Warnings for arguments that are ignored (`trials` or `forecast_trials` for
  the Poisson family, `baseline` for non-multinomial families, forecast
  inputs with `horizon = 0`) and for forecasting a model with an offset
  without a `forecast_offset`.
* `simulate_dynamic_binomial()` validates `trials`.

# DynCount 0.1.0

* Initial release.
