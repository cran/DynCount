## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width = 7,
  fig.height = 4.2
)
set.seed(1)
library(DynCount)
# Short MCMC runs keep the vignette fast to build; use longer runs in practice.
NSAVE <- 1000L
NBURN <- 1000L

## ----simulate-----------------------------------------------------------------
sim <- simulate_dynamic_poisson(n = 80, sigma = 0.18, log_rate0 = 2.5, seed = 1)
str(sim, max.level = 1)
plot(sim$y, type = "h", xlab = "time", ylab = "count",
     main = "Simulated Poisson random walk")
lines(sim$rate, col = "steelblue", lwd = 2)

## ----fit-poisson--------------------------------------------------------------
fit <- fit_dynamic_model(sim$y, family = "poisson",
                         nsave = NSAVE, nburn = NBURN, seed = 1)
fit
summary(fit)

## ----plot-fitted--------------------------------------------------------------
plot_fitted(fit)

## ----plot-latent--------------------------------------------------------------
plot_latent(fit)

## ----forecast-----------------------------------------------------------------
fit_fc <- fit_dynamic_model(sim$y, family = "poisson", horizon = 8,
                            nsave = NSAVE, nburn = NBURN, seed = 1)
fc <- forecast(fit_fc)
fc                                 # prints the forecast path
fc$final                           # the single 8-step-ahead forecast
plot_forecast(fit_fc)

## ----ar1----------------------------------------------------------------------
# a genuinely stationary AR(1) log-rate: stationary mean 4, so mu = 4 * (1 - rho)
sim_ar <- simulate_dynamic_poisson(150, sigma = 0.2, log_rate0 = 4,
                                   rho = 0.9, mu = 0.4, seed = 3)
# no need to set include_mu: AR(1) enables the intercept automatically
fit_ar <- fit_dynamic_model(sim_ar$y, latent_dynamics = "ar1",
                            nsave = NSAVE, nburn = NBURN, seed = 3)
summary(fit_ar)                    # reports the posteriors of ar1_rho and mu

## ----offset-------------------------------------------------------------------
expo  <- log(runif(120, 50, 200))  # known exposure, e.g. population at risk
sim_o <- simulate_dynamic_poisson(120, sigma = 0.12, log_rate0 = -3.5,
                                  offset = expo, seed = 4)
fit_o <- fit_dynamic_model(sim_o$y, offset = expo, horizon = 6,
                           forecast_offset = log(120),
                           nsave = NSAVE, nburn = NBURN, seed = 4)
forecast(fit_o)$final

## ----med-t--------------------------------------------------------------------
data(med_weekly)
med <- tail(med_weekly$count, 120)
fit_med <- fit_dynamic_model(med, family = "poisson",
                             innovations = "t",
                             nsave = NSAVE, nburn = NBURN, seed = 2)
summary(fit_med)

## ----fit-zip------------------------------------------------------------------
data(uk_weekly)
uk <- uk_weekly$count[1:130]
mean(uk == 0)                         # many zeros
fit_zip <- fit_dynamic_model(uk, family = "poisson",
                             zero_inflation = TRUE,
                             nsave = NSAVE, nburn = NBURN, seed = 3)
summary(fit_zip)

## ----structural---------------------------------------------------------------
sz <- structural_zero_prob(fit_zip, zeros_only = TRUE)
head(sz, 10)
plot_zero_inflation(fit_zip)

## ----ppc-zip------------------------------------------------------------------
# zero proportion in the data vs both flavours of replicate
c(observed      = mean(uk == 0),
  yrep          = mean(fit_zip$draws$yrep == 0),        # gate applied: comparable
  yrep_open     = mean(fit_zip$draws$yrep_open == 0))   # gate-open only: too few zeros

## ----binomial-----------------------------------------------------------------
simb <- simulate_dynamic_binomial(n = 80, sigma = 0.12, trials = 50, seed = 4)
fit_bin <- fit_dynamic_model(simb$y, family = "binomial", trials = simb$trials,
                             horizon = 8, forecast_trials = 50,
                             nsave = NSAVE, nburn = NBURN, seed = 4)
summary(fit_bin)
plot_fitted(fit_bin)

## ----binomial-forecast--------------------------------------------------------
fc_bin <- forecast(fit_bin)
fc_bin$summary

## ----binomial-zip-------------------------------------------------------------
simz <- simulate_dynamic_binomial(n = 80, sigma = 0.1, trials = 40, logit0 = 1.5,
                                  zero_inflation = 0.2, seed = 7)
fit_bz <- fit_dynamic_model(simz$y, family = "binomial", trials = 40,
                            zero_inflation = TRUE,
                            nsave = NSAVE, nburn = NBURN, seed = 7)
head(structural_zero_prob(fit_bz))

## ----priors-------------------------------------------------------------------
dynamic_prior()
# a tighter prior favouring smoother latent paths
pr <- dynamic_prior(var_shape = 10, var_rate = 0.2)
fit_smooth <- fit_dynamic_model(sim$y, prior = pr,
                                nsave = NSAVE, nburn = NBURN, seed = 1)

