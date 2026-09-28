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

## ----predict------------------------------------------------------------------
head(predict(fit)$summary)

## ----forecast-----------------------------------------------------------------
fc <- forecast(fit, horizon = 8, seed = 1)
fc                                 # prints the forecast path
fc$final                           # the single 8-step-ahead forecast
plot_forecast(fit, horizon = 8, seed = 1)

## ----ar1----------------------------------------------------------------------
# a genuinely stationary AR(1) log-rate with stationary mean 4, so mu = 4 * (1 - rho)
sim_ar <- simulate_dynamic_poisson(n = 150, sigma = 0.2, log_rate0 = 4,
                                   rho = 0.9, mu = 0.4, seed = 3)
# no need to set include_mu, because AR(1) enables the intercept automatically
fit_ar <- fit_dynamic_model(sim_ar$y, family = "poisson", latent_dynamics = "ar1",
                            nsave = NSAVE, nburn = NBURN, seed = 3)
summary(fit_ar)                    # reports the posteriors of ar1_rho and intercept_mu

## ----offset-------------------------------------------------------------------
expo  <- log(runif(120, 50, 200))  # known exposure, e.g. population at risk
sim_o <- simulate_dynamic_poisson(n = 120, sigma = 0.12, log_rate0 = -3.5,
                                  offset = expo, seed = 4)
fit_o <- fit_dynamic_model(sim_o$y, family = "poisson", offset = expo,
                           nsave = NSAVE, nburn = NBURN, seed = 4)
forecast(fit_o, horizon = 6, forecast_offset = log(120), seed = 4)$final

## ----data---------------------------------------------------------------------
str(uk_weekly)
plot(med_weekly$date, med_weekly$count, type = "h", xlab = "week",
     ylab = "crossings", main = "Weekly Mediterranean crossings")

## ----med-t--------------------------------------------------------------------
med <- tail(med_weekly$count, 120)
fit_med <- fit_dynamic_model(med, family = "poisson",
                             innovations = "t",
                             nsave = NSAVE, nburn = NBURN, seed = 2)
summary(fit_med)

## ----med-mixture--------------------------------------------------------------
fit_mix <- fit_dynamic_model(med, family = "poisson", innovations = "mixture",
                             nsave = NSAVE, nburn = NBURN, seed = 2)
summary(fit_mix)

## ----med-sv, eval = requireNamespace("stochvol", quietly = TRUE)--------------
fit_sv <- fit_dynamic_model(med, family = "poisson", innovations = "sv",
                            nsave = NSAVE, nburn = NBURN, seed = 2)
summary(fit_sv)
vol <- sqrt(apply(fit_sv$draws$sig2, 2, median))
plot(vol, type = "l", xlab = "week", ylab = "increment SD (posterior median)",
     main = "Time-varying innovation SD")

## ----ess, eval = requireNamespace("coda", quietly = TRUE)---------------------
ess <- function(x) round(unname(coda::effectiveSize(x)))
c(innov_sd = ess(sqrt(fit$draws$innov_var)),
  z_40     = ess(fit$draws$z[, 40]),
  t_df     = ess(fit_med$draws$nu))

## ----fit-zip------------------------------------------------------------------
uk <- uk_weekly$count[1:130]
mean(uk == 0)                         # many zeros
fit_zip <- fit_dynamic_model(uk, family = "poisson",
                             zero_inflation = TRUE,
                             nsave = NSAVE, nburn = NBURN, seed = 3)
summary(fit_zip)

## ----structural---------------------------------------------------------------
sz <- structural_zero_prob(fit_zip)
head(sz, 10)
plot_zero_inflation(fit_zip)

## ----ppc-zip------------------------------------------------------------------
# zero proportion in the data vs both flavours of replicate
c(observed      = mean(uk == 0),
  yrep          = mean(fit_zip$draws$yrep == 0),        # gate applied, so comparable
  yrep_open     = mean(fit_zip$draws$yrep_open == 0))   # gate open only, so too few zeros

## ----binomial-----------------------------------------------------------------
simb <- simulate_dynamic_binomial(n = 80, sigma = 0.12, trials = 50, seed = 4)
fit_bin <- fit_dynamic_model(simb$y, family = "binomial", trials = simb$trials,
                             nsave = NSAVE, nburn = NBURN, seed = 4)
summary(fit_bin)
plot_fitted(fit_bin)

## ----binomial-forecast--------------------------------------------------------
fc_bin <- forecast(fit_bin, horizon = 8, forecast_trials = 50, seed = 4)
fc_bin$summary

## ----binomial-zip-------------------------------------------------------------
simz <- simulate_dynamic_binomial(n = 80, sigma = 0.1, trials = 40, logit0 = 1.5,
                                  zero_inflation = 0.2, seed = 7)
fit_bz <- fit_dynamic_model(simz$y, family = "binomial", trials = 40,
                            zero_inflation = TRUE,
                            nsave = NSAVE, nburn = NBURN, seed = 7)
summary(fit_bz)$params
head(structural_zero_prob(fit_bz))

## ----multinomial--------------------------------------------------------------
sim_m <- simulate_dynamic_multinomial(n = 80, sigma = c(0.15, 0.08), trials = 250,
                                      alr0 = c(-1, 0.3), baseline = 3,
                                      categories = c("A", "B", "C"), seed = 5)
head(sim_m$y)
colSums(sim_m$y)
fit_m <- fit_dynamic_model(sim_m$y, family = "multinomial",
                           baseline = sim_m$baseline,
                           nsave = NSAVE, nburn = NBURN, seed = 5)
fit_m
summary(fit_m)$params

## ----multinomial-plots, fig.height = 6----------------------------------------
head(predict(fit_m, type = "prob")$summary)
forecast(fit_m, horizon = 6, forecast_trials = 250, seed = 5)$final
plot_fitted(fit_m)
plot_latent(fit_m, category = "A")

## ----priors-------------------------------------------------------------------
dynamic_prior()
# an informative prior favouring rougher latent paths
pr <- dynamic_prior(var_shape = 2.5, var_rate = 0.5)
fit_inf <- fit_dynamic_model(sim$y, family = "poisson", prior = pr,
                             nsave = NSAVE, nburn = NBURN, seed = 1)
rbind(default = summary(fit)$params["innov_sd", ],
      informative = summary(fit_inf)$params["innov_sd", ])

