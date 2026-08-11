# FIT STATE-LEVEL TIME SERIES SMR MODEL FOR THE FOUR DRUG CLASSES #
# model specified in section 1.1 of the statistical supplement #
# BRIAN N. WHITE #
# 2026-08-11 #

# LOAD R PACKAGES

library(tidyverse) # data manipulation
library(nimble)    # nimbleCode, nimbleMCMC

# LOAD ANALYSIS DATA

load('data/data_for_analysis.Rda')

# PREPARE MODEL INPUTS

  # order as White block then Black block with years ascending within race
  state_data <- state_data %>%
    mutate(race = factor(race, levels = c('White', 'Black'))) %>%
    arrange(race, year)

  # design matrix: intercept, years since 2010, race (White reference) and interaction
  X <- model.matrix(~ race*year_int, data = state_data)
  p <- ncol(X)
  N <- nrow(state_data)
  T <- length(years)

# DEFINE NIMBLE MODEL

code_timeseries <- nimbleCode({

  # observed counts and log-linear SMR model (supplement eq. 1 and 4)
  for(i in 1:N) {

    O[i] ~ dpois(E[i]*lambda[i])

    log(lambda[i]) <- inprod(X[i, 1:p], beta[1:p]) + eps[race_idx[i], year_idx[i]]

  }

  # race-specific residual heterogeneity with AR(1) temporal structure (supplement eq. 5)
  for(r in 1:2) {

    eps[r, 1] ~ dnorm(0, tau = tau_eps[r])

    for(t in 2:T) {

      eps[r, t] ~ dnorm(phi[r]*eps[r, t-1], tau = tau_eps[r])

    }

  }

  # priors: flat on fixed effects, uniform(0,1) on AR parameters and
  # gamma(0.5, 0.5) on precisions, i.e. inverse gamma(0.5, 0.5) on variances
  for(j in 1:p) {

    beta[j] ~ dflat()

  }

  for(r in 1:2) {

    phi[r] ~ dunif(0, 1)
    tau_eps[r] ~ dgamma(0.5, 0.5)

  }

})

# SPECIFY MCMC DETAILS

  n_iter   <- 1000000
  n_burnin <- 500000
  n_thin   <- 100
  n_chains <- 3

  constants <- list(N        = N,
                    p        = p,
                    T        = T,
                    X        = X,
                    race_idx = as.integer(state_data$race), # 1 = White, 2 = Black
                    year_idx = state_data$year_int + 1)

  inits <- list(beta    = rep(0, p),
                eps     = matrix(0, 2, T),
                phi     = rep(0.5, 2),
                tau_eps = rep(1, 2))

  monitors <- c('beta', 'lambda', 'phi', 'tau_eps')

  # each drug class pairs its observed counts with the matching expected counts;
  # the all-involvement classes share the no-fentanyl expected counts because the
  # 2010 standard predates fentanyl
  classes <- list(cocaine_all        = c(O = 'O_cocaine_all',        E = 'E_cocaine'),
                  cocaine_no_fent    = c(O = 'O_cocaine_no_fent',    E = 'E_cocaine'),
                  psychostim_all     = c(O = 'O_psychostim_all',     E = 'E_psychostim'),
                  psychostim_no_fent = c(O = 'O_psychostim_no_fent', E = 'E_psychostim'))

# FIT MODEL FOR EACH DRUG CLASS

  dir.create('output/mcmc', recursive = T, showWarnings = F)

  for(class_name in names(classes)) {

    data_class <- list(O = state_data[[classes[[class_name]]['O']]],
                       E = state_data[[classes[[class_name]]['E']]])

    set.seed(1234)

    start <- Sys.time()

    samples <- nimbleMCMC(code      = code_timeseries,
                          data      = data_class,
                          constants = constants,
                          inits     = inits,
                          monitors  = monitors,
                          niter     = n_iter,
                          nburnin   = n_burnin,
                          thin      = n_thin,
                          nchains   = n_chains,
                          setSeed   = 1:n_chains,
                          progressBar = T)

    print(paste0(class_name, ' run time:'))
    print(Sys.time() - start)

    save(samples, file = paste0('output/mcmc/mcmc_timeseries_', class_name, '.Rda'))

  }
