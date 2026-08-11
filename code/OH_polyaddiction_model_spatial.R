# FIT COUNTY-LEVEL SPATIO-TEMPORAL SMR MODEL FOR THE FOUR DRUG CLASSES #
# model specified in section 1.2 of the statistical supplement, adapted from #
# Kline, Pan and Hepler (2021, Epidemiology 32:295-302) #
# BRIAN N. WHITE #
# 2026-08-11 #

# LOAD R PACKAGES

library(tidyverse) # data manipulation
library(nimble)    # nimbleCode, nimbleMCMC

# LOAD ANALYSIS DATA

load('data/data_for_analysis.Rda')

# PREPARE MODEL INPUTS

  C <- length(counties)
  T <- length(years)

  # order as White block then Black block, counties alphabetical within year within
  # race, matching the county x year x race array construction below
  county_data <- county_data %>%
    mutate(race = factor(race, levels = c('White', 'Black'))) %>%
    arrange(race, year, county)

  # function to reshape a count column into a county x year x race array (1 = White, 2 = Black)
  make_array <- function(col) array(county_data[[col]], dim = c(C, T, 2))

# DEFINE NIMBLE MODEL

code_spatial <- nimbleCode({

  for(r in 1:2) {

    # observed counts and log-linear SMR model (supplement eq. 6 and 8)
    for(t in 1:T) {

      for(c in 1:C) {

        O[c, t, r] ~ dpois(E[c, t, r]*lambda[c, t, r])

        log(lambda[c, t, r]) <- gamma0[r] + gamma1[r]*(t - 1) + delta[r]*U[c, t] + eps[c, t, r]

      }

    }

    # race-specific residual variation with AR(1) temporal structure within county
    # (supplement eq. 10)
    for(c in 1:C) {

      eps[c, 1, r] ~ dnorm(0, tau = tau_eps[r])

      for(t in 2:T) {

        eps[c, t, r] ~ dnorm(theta[r]*eps[c, t-1, r], tau = tau_eps[r])

      }

    }

    # race-specific fixed effects: intercept and years since 2010
    gamma0[r] ~ dflat()
    gamma1[r] ~ dflat()

    theta[r] ~ dunif(0, 1)
    tau_eps[r] ~ dgamma(0.5, 0.5)

  }

  # shared spatio-temporal component: ICAR innovations propagated by a shared AR(1)
  # (supplement eq. 9); zero_mean enforces the yearly sum-to-zero constraint
  uu[1:C, 1] ~ dcar_normal(adj[1:L], weights[1:L], num[1:C], tau = tau_u, zero_mean = 1)

  for(c in 1:C) {

    U[c, 1] <- uu[c, 1]

  }

  for(t in 2:T) {

    uu[1:C, t] ~ dcar_normal(adj[1:L], weights[1:L], num[1:C], tau = tau_u, zero_mean = 1)

    for(c in 1:C) {

      U[c, t] <- eta*U[c, t-1] + uu[c, t]

    }

  }

  # race-specific loading on the shared component; White fixed to 1 for identifiability
  delta[1] <- 1
  delta[2] ~ dflat()

  eta ~ dunif(0, 1)
  tau_u ~ dgamma(0.5, 0.5)

})

# SPECIFY MCMC DETAILS

  n_iter   <- 1000000
  n_burnin <- 500000
  n_thin   <- 100
  n_chains <- 3

  constants <- list(C       = C,
                    T       = T,
                    L       = length(W$adj),
                    adj     = W$adj,
                    weights = W$weights,
                    num     = W$num)

  inits <- list(gamma0  = rep(0, 2),
                gamma1  = rep(0, 2),
                eps     = array(0, dim = c(C, T, 2)),
                uu      = matrix(0, C, T),
                delta   = c(NA, 1),
                eta     = 0.9,
                theta   = rep(0.9, 2),
                tau_u   = 1,
                tau_eps = rep(1, 2))

  monitors <- c('lambda', 'U', 'gamma0', 'gamma1', 'delta', 'eta', 'theta', 'tau_u', 'tau_eps')

  # each drug class pairs its observed counts with the matching expected counts;
  # the all-involvement classes share the no-fentanyl expected counts because the
  # 2010 standard predates fentanyl
  classes <- list(cocaine_all        = c(O = 'O_cocaine_all',        E = 'E_cocaine'),
                  cocaine_no_fent    = c(O = 'O_cocaine_no_fent',    E = 'E_cocaine'),
                  psychostim_all     = c(O = 'O_psychostim_all',     E = 'E_psychostim'),
                  psychostim_no_fent = c(O = 'O_psychostim_no_fent', E = 'E_psychostim'))

# FIT MODEL FOR EACH DRUG CLASS

  # guard so OH_polyaddiction_model_spatial_parallel.R can source this file
  # build-only (skips the fits); running this file directly fits the four classes
  # sequentially (roughly an hour per chain per class), use the parallel script
  # to spread the class/chain combinations over multiple cores
  if(!exists('BUILD_ONLY')) {

  dir.create('output/mcmc', recursive = T, showWarnings = F)

  for(class_name in names(classes)) {

    data_class <- list(O = make_array(classes[[class_name]]['O']),
                       E = make_array(classes[[class_name]]['E']))

    set.seed(1234)

    start <- Sys.time()

    samples <- nimbleMCMC(code      = code_spatial,
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

    save(samples, file = paste0('output/mcmc/mcmc_spatial_', class_name, '.Rda'))

  }

  }
