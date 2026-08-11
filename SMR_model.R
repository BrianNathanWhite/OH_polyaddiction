#load required packages
library(tidyverse)
library(reshape2)
library(nimble)
library(sf)
library(maps)
library(maptools) # map2SpatialPolygon
library(sp) # maptools dependency
library(spdep)

setwd("C:/Users/bnwhite/OneDrive - Wake Forest Baptist Health/projects/Kline")
# mortality counts for drug w/o fentanyl for 2010 in Ohio
baseline_coca <- read.csv("./polyaddiction/Ohio Stimulant OD/Brian/data/baseline_cocaine_no_fent.csv") # cocaine
baseline_stim <- read.csv("./polyaddiction/Ohio Stimulant OD/Brian/data/baseline_psychostim_no_fent.csv")  # psychostimulants

# population counts for year/county/race/age combinations in Ohio
year_county_race_age <- read.csv('./polyaddiction/Ohio Stimulant OD/Brian/data/year_county_race_age_2023_05_18.csv')

# mortality counts for drug overdose (any instance of drug and only drug/fentanyl combo)
od_coca <- read.csv('./polyaddiction/Ohio Stimulant OD/Brian/data/ODRaceCocaine_2023_05_18.csv')
od_fent_coca <- read.csv('./polyaddiction/Ohio Stimulant OD/Brian/data/ODRaceFentanylCocaine_2023_05_18.csv') 
od_stim <- read.csv('./polyaddiction/Ohio Stimulant OD/Brian/data/ODRacePsychostimulant_2023_05_18.csv')
od_fent_stim <- read.csv('./polyaddiction/Ohio Stimulant OD/Brian/data/ODRaceFentanylPsychostimulant_2023_05_18.csv')

# read shape files into sf files for use in R
shape_county <- st_read("./data/Shapefiles/2020 County/cb_2020_us_county_500k.shp", stringsAsFactors = F)

# filter out data for OH
shape_county_OH <- shape_county %>% filter(STUSPS == 'OH')


# rename columns
names(baseline_coca) <- c('year', 'age', 'c_only_deaths', 'population')
names(baseline_stim) <- c('year', 'age', 's_only_deaths', 'population')

# remove SORT column
baseline_coca <- baseline_coca[, -5]
baseline_stim <- baseline_stim[, -5]

# combine baseline drug data sets, reorder columns, filter out 'Total' years and 'Unk' age
baseline_drug <- baseline_coca %>%
                    left_join(baseline_stim, by = c('year', 'age', 'population')) %>%
                    select(year, age, c_only_deaths, s_only_deaths, population) %>%
                    filter(year != 'Total') %>%
                    filter(age!= 'Unk') %>%
                    select(-year)

# rename columns
names(od_coca) <- c('race','county','year','c_deaths','population')
names(od_fent_coca) <- c('race','county','year','fc_deaths','population')
names(od_stim) <- c('race','county','year','s_deaths','population')
names(od_fent_stim) <- c('race','county','year','fs_deaths','population')

# remove SORT column
od_coca <- od_coca[, -6]
od_fent_coca <- od_fent_coca[, -6]
od_stim <- od_stim[, -6]
od_fent_stim <- od_fent_stim[, -6]

# combine drug data sets, reorder columns, filter for black/white race, filter out 'NonOH' and 'Unknown' counties, filter years from 2010-2019, make sure each variable has the right class
od_drug <- od_coca %>%
              left_join(od_fent_coca, by = c('race','year','county','population')) %>%
              left_join(od_stim, by = c('race','year','county','population')) %>%
              left_join(od_fent_stim, by = c('race','year','county','population')) %>%
              select(race, county, year, c_deaths, fc_deaths, s_deaths, fs_deaths, population) %>%
              filter(year >= 2010 & year <= 2020) %>%
              filter(race %in% c('Black', 'White')) %>%
              filter(!county %in% c('NonOH', 'Unknown')) %>%
              mutate_at(vars('race', 'county'), as.factor) %>%
              mutate(year = as.integer(year))


# rename columns, filter years from 2010-2020, remove rows with unknown or total entries (these all have 0 mortality counts)
year_county_race_age <- year_county_race_age %>%
                            rename(race = RacePopRaceDesc,
                                   age = AgeGroupNCHSAgeNCHS,
                                   year = PopulationYearYear,
                                   county = CountyPopCountyName,
                                   population = Count) %>%
                            select(-SORT) %>%
                            filter(!age %in% c('Unk', 'Unknown')) %>%
                            filter(population != 'Unknown' & county!= 'Unknown') %>%
                            filter(race != 'Total' &
                                     age != 'Total' &
                                     year != 'Total' & 
                                     !county %in% c('Total', 'NonOH')
                            ) %>%
                            filter(year >= 2010 & year <= 2020) %>%
                            mutate_at(vars('race', 'age', 'county'), as.factor) %>%
                            mutate(year = as.integer(year)) %>%
                            arrange(county, year, race)

# compute overdose death rate per 100k
baseline_drug <- baseline_drug %>%
                    mutate(c_only_death_rate = c_only_deaths/population,
                           c_only_death_rate_per100k = c_only_death_rate*10^5,
                           s_only_death_rate = s_only_deaths/population,
                           s_only_death_rate_per100k = s_only_death_rate*10^5
                    )

# number of rows excluding age
n <- nrow(year_county_race_age)/nrow(baseline_drug)

# compute the expected death count for each age/county/year/race
year_county_race_age <- year_county_race_age %>%
                            mutate(s_only_death_rate = rep(baseline_drug$s_only_death_rate, n), 
                                   c_only_death_rate = rep(baseline_drug$c_only_death_rate, n)) %>%
                            mutate(E_s = population*s_only_death_rate,
                                   E_c = population*c_only_death_rate)

# sum the expected death counts over the age categories
death_count <- year_county_race_age %>%
                    group_by(race, year, county) %>%
                    summarise(E_s = sum(E_s), E_c = sum(E_c))

# subtract fentanyl laced drug overdose deaths from larger group to obtain cocaine-only overdose deaths and psychostim only overdose deaths, respectively
od_drug <- od_drug %>%
              mutate(c_only_deaths = c_deaths - fc_deaths,
                     s_only_deaths = s_deaths - fs_deaths) %>%
              select(race, county, year, c_deaths, fc_deaths, c_only_deaths, s_deaths, fs_deaths, s_only_deaths, population)

# observed death count for 4 (non-disjoint) groups: x, fx where x = cocaine, psychostimulants. Here, x refers to drug-x-only whereas fx is drug-x-only + fentanyl-laced-drug-x.
death_count <- death_count %>%
                  left_join(od_drug, by = c('race', 'year', 'county')) %>%
                  select(race, 
                         year, 
                         county, 
                         E_s, 
                         s_deaths, 
                         s_only_deaths, 
                         E_c, 
                         c_deaths, 
                         c_only_deaths, 
                         population) %>%
                  rename(O_s = s_only_deaths, O_fs = s_deaths, O_c = c_only_deaths, O_fc = c_deaths)

### STATE-LEVEL TIME-SERIES ANALYSIS ###

# sum over counties to get expected/observed counts for time-series model
death_count_ts <- death_count %>%
                    group_by(race, year) %>%
                    summarise(E_s = sum(E_s),
                              O_s = sum(O_s),
                              E_fs = sum(E_s), # this is the same as E_s as fentanyl wasn't around in 2010
                              O_fs = sum(O_fs),
                              E_c = sum(E_c),
                              O_c = sum(O_c),
                              E_fc = sum(E_c), # this is the same as E_c for the same reason as above
                              O_fc = sum(O_fc)
                    ) %>%
                    mutate(year_int = year - 2010)

# set 'White' to be reference level for main effect (race)
death_count_ts$race <- relevel(death_count_ts$race, 'White')

# set seed for reproducibility of nimble output
set.seed(1234)

# construct design matrix
X <- model.matrix(~ race*year_int, data = death_count_ts) # use year_int = year - 2010
p <- ncol(X)

# specify nimble model
code <- nimbleCode({
  
  
  for(i in 1:N) {
    
    O[i] ~ dpois(E[i]*lambda[i]) # observed mortality counts
    
    log(lambda[i]) <- inprod(X[i, 1:p], beta[1:p]) + b_star[i]
    
  }
  
  for(i in 1:p) {
    
    beta[i] ~ dflat() # uninformative improper prior for fixed effects
    
  }
  
  # enforce sum to 0 constraint over time within race group                   
  b_star[1] <- b[1] - mean(b[1:(N/2)])
  b_star[N/2 + 1] <- b[N/2 + 1] - mean(b[(N/2 + 1):N])
  
  for(i in 2:(N/2)){
    
    b_star[i] <- b[i] - mean(b[1:(N/2)]) # black
    b_star[N/2 + i] <- b[N/2 + i] - mean(b[(N/2 + 1):N]) # white
    
  }
  
  # random effect for race with AR(1) temporal structure
  b[1] ~ dnorm(0,sigma1)
  b[N/2+1] ~ dnorm(0, sigma2)
  
  for (i in 2:(N/2)){ 
    
    b[i] ~ dnorm(eta1*b[i-1], sigma1) # black
    b[N/2+i] ~ dnorm(eta2*b[N/2+i-1], sigma2) # white
    
  }
  
  # priors for random effect
  eta1 ~ dunif(-1,1)
  eta2 ~ dunif(-1,1)
  
  sigma1 ~ dgamma(0.1, 0.1)
  sigma2 ~ dgamma(0.1, 0.1)
  
  
})

# data
data_s <-  list(O = death_count_ts$O_s,
                E = death_count_ts$E_s)

data_fs <-  list(O = death_count_ts$O_fs,
                 E = death_count_ts$E_fs)

data_c <-  list(O = death_count_ts$O_c,
                E = death_count_ts$E_c)

data_fc <-  list(O = death_count_ts$O_fc,
                 E = death_count_ts$E_fc)

# constants 
constants <- list(N = nrow(death_count_ts),
                  X = X,
                  p = p)

# initial values of parameters
inits <- list(beta = rep(0, p), sigma1 = 0.5, sigma2 = 0.5, eta1 = 0.5, eta2 = 0.5)

# specify parameters to monitor
monitors <- c("lambda", "beta", "sigma1", "sigma2", "eta1", "eta2")

# specify MCMC details
n_iter <- 1000000
n_burnin <- 500000
n_chains <- 3
thin <- 100
n_samples <- (n_iter - n_burnin)/thin


# fit time-series model
mcmc_s <- nimbleMCMC(code = code,
                     data = data_s,
                     inits = inits,
                     constants = constants,
                     monitors = monitors,
                     niter = n_iter,
                     nburnin = n_burnin,
                     nchains = n_chains,
                     thin = thin,
                     progressBar = T)

mcmc_fs <- nimbleMCMC(code = code,
                      data = data_fs,
                      inits = inits,
                      constants = constants,
                      monitors = monitors,
                      niter = n_iter,
                      nburnin = n_burnin,
                      nchains = n_chains,
                      thin = thin,
                      progressBar = T)

mcmc_c <- nimbleMCMC(code = code,
                     data = data_c,
                     inits = inits,
                     constants = constants,
                     monitors = monitors,
                     niter = n_iter,
                     nburnin = n_burnin,
                     nchains = n_chains,
                     thin = thin,
                     progressBar = T)

mcmc_fc <- nimbleMCMC(code = code,
                      data = data_fc,
                      inits = inits,
                      constants = constants,
                      monitors = monitors,
                      niter = n_iter,
                      nburnin = n_burnin,
                      nchains = n_chains,
                      thin = thin,
                      progressBar = T)

#save(mcmc_c, file = "../output/modeling/timeseries/runs/mcmc_c.Rda")
#save(mcmc_fc, file = "../output/modeling/timeseries/runs/mcmc_fc.Rda")
#save(mcmc_s, file = "../output/modeling/timeseries/runs/mcmc_s.Rda")
#save(mcmc_fs, file = "../output/modeling/timeseries/runs/mcmc_fs.Rda")
#load("../output/modeling/timeseries/runs/mcmc_c.Rda")
#load("../output/modeling/timeseries/runs/mcmc_fc.Rda")
#load("../output/modeling/timeseries/runs/mcmc_s.Rda")
#load("../output/modeling/timeseries/runs/mcmc_fs.Rda")

### COUNTY-LEVEL SPATIO-TEMPORAL ANALYSIS ###

# spatial data for ICAR prior (i.e. incidince matrix, spatial weights matrix)
OH.county <- map("county", "Ohio", plot = FALSE, fill = TRUE)
OH.map <- map2SpatialPolygons(OH.county, IDs = OH.county$names)
W1 <- nb2WB(poly2nb(OH.map))

# observed and expected overdose death counts
death_count <- death_count %>%
                  mutate(year_int = year - 2010, # shift year column to begin in year 0 and end in year 9 
                         E_fc = E_c, # add expected death counts for fc and fs. These are same as c and s as they are for 2010
                         E_fs = E_s) %>%
                  select(year_int, county, race, O_c, O_fc, O_s, O_fs, E_c, E_fc, E_s, E_fs)

# split data by race
death_count_b <- death_count %>% filter(race == 'Black')
death_count_w <- death_count %>% filter(race == 'White')

# set seed to make nimble model output reproducible
set.seed(1234)

code_spatial <- nimbleCode({
  
  #first year
  for (t in 2010:2010){
    
    for (i in 1:88){
      
      #data model
      Ow[88*(t-2010)+i] ~ dpois(Ew[88*(t-2010)+i]*lambda_w[88*(t-2010)+i])
      Ob[88*(t-2010)+i] ~ dpois(Eb[88*(t-2010)+i]*lambda_b[88*(t-2010)+i])
      
      #GLM
      log(lambda_w[88*(t-2010)+i]) <- alpha0 + alpha1*(t-2010) + U[88*(t-2010)+i] + Vw[88*(t-2010)+i]
      log(lambda_b[88*(t-2010)+i]) <- beta0 + beta1*(t-2010) + delta*U[88*(t-2010)+i] + Vb[88*(t-2010)+i]
      
      #race specific error
      Vw[88*(t-2010)+i] ~ dnorm(0, tau_Vw)
      Vb[88*(t-2010)+i] ~ dnorm(0, tau_Vb)
      
      #shared component
      U[88*(t-2010)+i] <- uu[88*(t-2010)+i]
      
    }
    
    #shared component ICAR
    uu[(88*(t-2010)+1):(88*(t-2010)+88)] ~ dcar_normal(adj[], weights[], num[], tau_U,zero_mean=1)
    
  }  
  
  #other years
  for (t in 2011:2020){
    
    for (i in 1:88){
      
      #data model
      Ow[88*(t-2010)+i] ~ dpois(Ew[88*(t-2010)+i]*lambda_w[88*(t-2010)+i])
      Ob[88*(t-2010)+i] ~ dpois(Eb[88*(t-2010)+i]*lambda_b[88*(t-2010)+i])
      
      #GLM
      log(lambda_w[88*(t-2010)+i]) <- alpha0 + alpha1*(t-2010) +  U[88*(t-2010)+i] + Vw[88*(t-2010)+i]
      log(lambda_b[88*(t-2010)+i]) <- beta0 + beta1*(t-2010) + delta*U[88*(t-2010)+i] + Vb[88*(t-2010)+i]
      
      #race specific error
      Vw[88*(t-2010)+i] ~ dnorm(phi[1]*Vw[88*((t-2010)-1)+i], tau_Vw) 
      Vb[88*(t-2010)+i] ~ dnorm(phi[2]*Vb[88*((t-2010)-1)+i], tau_Vb) 
      
      #shared component
      U[88*(t-2010)+i] <- eta*U[88*((t-2010)-1)+i] + uu[88*(t-2010)+i]
      
    }
    
    #shared component ICAR
    uu[(88*(t-2010)+1):(88*(t-2010)+88)] ~ dcar_normal(adj[], weights[], num[], tau_U,zero_mean=1)
    
    
  }
  
  # priors
  
  # fixed effects
  alpha0 ~ dflat()
  alpha1 ~ dflat()
  beta0 ~ dflat()
  beta1 ~ dflat()
  
  # auto-regressive parameter for shared spatial component
  eta ~ dunif(0,1) 
  
  # auto-regressive parameters for race-specific variation
  phi[1] ~ dunif(0,1) # White
  phi[2] ~ dunif(0,1) # Black
  
  delta ~ dflat() # race specific loading for shared spatial component
  
  # variances
  tau_U ~ dgamma(0.5,0.5) # shared spatial component
  tau_Vw ~ dgamma(0.5,0.5) # White variation
  tau_Vb ~ dgamma(0.5,0.5) # Black variation
  
  
})

# data
data_c_spatial <-  list(Ow = death_count_w$O_c,
                        Ob = death_count_b$O_c,
                        Ew = death_count_w$E_c,
                        Eb = death_count_b$E_c)

data_fc_spatial <-  list(Ow = death_count_w$O_fc,
                         Ob = death_count_b$O_fc,
                         Ew = death_count_w$E_fc,
                         Eb = death_count_b$E_fc)

data_s_spatial <-  list(Ow = death_count_w$O_s,
                        Ob = death_count_b$O_s,
                        Ew = death_count_w$E_s,
                        Eb = death_count_b$E_s)

data_fs_spatial <-  list(Ow = death_count_w$O_fs,
                         Ob = death_count_b$O_fs,
                         Ew = death_count_w$E_fs,
                         Eb = death_count_b$E_fs)

# constants
constants_spatial <- list(adj = W1$adj,
                          weights = W1$weights,
                          num = W1$num)

# initial values for model parameters and hyper-parameters
inits_spatial <- list(alpha0 = 0,
                      alpha1 = 0,
                      beta0 = 0,
                      beta1 = 0,
                      Vw = rep(0, 88*11),
                      Vb = rep(0, 88*11),
                      uu = rep(0, 88*11),
                      tau_Vw = 1,
                      tau_Vb = 1,
                      tau_U = 1,
                      eta = 0.9,
                      phi = rep(0.9, 2),
                      delta = 1)

# specify parameters to monitor
monitors_spatial <- c("lambda_w",
                      "lambda_b",
                      "alpha0",
                      "alpha1",
                      "beta0",
                      "beta1",
                      "Vw",
                      "Vb",
                      "U",
                      "tau_Vw",
                      "tau_Vb",
                      'tau_U',
                      'eta',
                      'phi',
                      'delta')

# specify MCMC details
n_iter_spatial <- 10^6
n_burnin_spatial <- n_iter_spatial/2
n_chains_spatial <- 1
thin_spatial <- 100
n_samples_spatial <- (n_iter_spatial - n_burnin_spatial)/thin_spatial

# fit models
start <- Sys.time()

mcmc_c_spatial <- nimbleMCMC(code = code_spatial,
                             data = data_c_spatial,
                             inits = inits_spatial,
                             constants = constants_spatial,
                             monitors = monitors_spatial,
                             niter = n_iter_spatial,
                             nburnin = n_burnin_spatial,
                             nchains = n_chains_spatial,
                             thin = thin_spatial,
                             progressBar = T)

mcmc_fc_spatial <- nimbleMCMC(code = code_spatial,
                              data = data_fc_spatial,
                              inits = inits_spatial,
                              constants = constants_spatial,
                              monitors = monitors_spatial,
                              niter = n_iter_spatial,
                              nburnin = n_burnin_spatial,
                              nchains = n_chains_spatial,
                              thin = thin_spatial,
                              progressBar = T)

mcmc_s_spatial <- nimbleMCMC(code = code_spatial,
                             data = data_s_spatial,
                             inits = inits_spatial,
                             constants = constants_spatial,
                             monitors = monitors_spatial,
                             niter = n_iter_spatial,
                             nburnin = n_burnin_spatial,
                             nchains = n_chains_spatial,
                             thin = thin_spatial,
                             progressBar = T)

mcmc_fs_spatial <- nimbleMCMC(code = code_spatial,
                              data = data_fs_spatial,
                              inits = inits_spatial,
                              constants = constants_spatial,
                              monitors = monitors_spatial,
                              niter = n_iter_spatial,
                              nburnin = n_burnin_spatial,
                              nchains = n_chains_spatial,
                              thin = thin_spatial,
                              progressBar = T)

finish <- Sys.time()
finish - start

mcmc_spatial <- list(mcmc_c_spatial, mcmc_fc_spatial, mcmc_s_spatial, mcmc_fs_spatial)
#save(mcmc_spatial, file = "../output/OH_stimulant_output.Rda")