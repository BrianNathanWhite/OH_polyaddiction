# PARALLEL RUN OF OH_polyaddiction_model_spatial.R (4 drug classes x 3 chains) #
# each worker builds its own compiled model and runs one chain; launch from the repo root #
# needs up to 12 cores; per-chain seeds match the sequential script so results are identical #
# BRIAN N. WHITE #
# 2026-08-11 #

# LOAD R PACKAGES

library(parallel) # makeCluster, parLapply

# SPECIFY WORKER JOBS

  # one job per drug class/chain combination
  jobs <- expand.grid(class_name = c('cocaine_all', 'cocaine_no_fent',
                                     'psychostim_all', 'psychostim_no_fent'),
                      chain      = 1:3,
                      stringsAsFactors = F)

  n_workers <- min(nrow(jobs), detectCores() - 1)

# RUN CHAINS IN PARALLEL

  cluster <- makeCluster(n_workers)

  clusterExport(cluster, 'jobs')

  start <- Sys.time()

  chain_samples <- parLapply(cluster, seq_len(nrow(jobs)), function(job) {

    library(tidyverse)
    library(nimble)

    # source the model script build-only: defines code_spatial, constants, inits,
    # monitors, classes, make_array and the MCMC settings without running the fits;
    # local = T keeps BUILD_ONLY visible to the sourced guard
    BUILD_ONLY <- T
    source('code/OH_polyaddiction_model_spatial.R', local = T)

    class_name <- jobs$class_name[job]
    chain      <- jobs$chain[job]

    data_class <- list(O = make_array(classes[[class_name]]['O']),
                       E = make_array(classes[[class_name]]['E']))

    set.seed(1234)

    nimbleMCMC(code      = code_spatial,
               data      = data_class,
               constants = constants,
               inits     = inits,
               monitors  = monitors,
               niter     = n_iter,
               nburnin   = n_burnin,
               thin      = n_thin,
               nchains   = 1,
               setSeed   = chain,
               progressBar = F)

  })

  stopCluster(cluster)

  print(Sys.time() - start)

# COMBINE CHAINS AND SAVE ONE FILE PER DRUG CLASS

  dir.create('output/mcmc', recursive = T, showWarnings = F)

  for(class_name in unique(jobs$class_name)) {

    samples <- chain_samples[jobs$class_name == class_name]
    names(samples) <- paste0('chain', jobs$chain[jobs$class_name == class_name])

    save(samples, file = paste0('output/mcmc/mcmc_spatial_', class_name, '.Rda'))

  }
