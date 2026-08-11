# CHECK CONVERGENCE AND PRODUCE FIGURES 4-5 OF THE MANUSCRIPT #
# PLUS FIGURES S1-S7 AND TABLE S1 OF SUPPLEMENT 2 #
# BRIAN N. WHITE #
# 2026-08-11 #

# LOAD R PACKAGES

library(tidyverse) # data manipulation and visualization
library(sf)        # county shape file handling
library(MCMCvis)   # MCMCtrace
library(cowplot)   # plot_grid
library(flextable) # make pretty tables

# LOAD ANALYSIS DATA

  load('data/data_for_analysis.Rda')
  load('data/shape_county_OH.Rda')

  C <- length(counties)
  T <- length(years)

  class_names <- c('cocaine_all', 'cocaine_no_fent', 'psychostim_all', 'psychostim_no_fent')

  # function to load the saved posterior samples for one model/drug class
  load_samples <- function(model, class_name) {

    samples_env <- new.env()
    load(paste0('output/mcmc/mcmc_', model, '_', class_name, '.Rda'), envir = samples_env)
    samples_env$samples

  }

  # function to pool the chains of a nimbleMCMC samples list into one matrix
  pool_chains <- function(samples) do.call(rbind, samples)

  # column name helpers matching nimble's node naming
  u_cols      <- function() as.vector(outer(1:C, 1:T, function(c, t) paste0('U[', c, ', ', t, ']')))
  lambda_cols <- function(r) as.vector(outer(1:C, 1:T, function(c, t) paste0('lambda[', c, ', ', t, ', ', r, ']')))

  # county/year index frame in the same column order as the helpers above
  county_year <- expand_grid(year = years, county = counties) %>%
    arrange(year, county)

# EXAMINE MCMC CONVERGENCE AND SUMMARIZE FIXED EFFECTS (TABLE S1)

  dir.create('output/diagnostics', recursive = T, showWarnings = F)
  dir.create('output/tables', recursive = T, showWarnings = F)

  table_s1 <- list()

  for(class_name in class_names) {

    samples <- load_samples('spatial', class_name)

    # note: delta[1] is fixed to 1 so its trace is a flat line by construction
    MCMCtrace(samples,
              params   = c('gamma0', 'gamma1', 'delta', 'eta', 'theta', 'tau_u', 'tau_eps'),
              open_pdf = F,
              Rhat     = T,
              n.eff    = T,
              ind      = T,
              filename = paste0('trace_spatial_', class_name, '.pdf'),
              wd       = 'output/diagnostics')

    # annual-change rate ratios from the race-specific slopes (1 = White, 2 = Black)
    gamma1 <- pool_chains(samples)[, c('gamma1[1]', 'gamma1[2]')]

    table_s1[[class_name]] <- tibble(White            = exp(gamma1[, 1]),
                                     Black            = exp(gamma1[, 2]),
                                     'Black vs White' = exp(gamma1[, 2] - gamma1[, 1])) %>%
      pivot_longer(everything(), names_to = 'quantity', values_to = 'draw') %>%
      group_by(quantity) %>%
      summarise(rr    = median(draw),
                lwr95 = quantile(draw, 0.025),
                upr95 = quantile(draw, 0.975),
                .groups = 'drop') %>%
      mutate(class = class_name, .before = 1)

    rm(samples)

  }

  table_s1 <- bind_rows(table_s1) %>%
    mutate(across(c(rr, lwr95, upr95), ~round(.x, 2)))

  write.csv(table_s1, 'output/tables/tableS1_spatial_rate_ratios.csv', row.names = F)

  table_s1 %>%
    flextable() %>%
    set_caption('Table S1: Spatio-temporal cocaine and psychostimulant-involved overdose rate ratio by annual change and race') %>%
    save_as_docx(path = 'output/tables/tableS1_spatial_rate_ratios.docx')

# MAP HELPERS

  dir.create('output/figures', recursive = T, showWarnings = F)

  # function to attach county/year values to the shape file for mapping
  map_data <- function(values) {

    shape_county_OH %>%
      left_join(values, by = c('NAME' = 'county'))

  }

  # function to make a year-faceted diverging choropleth; fill is on the log or
  # probability scale with legend breaks/labels supplied by the caller
  choropleth <- function(values, fill_lab, breaks, labels, midpoint = 0, limits = NULL) {

    ggplot(map_data(values)) +
      geom_sf(aes(fill = value), color = 'grey35', linewidth = 0.05) +
      facet_wrap(~year, ncol = 6) +
      scale_fill_gradient2(low      = '#2166AC',
                           mid      = '#F7F7F7',
                           high     = '#B2182B',
                           midpoint = midpoint,
                           limits   = limits,
                           breaks   = breaks,
                           labels   = labels,
                           name     = fill_lab) +
      theme_void() +
      theme(strip.text  = element_text(size = 9, face = 'bold'),
            legend.title = element_text(size = 9))

  }

# CREATE FIGURES 4-5: SHARED SPATIAL COMPONENT MAPS (FOR MANUSCRIPT)
# CREATE FIGURES S2-S7: RATE RATIO, SMR AND POSTERIOR PROBABILITY MAPS (FOR SUPPLEMENT)

  # substance label -> drug class used for the manuscript and supplement maps
  map_classes <- c(cocaine = 'cocaine_all', psychostimulant = 'psychostim_all')

  # accumulators for the processed posterior estimates written after the loop
  shared_component_estimates <- list()
  rate_ratio_estimates       <- list()
  smr_estimates              <- list()

  for(substance in names(map_classes)) {

    samples <- pool_chains(load_samples('spatial', map_classes[substance]))

    fig_num  <- ifelse(substance == 'cocaine', 4, 5)
    s_offset <- ifelse(substance == 'cocaine', 0, 1) # S2/S4/S6 cocaine, S3/S5/S7 psychostimulant

    # figure 4/5: posterior median of the shared spatio-temporal component
    u_median <- county_year %>%
      mutate(value = apply(samples[, u_cols()], 2, median))

    ggsave(plot = choropleth(u_median, 'Shared\ncomponent', breaks = waiver(), labels = waiver()),
           filename = paste0('fig', fig_num, '_shared_component_', substance, '.png'),
           path = 'output/figures', width = 12, height = 5, dpi = 'retina', bg = 'white')

    # figure S6/S7: posterior probability that the shared component exceeds zero
    u_prob <- county_year %>%
      mutate(value = colMeans(samples[, u_cols()] > 0))

    shared_component_estimates[[substance]] <- county_year %>%
      mutate(substance     = substance,
             u_median      = u_median$value,
             prob_positive = u_prob$value)

    ggsave(plot = choropleth(u_prob, 'P(U > 0)', breaks = seq(0, 1, 0.25), labels = waiver(),
                             midpoint = 0.5, limits = c(0, 1)),
           filename = paste0('figS', 6 + s_offset, '_shared_component_prob_', substance, '.png'),
           path = 'output/figures', width = 12, height = 5, dpi = 'retina', bg = 'white')

    # figure S2/S3: White/Black rate ratio of the SMRs and P(rate ratio > 1)
    log_rr <- log(samples[, lambda_cols(1)]) - log(samples[, lambda_cols(2)])

    rr_median <- county_year %>%
      mutate(value = apply(log_rr, 2, median))

    rr_prob <- county_year %>%
      mutate(value = colMeans(log_rr > 0))

    rate_ratio_estimates[[substance]] <- county_year %>%
      mutate(substance     = substance,
             log_rr_median = rr_median$value,
             prob_rr_gt1   = rr_prob$value)

    rr_breaks <- log(c(0.14, 0.22, 0.37, 0.61, 1, 1.65))

    fig_s23 <- plot_grid(choropleth(rr_median, 'White/Black\nrate ratio',
                                    breaks = rr_breaks, labels = round(exp(rr_breaks), 2)),
                         choropleth(rr_prob, 'P(rate\nratio > 1)',
                                    breaks = seq(0, 1, 0.25), labels = waiver(),
                                    midpoint = 0.5, limits = c(0, 1)),
                         ncol = 2)

    ggsave(plot = fig_s23,
           filename = paste0('figS', 2 + s_offset, '_rate_ratio_', substance, '.png'),
           path = 'output/figures', width = 20, height = 5, dpi = 'retina', bg = 'white')

    # figure S4/S5: log SMR by county, year and race (1 = White, 2 = Black)
    smr_breaks <- -1:4

    smr_median <- map(1:2, function(r) county_year %>%
                        mutate(value = apply(log(samples[, lambda_cols(r)]), 2, median)))

    smr_estimates[[substance]] <- map2_dfr(smr_median, c('White', 'Black'),
                                           ~mutate(.x, race = .y)) %>%
      mutate(substance = substance) %>%
      rename(log_smr_median = value)

    smr_map <- function(r, race_label) {

      choropleth(smr_median[[r]], paste0('SMR\n(', race_label, ')'),
                 breaks = smr_breaks, labels = round(exp(smr_breaks), 2),
                 limits = range(c(smr_breaks, smr_median[[r]]$value)))

    }

    fig_s45 <- plot_grid(smr_map(1, 'White'), smr_map(2, 'Black'), ncol = 2)

    ggsave(plot = fig_s45,
           filename = paste0('figS', 4 + s_offset, '_smr_', substance, '.png'),
           path = 'output/figures', width = 20, height = 5, dpi = 'retina', bg = 'white')

    rm(samples, log_rr)

  }

# WRITE POSTERIOR ESTIMATES (PROCESSED MCMC OUTPUT FOR REUSE WITHOUT REFITTING)

  dir.create('output/estimates', recursive = T, showWarnings = F)

  bind_rows(shared_component_estimates) %>%
    write.csv('output/estimates/spatial_shared_component.csv', row.names = F)

  bind_rows(rate_ratio_estimates) %>%
    write.csv('output/estimates/spatial_rate_ratio_white_black.csv', row.names = F)

  bind_rows(smr_estimates) %>%
    write.csv('output/estimates/spatial_smr.csv', row.names = F)

# CREATE FIGURE S1: COUNTY REFERENCE MAP WITH PERTINENT CITIES (FOR SUPPLEMENT)

  # city center coordinates
  # source: US Geographic Names Information System
  cities <- tibble(city = c('Cleveland', 'Akron', 'Columbus', 'Dayton', 'Cincinnati', 'Portsmouth'),
                   lon  = c(-81.6944, -81.5190, -82.9988, -84.1916, -84.5120, -82.9977),
                   lat  = c(41.4993, 41.0814, 39.9612, 39.7589, 39.1031, 38.7317)) %>%
    st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>%
    st_transform(st_crs(shape_county_OH))

  fig_s1 <- ggplot() +
    geom_sf(data = shape_county_OH, fill = 'white', color = 'grey40', linewidth = 0.2) +
    geom_sf(data = cities, size = 1.5, color = 'navy') +
    geom_sf_text(data = cities, aes(label = city),
                 nudge_x = 0.08, nudge_y = 0.09, hjust = 0, fontface = 'bold', size = 4) +
    theme_void()

  ggsave(plot = fig_s1, filename = 'figS1_county_reference_map.png', path = 'output/figures',
         width = 7, height = 6.5, dpi = 'retina', bg = 'white')
