# CHECK CONVERGENCE AND PRODUCE FIGURES 1-3 AND TABLE 1 OF THE MANUSCRIPT #
# PLUS THE CRUDE STATE RATE TABLE REPORTED IN THE RESULTS #
# BRIAN N. WHITE #
# 2026-08-11 #

# LOAD R PACKAGES

library(tidyverse) # data manipulation and visualization
library(MCMCvis)   # MCMCtrace
library(cowplot)   # plot_grid
library(flextable) # make pretty tables

# LOAD ANALYSIS DATA AND MCMC OUTPUT

  load('data/data_for_analysis.Rda')

  # order must match the model input ordering in OH_polyaddiction_model_timeseries.R
  state_data <- state_data %>%
    mutate(race = factor(race, levels = c('White', 'Black'))) %>%
    arrange(race, year)

  class_names <- c('cocaine_all', 'cocaine_no_fent', 'psychostim_all', 'psychostim_no_fent')

  # function to load the saved posterior samples for one model/drug class
  load_samples <- function(model, class_name) {

    samples_env <- new.env()
    load(paste0('output/mcmc/mcmc_', model, '_', class_name, '.Rda'), envir = samples_env)
    samples_env$samples

  }

  mcmc_timeseries <- map(set_names(class_names), ~load_samples('timeseries', .x))

  # function to pool the chains of a nimbleMCMC samples list into one matrix
  pool_chains <- function(samples) do.call(rbind, samples)

# EXAMINE MCMC CONVERGENCE

  dir.create('output/diagnostics', recursive = T, showWarnings = F)

  for(class_name in class_names) {

    MCMCtrace(mcmc_timeseries[[class_name]],
              params   = c('beta', 'phi', 'tau_eps'),
              open_pdf = F,
              Rhat     = T,
              n.eff    = T,
              ind      = T,
              filename = paste0('trace_timeseries_', class_name, '.pdf'),
              wd       = 'output/diagnostics')

  }

# CREATE FIGURE 1: STATEWIDE CRUDE OVERDOSE DEATH RATES PER 100K (FOR MANUSCRIPT)

  dir.create('output/figures', recursive = T, showWarnings = F)

  crude_rates <- state_data %>%
    transmute(race, year,
              'All cocaine'                            = O_cocaine_all/population*10^5,
              'Cocaine not involving fentanyl'         = O_cocaine_no_fent/population*10^5,
              'All psychostimulant'                    = O_psychostim_all/population*10^5,
              'Psychostimulant not involving fentanyl' = O_psychostim_no_fent/population*10^5) %>%
    pivot_longer(-c(race, year), names_to = 'class', values_to = 'rate') %>%
    mutate(race = factor(race, levels = c('Black', 'White')))

  # function to make one substance panel of figure 1
  fig1_panel <- function(class_levels) {

    crude_rates %>%
      filter(class %in% class_levels) %>%
      mutate(class = factor(class, levels = class_levels)) %>%
      ggplot(aes(x = year, y = rate, linetype = class)) +
      geom_line(linewidth = 0.4) +
      facet_wrap(~race) +
      scale_linetype_manual(values = c('solid', 'dotted')) +
      scale_x_continuous(breaks = seq(min(years), max(years), 2)) +
      labs(x = 'Year', y = 'Mortality Rate', linetype = NULL) +
      theme_bw() +
      theme(legend.position        = 'inside',
            legend.position.inside = c(0.02, 0.98),
            legend.justification   = c(0, 1),
            legend.background      = element_blank(),
            legend.key.height      = unit(0.35, 'cm'),
            panel.grid             = element_blank())

  }

  fig1 <- plot_grid(fig1_panel(c('All cocaine', 'Cocaine not involving fentanyl')),
                    fig1_panel(c('All psychostimulant', 'Psychostimulant not involving fentanyl')),
                    ncol   = 1,
                    labels = c('(A)', '(B)'))

  ggsave(plot = fig1, filename = 'fig1_state_rates.png', path = 'output/figures',
         width = 7, height = 8, dpi = 'retina', bg = 'white')

# CREATE FIGURES 2-3: SMR CREDIBLE INTERVALS OVER TIME BY RACE (FOR MANUSCRIPT)

  # function to summarize the posterior log SMR by race and year for one drug class
  smr_summary <- function(samples, class_label) {

    log_lambda <- log(pool_chains(samples)[, paste0('lambda[', 1:nrow(state_data), ']')])

    tibble(race    = state_data$race,
           year    = state_data$year,
           class   = class_label,
           log_smr = apply(log_lambda, 2, median),
           lwr95   = apply(log_lambda, 2, quantile, probs = 0.025),
           upr95   = apply(log_lambda, 2, quantile, probs = 0.975))

  }

  # function to make the SMR credible interval figure for one substance
  smr_figure <- function(class_all, class_no_fent, label_all, label_no_fent) {

    bind_rows(smr_summary(mcmc_timeseries[[class_all]], label_all),
              smr_summary(mcmc_timeseries[[class_no_fent]], label_no_fent)) %>%
      mutate(class = factor(class, levels = c(label_all, label_no_fent)),
             race  = factor(race, levels = c('Black', 'White'))) %>%
      ggplot(aes(x = year, y = log_smr, linetype = class)) +
      geom_ribbon(aes(ymin = lwr95, ymax = upr95, group = class), fill = 'grey60', alpha = 0.4) +
      geom_line(linewidth = 0.4) +
      facet_wrap(~race) +
      scale_linetype_manual(values = c('solid', 'dotted')) +
      scale_x_continuous(breaks = seq(min(years), max(years), 2)) +
      scale_y_continuous(sec.axis = sec_axis(~., name   = 'SMR',
                                                 breaks = -4:6,
                                                 labels = round(exp(-4:6), 2))) +
      labs(x = 'Year', y = 'log(SMR)', linetype = NULL) +
      theme_bw() +
      theme(legend.position        = 'inside',
            legend.position.inside = c(0.02, 0.98),
            legend.justification   = c(0, 1),
            legend.background      = element_blank(),
            legend.key.height      = unit(0.35, 'cm'),
            panel.grid             = element_blank())

  }

  ggsave(plot = smr_figure('cocaine_all', 'cocaine_no_fent',
                           'All cocaine', 'Cocaine not involving fentanyl'),
         filename = 'fig2_smr_cocaine.png', path = 'output/figures',
         width = 8, height = 4.5, dpi = 'retina', bg = 'white')

  ggsave(plot = smr_figure('psychostim_all', 'psychostim_no_fent',
                           'All psychostimulant', 'Psychostimulant not involving fentanyl'),
         filename = 'fig3_smr_psychostimulant.png', path = 'output/figures',
         width = 8, height = 4.5, dpi = 'retina', bg = 'white')

# CREATE TABLE 1: RATE RATIOS BY ANNUAL CHANGE AND RACE (FOR MANUSCRIPT)

  dir.create('output/tables', recursive = T, showWarnings = F)

  # function to summarize the annual-change rate ratios from the fixed effects
  # design matrix columns: intercept, race (Black), years since 2010, interaction
  rr_summary <- function(samples, class_label) {

    beta <- pool_chains(samples)[, paste0('beta[', 1:4, ']')]

    tibble(White            = exp(beta[, 3]),
           Black            = exp(beta[, 3] + beta[, 4]),
           'Black vs White' = exp(beta[, 4])) %>%
      pivot_longer(everything(), names_to = 'quantity', values_to = 'draw') %>%
      group_by(quantity) %>%
      summarise(rr    = median(draw),
                lwr95 = quantile(draw, 0.025),
                upr95 = quantile(draw, 0.975),
                .groups = 'drop') %>%
      mutate(class = class_label, .before = 1)

  }

  table1 <- imap_dfr(mcmc_timeseries, ~rr_summary(.x, .y)) %>%
    mutate(across(c(rr, lwr95, upr95), ~round(.x, 2)))

  write.csv(table1, 'output/tables/table1_rate_ratios.csv', row.names = F)

  table1 %>%
    flextable() %>%
    set_caption('Table 1: Cocaine- and psychostimulant-involved overdose rate ratios by annual change and race') %>%
    save_as_docx(path = 'output/tables/table1_rate_ratios.docx')

# CREATE CRUDE STATE RATE TABLE (FOR RESULTS TEXT)

  crude_table <- state_data %>%
    left_join(select(state_fentanyl, race, year, fentanyl_deaths), by = c('race', 'year')) %>%
    transmute(race, year,
              cocaine_all_rate        = round(O_cocaine_all/population*10^5, 2),
              cocaine_no_fent_rate    = round(O_cocaine_no_fent/population*10^5, 2),
              psychostim_all_rate     = round(O_psychostim_all/population*10^5, 2),
              psychostim_no_fent_rate = round(O_psychostim_no_fent/population*10^5, 2),
              fentanyl_rate           = round(fentanyl_deaths/population*10^5, 2))

  write.csv(crude_table, 'output/tables/crude_state_rates.csv', row.names = F)
