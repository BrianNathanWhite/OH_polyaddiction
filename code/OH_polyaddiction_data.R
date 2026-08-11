# BUILD ANALYSIS DATA SET FOR OHIO STIMULANT OVERDOSE SMR MODELS #
# BRIAN N. WHITE #
# 2026-08-11 #

# LOAD R PACKAGES

library(tidyverse) # data manipulation
library(sf)        # county shape file handling
library(spdep)     # poly2nb, nb2WB
library(maps)      # county polygons for the adjacency structure
library(tigris)    # pull OH county shape file from US Census

# LOAD RAW WAREHOUSE EXTRACTS

  # county-level unintentional drug overdose death counts by race and year,
  # extracted 2023-05-18 (the pull used for the published analysis)
  # source: Ohio Public Health Information Warehouse (since migrated to https://data.ohio.gov)
  # extraction procedure: data/metadata/DataPull_081621.pdf
  # cocaine deaths are ICD-10 T40.5, psychostimulant deaths are ICD-10 T43.6, and
  # fentanyl involvement is a positive mention of select text strings on the death certificate
  od_cocaine_all     <- read.csv('data/ODRaceCocaine_2023_05_18.csv')                 # any cocaine involvement
  od_cocaine_fent    <- read.csv('data/ODRaceFentanylCocaine_2023_05_18.csv')         # cocaine and fentanyl jointly
  od_psychostim_all  <- read.csv('data/ODRacePsychostimulant_2023_05_18.csv')         # any psychostimulant involvement
  od_psychostim_fent <- read.csv('data/ODRaceFentanylPsychostimulant_2023_05_18.csv') # psychostimulant and fentanyl jointly

  # any fentanyl involvement, from the earlier 2021-08-16 pull; used only for the
  # crude context rates in the results text and ends in 2019
  od_fentanyl <- read.csv('data/ODRaceFentanyl_07-19.csv')

  # 2010 statewide overdose death counts not involving fentanyl by NCHS age group (all races)
  baseline_cocaine    <- read.csv('data/baseline_cocaine_no_fent.csv')
  baseline_psychostim <- read.csv('data/baseline_psychostim_no_fent.csv')

  # population counts by year, county, race and NCHS age group
  population <- read.csv('data/year_county_race_age_2023_05_18.csv')

# CLEAN WAREHOUSE EXTRACTS

  # function to standardize an overdose extract: rename columns, restrict to
  # Black/White residents of the 88 Ohio counties, drop total/unknown rows
  clean_od <- function(od, deaths_name) {

    od %>%
      rename(race       = DeathRaceRaceDesc,
             county     = DeathCountyCountyName,
             year       = DeathYearYear,
             deaths     = Deaths,
             population = Population) %>%
      select(-SORT) %>%
      filter(race %in% c('Black', 'White'),
             !county %in% c('NonOH', 'Unknown', 'Total'),
             !year %in% c('Total', 'Unknown')) %>%
      mutate(year = as.integer(year)) %>%
      rename(!!deaths_name := deaths)

  }

  # combine the four drug classes; the joins include population because every
  # extract carries the same county/race/year denominators
  od_drug <- clean_od(od_cocaine_all, 'cocaine_all') %>%
    left_join(clean_od(od_cocaine_fent, 'cocaine_fent'), by = c('race', 'county', 'year', 'population')) %>%
    left_join(clean_od(od_psychostim_all, 'psychostim_all'), by = c('race', 'county', 'year', 'population')) %>%
    left_join(clean_od(od_psychostim_fent, 'psychostim_fent'), by = c('race', 'county', 'year', 'population'))

  # analysis years: the study period is 2010 (the standardization year) through
  # 2020, matching the manuscript; the upper cap guards against future extracts
  # that extend past the study period
  years <- od_drug %>% filter(year >= 2010 & year <= 2020) %>% pull(year) %>% unique() %>% sort()

  # deaths not involving fentanyl are the class total minus the joint fentanyl counts
  od_drug <- od_drug %>%
    filter(year %in% years) %>%
    mutate(cocaine_no_fent    = cocaine_all - cocaine_fent,
           psychostim_no_fent = psychostim_all - psychostim_fent)

  # population denominators by age group for the same counties, races and years
  population <- population %>%
    rename(race       = RacePopRaceDesc,
           age        = AgeGroupNCHSAgeNCHS,
           year       = PopulationYearYear,
           county     = CountyPopCountyName,
           population = Count) %>%
    select(-SORT) %>%
    filter(race %in% c('Black', 'White'),
           !age %in% c('Unk', 'Unknown', 'Total'),
           !county %in% c('NonOH', 'Unknown', 'Total'),
           !year %in% c('Total', 'Unknown')) %>%
    mutate(year       = as.integer(year),
           population = as.numeric(population)) %>%
    filter(year %in% years)

# COMPUTE EXPECTED DEATH COUNTS BY INDIRECT AGE STANDARDIZATION

  # age-specific standard rates from the 2010 statewide baseline
  # note: the baseline denominators cover all races, not only Black/White
  clean_baseline <- function(baseline, rate_name) {

    baseline %>%
      rename(year       = DeathYearYear,
             age        = DeathAgeGroupNCHSAgeNCHS,
             deaths     = Deaths,
             population = Population) %>%
      filter(year == '2010', age != 'Unk') %>%
      transmute(age, !!rate_name := deaths/population)

  }

  baseline_rates <- clean_baseline(baseline_cocaine, 'cocaine_rate') %>%
    left_join(clean_baseline(baseline_psychostim, 'psychostim_rate'), by = 'age')

  # apply the 2010 age-specific rates to each county/year/race population and sum
  # over age groups; the standard is the no-fentanyl rate, and because fentanyl
  # involvement was negligible in 2010 (3 deaths among Black Ohioans statewide)
  # the same standard serves both the all-involvement and no-fentanyl classes
  expected <- population %>%
    left_join(baseline_rates, by = 'age') %>%
    group_by(race, county, year) %>%
    summarise(E_cocaine    = sum(population*cocaine_rate),
              E_psychostim = sum(population*psychostim_rate),
              .groups = 'drop')

  # check: every age group in the population data matched a standard rate
  stopifnot(!anyNA(expected))

# ASSEMBLE ANALYSIS DATA SET

  # county-level observed and expected counts for the four drug classes
  county_data <- expected %>%
    left_join(od_drug, by = c('race', 'county', 'year')) %>%
    transmute(race, county, year,
              year_int             = year - min(years),
              O_cocaine_all        = cocaine_all,
              O_cocaine_no_fent    = cocaine_no_fent,
              O_psychostim_all     = psychostim_all,
              O_psychostim_no_fent = psychostim_no_fent,
              E_cocaine,
              E_psychostim,
              population) %>%
    arrange(race, year, county)

  # check: complete panel of 88 counties x 2 races x all years with no missing values
  stopifnot(nrow(county_data) == 88*2*length(years), !anyNA(county_data))

  # state-level counts for the time series model and crude rates
  state_data <- county_data %>%
    group_by(race, year, year_int) %>%
    summarise(across(c(starts_with('O_'), starts_with('E_'), population), sum),
              .groups = 'drop')

  # state-level fentanyl-involved deaths for the crude context rates reported in the results
  state_fentanyl <- clean_od(od_fentanyl, 'fentanyl_deaths') %>%
    filter(year %in% years) %>%
    group_by(race, year) %>%
    summarise(fentanyl_deaths = sum(fentanyl_deaths),
              population      = sum(population),
              .groups = 'drop')

# COUNTY SHAPE FILE AND ADJACENCY STRUCTURE

  # un-comment these two lines if you need to re-pull the shape file
  # source: US Census Bureau 2020 cartographic boundary files via tigris
  #shape_county_OH <- counties(state = 'OH', cb = T, year = 2020) %>% arrange(NAME)
  #save(shape_county_OH, file = 'data/shape_county_OH.Rda')
  load('data/shape_county_OH.Rda')

  # check: shape file county names and order match the analysis data
  counties <- sort(unique(county_data$county))
  stopifnot(identical(shape_county_OH$NAME, counties))

  # queen-contiguity neighborhood structure in the format nimble's dcar_normal expects
  # approach: built from the maps-package county polygons, matching the adjacency used
  # for the published analysis (via the now-retired maptools); the Census cartographic
  # boundaries differ at the Auglaize/Darke/Mercer/Shelby corner point, so the Census
  # shape file above is used for mapping only
  county_polygons <- st_as_sf(map('county', 'Ohio', plot = F, fill = T)) %>%
    mutate(county = str_remove(ID, '^ohio,')) %>%
    arrange(county)

  # check: polygon order matches the analysis data
  stopifnot(identical(county_polygons$county, tolower(counties)))

  W <- nb2WB(poly2nb(county_polygons))

  # check: all counties have at least one neighbor
  stopifnot(all(W$num > 0))

# SAVE ANALYSIS DATA SET

  save(county_data, state_data, state_fentanyl, counties, years, W,
       file = 'data/data_for_analysis.Rda')
