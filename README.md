# OH_polyaddiction

Analysis code and data for "The impact of fentanyl on state- and county-level psychostimulant and cocaine overdose death rates by race in Ohio from 2010 to 2020: a time series and spatiotemporal analysis" (Estadt et al., *Harm Reduction Journal*, 2024, https://doi.org/10.1186/s12954-024-00936-9): Bayesian time series and spatiotemporal models of standardized mortality ratios (SMRs) for unintentional overdose deaths involving cocaine and psychostimulants, with and without fentanyl involvement, for Black and White populations across Ohio's 88 counties.

## Code

Run in order:

1. `OH_polyaddiction_data.R` builds the analysis data set: county overdose death counts for the four drug classes (all cocaine, cocaine not involving fentanyl, all psychostimulant, psychostimulant not involving fentanyl), expected counts via indirect age standardization to the 2010 Ohio rates, and the county adjacency structure.
2. `OH_polyaddiction_model_timeseries.R` fits the state-level time series SMR model (section 1.1 of the statistical supplement) for each drug class and saves the posterior samples; roughly half an hour per drug class.
3. `OH_polyaddiction_model_spatial.R` fits the county-level spatio-temporal SMR model (section 1.2 of the statistical supplement) for each drug class; roughly an hour per chain per class run sequentially. `OH_polyaddiction_model_spatial_parallel.R` spreads the twelve class/chain combinations over multiple cores and produces identical results.
4. `OH_polyaddiction_output_timeseries.R` checks convergence and produces Figures 1-3 and Table 1 of the manuscript plus the crude state rate table reported in the results.
5. `OH_polyaddiction_output_spatial.R` checks convergence and produces Figures 4-5 of the manuscript and Figures S1-S7 and Table S1 of Supplement 2.

Both model scripts run 1,000,000 MCMC iterations per chain (500,000 burn-in, thinned every 100th) across 3 chains with fixed seeds.

## Data

- `ODRace*_2023_05_18.csv` county-level unintentional drug overdose death counts by race and year for 2010-2020, extracted 2023-05-18 from the Ohio Public Health Information Warehouse (since migrated to https://data.ohio.gov); this is the extract behind the published results. Cocaine deaths are ICD-10 code T40.5, psychostimulant deaths are ICD-10 code T43.6, and fentanyl involvement is a positive mention of select text strings on the death certificate.
- `ODRaceFentanyl_07-19.csv` fentanyl-involved death counts from the earlier 2021-08-16 extract; used only for the crude fentanyl context rates and ends in 2019.
- `baseline_cocaine_no_fent.csv`, `baseline_psychostim_no_fent.csv` 2010 statewide overdose death counts not involving fentanyl by NCHS age group; these provide the standard rates for the indirect age standardization.
- `year_county_race_age_2023_05_18.csv` population denominators by year, county, race and NCHS age group.
- `shape_county_OH.Rda` Ohio county boundaries for mapping (US Census 2020 cartographic boundary file, pulled via tigris); the model adjacency structure is instead built from the maps-package county polygons, matching the published analysis.
- `data_for_analysis.Rda` the analysis data set produced by `OH_polyaddiction_data.R`.
- `metadata/` warehouse documentation: the extraction procedure, mortality data dictionary and warehouse data guide.

## Output

- `output/mcmc/` posterior samples, one file per model and drug class (not tracked; the spatial files run to a few hundred MB each)
- `output/diagnostics/` trace plots with Rhat and effective sample sizes
- `output/figures/` manuscript figures (`fig1`-`fig5`) and supplement figures (`figS1`-`figS7`)
- `output/tables/` Table 1, Table S1 and the crude state rates as csv, with docx versions of the tables
