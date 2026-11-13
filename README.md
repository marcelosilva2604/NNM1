# NNM1: maternal composition and area-level variation in neonatal mortality, Brazil, 2014-2024

Analysis code, aggregated municipal data and model outputs for the manuscript **"Maternal composition and area-level variation in
neonatal mortality in Brazil, 2014-2024: an ecological study with Bayesian spatial analysis"**.

The study uses registered neonatal deaths (0-27 days) and live births of 5,565 Brazilian municipalities, pooled over 2014-2024 and
analysed by municipality of maternal residence. It compares Poisson models with a restricted BYM2 spatial effect with and without
maternal composition, describes the residual variation and its persistence between 2014-2019 and 2020-2024, and evaluates by
simulation a classification rule based on P(RR > 1.10). It describes variation conditional on the fitted models; it does not
estimate causal effects.

## What is in this repository

| Path | Content |
|---|---|
| `v2/src/` | Analysis scripts, steps 01 to 14 (see the table below) |
| `v2/data/` | Municipal aggregates used by the models: composition of live births, deaths by municipality and year, adjacency pairs, analysis tables |
| `v2/outputs/` | Model summaries and results: `NUMBERS.json`, tables, simulation replicates, figures |
| `v2/supplement/` | Supplementary material of the manuscript (methods, specification, priors, diagnostics and tables) |
| `outputs/tables/excess_mortality.csv` | Municipal comparison columns read by step 03 (earlier municipal analysis, comparison only) |
| `DATA_DICTIONARY.md` | Columns of `v2/outputs/tables/municipal_estimates.csv` |

Posterior samples (`*.nc`), individual-level records and the IBGE meshes are **not** included.

## Data sources (public)

Neonatal deaths: Mortality Information System (SIM). Live births and maternal characteristics: Live Birth Information System (Sinasc).
Both are maintained by the Brazilian Ministry of Health (DATASUS). Neonatal intensive care beds and facility coordinates: National Registry
of Health Facilities (CNES). Resident population: Brazilian Institute of Geography and Statistics (IBGE). Municipal Human Development Index
(2010): Atlas of Human Development in Brazil. Municipal meshes: IBGE.

## Pipeline

| Step | Script | What it does | Needs files not distributed |
|---|---|---|---|
| 01 | `01_maternal_composition.py` | Aggregates maternal composition (age, schooling, race/skin colour, antenatal visits, gestation, birthweight, place of birth) by municipality and year from the Sinasc files. Ages outside 10-60 years are excluded from both age-category numerators and from the known-age denominator. | Sinasc files |
| 02 | `02_adjacency.py` | Queen-contiguity graph of the municipal mesh and the BYM2 scaling factor | IBGE mesh |
| 03 | `03_build_dataset.py` | Municipal analysis table (deaths, births, covariates, access indicators, health regions) | SIM, CNES, population, MHDI |
| 04 | `04_fit_models.py` | Fits the models (PyMC) | no |
| 05 | `05_results.py` | Residual relative risks, exceedance probabilities, count departures, comparisons, `NUMBERS.json` | no |
| 06 | `06_manuscript_tables.py` | Manuscript tables from the outputs | no |
| 07 | `07_period_persistence.py` | Period-specific fits (2014-2019 and 2020-2024) and persistence | SIM/Sinasc files |
| 08 | `08_simulation_calibration.py` | Simulation with known residual risk (40 datasets) and without residual (12 datasets) | no |
| 09 | `09_added_value.py` | Comparisons between models | no |
| 10 | `10_missingness_scenarios.py` | Missing-input scenarios for schooling and race/skin colour; refits the primary model | no |
| 11 | `11_compare_scenarios.py` | Posterior summaries of the scenarios against the primary fit | no |
| 12-13 | `12_figures_flow_dag.py`, `13_figure_map.py` | Figures 1 to 3 | IBGE mesh (map) |
| 14 | `14_evalues.py`, `14b_check_evalues_with_EValue_package.R` | E-values for the covariate rate ratios and their check against the R package EValue | no |

The aggregated files in `v2/data/` allow steps 04 to 11 and 14 to be rerun without the files listed above. Settings, seeds and diagnostics
are in the supplement. The model fits can be run with `NNM_DRAWS` and `NNM_TUNE` (retained draws and warm-up per chain) and `NNM1_ROOT`
(root of a working copy of `v2/`).

## Environment

`v2/requirements.txt` lists the versions used (PyMC 6.2.0, ArviZ 1.2.0, PyTensor 3.2.4, NumPy 2.4.6, pandas 3.0.5, SciPy 1.18.0,
xarray 2026.7.0, h5netcdf 1.8.1, GeoPandas 1.1.4, libpysal 4.15.0, Polars 1.44.2 for step 01). Posterior files are written with
`nc_io.py`, which falls back to NetCDF-3 when `h5py` is not installed.

## Limits

The analysis is ecological and describes municipal aggregates. Residual variation can reflect unmeasured risk, model
misspecification, registration error and sampling variation. The registration indicators describe recorded fields, not the completeness
of registration. See the manuscript and the supplement for the diagnostics and sensitivity analyses.
