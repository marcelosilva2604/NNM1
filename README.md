# Posterior exceedance classification of residual area-level risk under sparse counts

Code, aggregated municipal data, and replicate-level outputs for a simulation study of the rule P(RR > tau) >= c in a Poisson model with a BYM2 municipal effect and restricted spatial regression, with an application to neonatal mortality in 5,565 Brazilian municipalities, 2014-2024.

## Layout

```
primary_fit/            the fit to the registry data (generating values for the simulation)
  src/04_fit_models.py  model definitions and fitting (spec bym2_rsr_full is the primary model)
  src/nc_io.py          helpers to save and load the posterior
  data/                 aggregated municipal table (analysis_muni.csv), adjacency pairs, centroids, scaling factor, adjacency report
  model/                hyperparameter summary and municipal estimates of the primary fit
simulation/
  src/11_simulation_grid.py      one replicate (--k, --rep) or aggregation (--aggregate)
  src/12_sensitivity_W.py        refits with alternative neighborhood graphs
  src/13_sensW_diagnostics.py    graph and divergence diagnostics for those refits
  src/14_make_supplement_tables.py, src/15_make_figure1.py
  jobs.txt, run_grid.sh, run.sh  resumable grid driver (4 fits at a time)
  outputs/grid/                  one JSON per replicate (k = 0 has 40, the other scenarios 100)
  outputs/grid_null_initial/     the initial null replicates run before the plan amendment
  outputs/grid_summary.csv, grid_calibration.csv, grid_replicate_meta.csv, supp_tables.md
  outputs/sens_W/                neighborhood-graph refits (fit_*.json, muni_*.csv, diag/)
ANALYSIS_PLAN.md       plan written before the simulation was run, with its single amendment
```

## Reproducing

1. Install the versions in `requirements.txt` (Python 3.12).
2. Fit the primary model: `python primary_fit/src/04_fit_models.py bym2_rsr_full` (4 chains, 4000 tuning and 4000 retained draws; about 3 minutes on 16 cores). The posterior file is about 340 MB and is not in this repository. A refit may not return the same draws as the one used for the reported results; the replicate JSON files record the index of the posterior draw used (`draw`).
3. Run the grid: `simulation/run_grid.sh` (resumable), then `python simulation/src/11_simulation_grid.py --aggregate`.
4. Tables and figure: `python simulation/src/14_make_supplement_tables.py` and `python simulation/src/15_make_figure1.py`.

The scripts read the primary fit from `primary_fit/` by default; set `NNM1_ROOT` to use another location.

## Not included

- The scripts that build the municipal table from individual-level registry extracts, because they read individual-level files. The resulting aggregated table is provided.
- The posterior file of the primary fit (size).

## Notes on the outputs

- The refits with alternative graphs in `simulation/outputs/sens_W/` (first run, script 12) and `simulation/outputs/sens_W/diag/` (diagnostic refits, script 13) did not give identical numbers for the same nominal settings (for example, 2,220 and 172 divergent transitions for the eight-neighbor graph at target acceptance 0.90); the reason is not known. The diagnostic refits are the ones summarized in the paper's supplement.
- R-hat of the hyperparameters and of the municipal effects was recorded separately only for replicates run after the plan amendment.

## Data sources

Public records of the Brazilian Ministry of Health (Mortality Information System and Live Birth Information System), the Brazilian Institute of Geography and Statistics, and the Atlas of Human Development in Brazil. The tables contain municipal aggregates only.

## Use of large language model tools

During the preparation of this work the authors used large language model tools to assist with writing and debugging analysis code and with language editing of the text. After using these tools, the authors reviewed and edited the content as needed and take full responsibility for the content of the published article.
