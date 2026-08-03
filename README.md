# NNM1 — Neonatal mortality in excess of structure, Brazilian municipalities (2014–2024)

Analysis code for the study *"Beyond the North-South divide: neonatal mortality in
excess of structure across Brazilian municipalities (2014–2024)"*.

For each Brazilian municipality, we model the neonatal deaths expected from its
socioeconomic structure (and then access to care), read the residual excess as a
standardised mortality ratio (SMR) with empirical-Bayes shrinkage, decompose the
excess by avoidable cause, and quantify the avoidable burden across benchmarks.
The approach extends the residual-excess method of Martinelli et al. (2025,
*Public Health*) from adult mortality in a single metropolis to the neonatal
period across an entire country.

## Repository contents

```
notebooks/
  01_cohort_indicators.ipynb   neonatal cohort + avoidable-cause groups + municipal indicators
  02_excess_mortality.ipynb    nested negative binomial model (socioeconomic → +access) → SMR
                               with Byar CIs + empirical-Bayes shrinkage; cause decomposition;
                               avoidable burden ladder; robustness and data-quality checks
outputs/
  figures/                     final figures (excess map, burden map)
  tables/                      derived tables (municipal SMRs, avoidable ladder, coefficients)
requirements.txt               Python dependencies (Python 3.9)
```

Raw data are **not** versioned here — they are large and publicly available (see below).

## Data sources (all public)

- **SIM** (mortality) and **SINASC** (live births) — Brazilian Ministry of Health, [DATASUS](https://datasus.saude.gov.br/transferencia-de-arquivos/)
- **CNES** (neonatal ICU beds, primary care units) — DATASUS
- **IBGE/SIDRA** — municipal resident population
- **Atlas Brasil** — municipal Human Development Index (IDHM, 2010)
- Avoidable-cause list — Malta et al. (2007) / DATASUS

## Reproducing

```
pip install -r requirements.txt
```

Download the SIM, SINASC and CNES files for 2014–2024 from DATASUS into `data/raw/`,
then run the notebooks in order (`01` → `02`). The notebooks are self-documenting;
paths are relative to the project root.

## Data availability of results

`outputs/tables/excess_mortality.csv` contains the full municipal results
(observed and expected deaths, raw and empirical-Bayes SMR, confidence intervals).
