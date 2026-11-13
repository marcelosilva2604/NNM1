# Data dictionary

`v2/outputs/tables/municipal_estimates.csv`: one row per municipality (5,565 rows, 26 columns). CrI = credible interval.

| Column | Description |
|---|---|
| `CODMUNRES` | Six-digit IBGE code of the municipality of residence |
| `municipio` | Municipality name |
| `state` | State abbreviation (UF) |
| `REG` | Macro-region (North, Northeast, Central-West, Southeast, South) |
| `deaths` | Neonatal deaths (0 to 27 days), 2014-2024 pooled |
| `births` | Live births, 2014-2024 pooled |
| `nmr` | Crude neonatal mortality rate per 1,000 live births |
| `expected_median`, `expected_lo95`, `expected_hi95` | Expected deaths given municipal structure and maternal composition, with the municipal effect set to zero: posterior median and 95% CrI limits |
| `rr_median`, `rr_lo95`, `rr_hi95` | Residual relative risk, exp(u): posterior median and 95% CrI limits |
| `p_rr_gt1` | Posterior probability that the residual RR exceeds 1 |
| `p_rr_gt1_10` | Posterior probability that the residual RR exceeds 1.10 |
| `p_rr_gt1_25` | Posterior probability that the residual RR exceeds 1.25 |
| `p_rr_lt0_90` | Posterior probability that the residual RR is below 0.90 |
| `excess_median`, `excess_lo95`, `excess_hi95` | Count departure from the covariate-based expectation, E x (RR - 1): posterior median and 95% CrI limits (can be negative) |
| `v1_smr`, `v1_sig_excess` | Standardised mortality ratio and flag of an earlier municipal analysis, kept only for comparison |
| `sh_hospital` | Share of hospital births among births with known place of birth (0 to 1) |
| `sh_miss_esc` | Share of births with missing maternal schooling (0 to 1) |
| `sh_miss_race` | Share of births with missing maternal race/skin colour (0 to 1) |
| `class` | `above`: P(RR > 1.10) at least 0.95; `below`: P(RR < 0.90) at least 0.95; `undetermined`: otherwise |

Notes:

- The registration indicators describe record quality; they do not measure under-registration of births or deaths.
- The classification is the output of an analytical rule, with the operating characteristics reported in the manuscript; it is not a count of deaths attributable to any cause.
- Column definitions follow `v2/src/05_results.py` and `v2/src/03_build_dataset.py`.
