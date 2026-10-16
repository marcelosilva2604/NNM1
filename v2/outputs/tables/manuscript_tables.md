**Table 1.** Covariate effects and residual heterogeneity under three specifications. Incidence-rate ratios per standard deviation of the covariate, with 95% credible intervals; the standard deviation is given in parentheses after each covariate. All models are Poisson with a log-births offset; the spatial models carry a BYM2 random effect restricted to the orthogonal complement of the covariates.

| | Structure, exchangeable | Structure, spatial | Structure and composition, spatial (primary) |
|---|---|---|---|
| Income sub-index of the municipal HDI (SD 0.0807) | 0.922 (0.909–0.935) | 0.913 (0.902–0.923) | 0.983 (0.970–0.996) |
| Education sub-index of the municipal HDI (SD 0.0933) | 0.960 (0.947–0.974) | 0.971 (0.960–0.982) | 1.010 (0.996–1.026) |
| Log resident population (SD 1.18) | 1.030 (1.022–1.038) | 1.029 (1.022–1.035) | 1.003 (0.996–1.010) |
| Mothers younger than 20 years (SD 0.0555) | NA | NA | 1.061 (1.044–1.078) |
| Mothers aged 35 years or older (SD 0.038) | NA | NA | 1.014 (1.000–1.027) |
| Mothers with 7 or fewer years of schooling (SD 0.0975) | NA | NA | 1.026 (1.013–1.038) |
| Non-white mothers (black, brown or indigenous) (SD 0.308) | NA | NA | 1.072 (1.060–1.084) |
| SD of log relative risk across municipalities | 0.169 (0.162–0.176) | 0.169 (0.162–0.175) | 0.157 (0.151–0.164) |
| BYM2 spatial share (ρ, prior parameter) | NA | 0.96 (0.93–0.97) | 0.95 (0.93–0.97) |
| Difference in expected log predictive density (WAIC) from the primary model | -229 | -15 | 0 (reference) |
| Municipalities above expectation (P[RR > 1.10] ≥ 0.95) | 74 | 153 | 112 |
| Municipalities below expectation (P[RR < 0.90] ≥ 0.95) | 32 | 57 | 48 |

**Table 2.** Neonatal mortality and the residual relative risk by macro-region, 2014–2024. Excess deaths are the posterior median of Σ E_i(RR_i − 1) over the municipalities in each set, with 95% credible intervals; the net excess sums over all municipalities of the region and can be negative.

| Region | Municipalities | Neonatal deaths | Rate per 1,000 live births | Above expectation, n | Below expectation, n | Excess deaths in municipalities above expectation | Net excess deaths, all municipalities |
|---|---|---|---|---|---|---|---|
| North | 449 | 33,562 | 10.04 | 21 | 4 | 1,581 (1,414–1,746) | 863 (389–1,322) |
| Northeast | 1,794 | 82,490 | 9.65 | 76 | 29 | 3,607 (3,327–3,886) | 2,644 (2,010–3,288) |
| Central-West | 466 | 21,562 | 8.36 | 3 | 2 | 190 (127–251) | 18 (-312 to 341) |
| Southeast | 1,668 | 92,634 | 7.84 | 8 | 9 | 586 (488–680) | -1,093 (-2,137 to -106) |
| South | 1,188 | 29,717 | 7.14 | 4 | 4 | 145 (105–185) | -758 (-1,099 to -416) |
| Brazil | 5,565 | 259,965 | 8.54 | 112 | 48 | 6,106 (5,740–6,480) | 1,669 (-137 to 3,341) |

**Table 3.** The twelve municipalities with the largest excess deaths among those above expectation with at least 95% probability. Expected deaths come from structure and maternal composition; the relative risk is the residual beyond them. Hospital births is the share of the municipality's live births that took place in a hospital, a record-quality and access indicator.

| Municipality (state) | Neonatal deaths | Expected deaths | Relative risk (95% CrI) | Excess deaths (95% CrI) | Hospital births (%) |
|---|---|---|---|---|---|
| Salvador (BA) | 3,938 | 3,017 | 1.30 (1.25–1.36) | 919 (778–1,062) | 98 |
| São Luís (MA) | 1,811 | 1,419 | 1.27 (1.21–1.33) | 380 (295–467) | 99 |
| Aracaju (SE) | 1,140 | 806 | 1.40 (1.32–1.48) | 322 (258–388) | 99 |
| Macapá (AP) | 1,225 | 887 | 1.35 (1.27–1.42) | 308 (240–374) | 98 |
| Belém (PA) | 2,041 | 1,770 | 1.14 (1.09–1.20) | 255 (164–347) | 100 |
| Santarém (PA) | 930 | 680 | 1.35 (1.27–1.45) | 241 (184–302) | 94 |
| Campos dos Goytacazes (RJ) | 905 | 665 | 1.33 (1.24–1.41) | 217 (163–274) | 100 |
| Itabuna (BA) | 438 | 263 | 1.64 (1.50–1.78) | 167 (132–205) | 99 |
| Petrolina (PE) | 789 | 629 | 1.26 (1.18–1.34) | 162 (111–215) | 97 |
| Feira de Santana (BA) | 1,083 | 916 | 1.17 (1.10–1.24) | 156 (91–222) | 100 |
| Ilhéus (BA) | 388 | 243 | 1.56 (1.42–1.72) | 137 (103–174) | 99 |
| Guarujá (SP) | 506 | 372 | 1.33 (1.22–1.45) | 124 (83–166) | 99 |

**Table S1.** All fitted models: residual heterogeneity, spatial share, health-region standard deviation, fit and classification counts.

| Model | Covariates | Prior on u | SD of log RR (95% CrI) | ρ | Health-region effect, SD on the log scale | Δelpd (WAIC) vs primary | Above | Below | Spearman with primary |
|---|---|---|---|---|---|---|---|---|---|
| iid_struct | structure | exchangeable | 0.169 (0.162–0.176) | NA | NA | -229 | 74 | 32 | 0.76 |
| iid_full | structure + composition | exchangeable | 0.157 (0.150–0.164) | NA | NA | -187 | 59 | 26 | 0.78 |
| bym2_rsr_struct | structure | BYM2, restricted | 0.169 (0.162–0.175) | 0.96 (0.93–0.97) | NA | -15 | 153 | 57 | 0.89 |
| bym2_rsr_full | structure + composition | BYM2, restricted (primary) | 0.157 (0.151–0.164) | 0.95 (0.93–0.97) | NA | +0 | 112 | 48 | 1.00 |
| bym2_rsr_full_hr | structure + composition | BYM2, restricted, + health region | 0.156 (0.150–0.163) | 0.92 (0.87–0.95) | 0.068 (0.055–0.080) | +13 | 128 | 49 | 0.97 |
| bym2_struct | structure | BYM2, unrestricted | 0.176 (0.168–0.184) | 0.96 (0.93–0.98) | NA | -14 | 209 | 51 | 0.80 |
| bym2_full | structure + composition | BYM2, unrestricted | 0.158 (0.152–0.166) | 0.95 (0.92–0.97) | NA | -2 | 106 | 47 | 0.99 |
| bym2_full_hr | structure + composition | BYM2, unrestricted, + health region | 0.158 (0.151–0.165) | 0.92 (0.87–0.95) | 0.067 (0.055–0.080) | +13 | 119 | 53 | 0.96 |
| iid_full_access | + access to care | exchangeable | 0.156 (0.149–0.162) | NA | NA | -185 | 57 | 26 | 0.77 |
| iid_full_intermediate | + antenatal visits, preterm, low birth weight | exchangeable | 0.148 (0.141–0.155) | NA | NA | -115 | 46 | 19 | 0.76 |
| iid_full_completeness | + registration indicators (hospital births, missing schooling and race/colour) | exchangeable | 0.156 (0.150–0.163) | NA | NA | -186 | 56 | 26 | 0.78 |

**Table S2.** Sensitivity covariates (exchangeable models, incidence-rate ratio per SD with 95% CrI).

| Model | Covariate (SD) | IRR (95% CrI) |
|---|---|---|
| iid_full_access | IDHM income sub-index (0.0807) | 0.984 (0.968–1.000) |
| iid_full_access | IDHM education sub-index (0.0933) | 1.006 (0.988–1.023) |
| iid_full_access | log population (1.18) | 0.993 (0.981–1.005) |
| iid_full_access | mothers < 20 y (share) (0.0555) | 1.064 (1.044–1.084) |
| iid_full_access | mothers >= 35 y (share) (0.038) | 1.032 (1.016–1.047) |
| iid_full_access | mothers with <= 7 y schooling (share) (0.0975) | 1.024 (1.009–1.039) |
| iid_full_access | non-white mothers (share) (0.308) | 1.075 (1.061–1.088) |
| iid_full_access | distance to nearest NICU (per 100 km) (0.611) | 1.010 (1.002–1.018) |
| iid_full_access | NICU beds in the municipality (0.246) | 1.013 (1.005–1.021) |
| iid_full_access | births outside the municipality (share) (0.336) | 1.001 (0.990–1.012) |
| iid_full_access | births in a NICU municipality (share) (0.339) | 0.987 (0.976–0.997) |
| iid_full_intermediate | IDHM income sub-index (0.0807) | 0.986 (0.971–1.001) |
| iid_full_intermediate | IDHM education sub-index (0.0933) | 0.993 (0.976–1.010) |
| iid_full_intermediate | log population (1.18) | 0.992 (0.984–1.001) |
| iid_full_intermediate | mothers < 20 y (share) (0.0555) | 1.053 (1.033–1.074) |
| iid_full_intermediate | mothers >= 35 y (share) (0.038) | 1.004 (0.989–1.020) |
| iid_full_intermediate | mothers with <= 7 y schooling (share) (0.0975) | 1.016 (1.001–1.031) |
| iid_full_intermediate | non-white mothers (share) (0.308) | 1.077 (1.063–1.091) |
| iid_full_intermediate | < 7 antenatal visits (share) (0.124) | 1.017 (1.007–1.027) |
| iid_full_intermediate | preterm births (share) (0.0219) | 1.005 (0.997–1.014) |
| iid_full_intermediate | low birth weight (share) (0.0159) | 1.077 (1.065–1.089) |
| iid_full_completeness | IDHM income sub-index (0.0807) | 0.987 (0.972–1.003) |
| iid_full_completeness | IDHM education sub-index (0.0933) | 1.005 (0.987–1.023) |
| iid_full_completeness | log population (1.18) | 0.997 (0.988–1.006) |
| iid_full_completeness | mothers < 20 y (share) (0.0555) | 1.066 (1.046–1.085) |
| iid_full_completeness | mothers >= 35 y (share) (0.038) | 1.031 (1.016–1.047) |
| iid_full_completeness | mothers with <= 7 y schooling (share) (0.0975) | 1.021 (1.007–1.037) |
| iid_full_completeness | non-white mothers (share) (0.308) | 1.079 (1.065–1.093) |
| iid_full_completeness | hospital births (share) (0.052) | 0.991 (0.984–0.997) |
| iid_full_completeness | births with missing schooling (share) (0.032) | 1.001 (0.994–1.009) |
| iid_full_completeness | births with missing race/colour (share) (0.0759) | 0.991 (0.983–0.999) |

**Table S3.** Residual heterogeneity by cause group (exchangeable model with structure and composition, one fit per group) and the cause-specific relative risk in the municipalities classified above expectation on all causes.

| Cause group | Deaths | SD of log RR (95% CrI) | Median cause-specific RR in the 'above' municipalities | Median in the others | Municipalities above expectation on this cause |
|---|---|---|---|---|---|
| Antenatal-care causes | 73,585 | 0.377 (0.359–0.394) | 1.29 | 0.97 | 207 |
| Delivery-care causes | 19,966 | 0.363 (0.339–0.387) | 1.20 | 0.98 | 55 |
| Newborn-care causes | 97,583 | 0.298 (0.286–0.310) | 1.34 | 0.99 | 160 |
| Congenital malformations | 54,606 | 0.127 (0.111–0.144) | 1.03 | 1.00 | 4 |

**Table S4.** Minimum residual relative risk detectable with 80% power (one-sided α = 0.05) from sampling variation alone, by expected deaths, and the number of municipalities at or above each expected count.

| Expected deaths | Minimum detectable RR | Municipalities with at least this expectation | Share of all neonatal deaths in them |
|---|---|---|---|
| 10 | 1.79 | 3,265 | 95% |
| 20 | 1.56 | 2,153 | 89% |
| 50 | 1.35 | 942 | 74% |
| 100 | 1.25 | 448 | 61% |
| 200 | 1.18 | 196 | 48% |
| 500 | 1.11 | 63 | 32% |
| 1,000 | 1.08 | 21 | 21% |
| 2,000 | 1.06 | 7 | 13% |
| 4,000 | 1.04 | 2 | 7% |

**Posterior predictive check (primary model).** 1.2% of municipalities fall outside their 95% posterior predictive interval; share of municipalities with zero deaths 1.6% observed versus 1.8% predicted; largest count 12,293 observed versus 12,294 predicted (median); variance-to-mean ratio 1,195 versus 1,195.
