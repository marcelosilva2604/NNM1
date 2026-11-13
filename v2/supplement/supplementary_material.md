# Supplementary material

**Maternal composition and area-level variation in neonatal mortality in Brazil, 2014-2024: an ecological study with Bayesian spatial analysis**

## eMethods 1. Adjacency graph and BYM2 scaling

Queen contiguity was computed on the IBGE municipal mesh (46 invalid geometries repaired before use) for the 5,565 analysed municipalities, giving 16,129 undirected neighbour pairs (mean 5.8 neighbours; minimum 1, maximum 23). Four municipalities had no contiguous neighbour (offshore or coastal polygons that do not touch) and the graph had several components; each isolated component was joined to its nearest municipality by centroid distance (Brazil Polyconic projection, EPSG:5880) until the graph was connected, as an intrinsic conditional autoregressive prior is improper on every disconnected component. Nazária (PI), emancipated in 2008 and absent from the mesh, was placed at the median coordinates of its eight CNES establishments and joined to its three nearest neighbours. The artificial joins numbered 10, with a median length of 19.1 km; the longest was 374.8 km. Component-specific constraints and scaling are an alternative that was not fitted. The BYM2 scaling factor, the geometric mean of the diagonal of the generalised inverse of the graph Laplacian (reference 8 of the main text), was 3.58.

## eMethods 2. Model specification, priors, sampling and diagnostics

For municipality $i$ with $b_i$ live births, standardised covariates $x_i$, and $y_i$ neonatal deaths:

$$y_i \sim \mathrm{Poisson}\big(b_i \exp(\beta_0 + x_i'\beta + u_i)\big)$$

$$u^{*} = \sigma_u \big( (\rho/s)^{1/2} \phi + (1-\rho)^{1/2} \theta \big)$$

with $\phi$ an intrinsic conditional autoregressive field on $W$ (ICAR) and $\theta_i \sim N(0, 1)$ exchangeable, and

$$u = \big(I - X_a (X_a' X_a)^{-1} X_a'\big)\, u^{*}, \qquad X_a = [\mathbf{1}, X]$$

with $s$ the scaling factor and $W$ the adjacency matrix.

The projection in the last line is restricted spatial regression (references 10 and 11 of the main text): the random effect is confined to the orthogonal complement of the covariate space, so that $u$ is orthogonal to the covariates and represents the part of the municipal variation not aligned with them. Centring $u$ on the log scale does not centre $\exp(u)$, so the national net count departure need not be zero. Restricted methods can understate the uncertainty of the covariate coefficients (reference 12); the coefficients are not the primary estimand here.

Without the projection (unrestricted BYM2), the smooth spatial field competes with the covariates: the income mortality rate ratio per standard deviation was 0.95 with structure only and 0.98 with composition, against 0.91 and 0.98 under the projection. The unrestricted model with composition produced a similar residual map (Spearman correlation 0.99 with the restricted fit; 108 versus 110 municipalities above expectation) and a similar fit (difference in expected log predictive density +0.3). The coefficients therefore depend on the specification more than the map does.

| Model | Income, mortality rate ratio per SD (95% CrI) | Education, mortality rate ratio per SD (95% CrI) |
|---|---|---|
| Structure, restricted spatial | 0.913 (95% CrI 0.902, 0.924) | 0.971 (95% CrI 0.960, 0.982) |
| Structure, unrestricted spatial | 0.955 (95% CrI 0.935, 0.974) | 0.977 (95% CrI 0.962, 0.993) |
| Structure and composition, restricted spatial (primary) | 0.983 (95% CrI 0.970, 0.996) | 1.011 (95% CrI 0.996, 1.026) |
| Structure and composition, unrestricted spatial | 0.975 (95% CrI 0.955, 0.997) | 1.018 (95% CrI 0.998, 1.037) |

Priors: $\beta_0 \sim N(\log(\bar y / \bar b), 1)$, $\beta_j \sim N(0, 0.5)$ on the standardised scale, $\sigma_u \sim$ half-normal(0.5), $\rho \sim$ Beta(1, 1), and in the health-region model $\gamma_r = \sigma_{hr} z_r$ with $z_r \sim N(0,1)$ and $\sigma_{hr} \sim$ half-normal(0.5), added to $u^{*}$ before the projection. The exchangeable models replace the convolution by $\sigma_u \theta$.

Sampling: PyMC 6.2 (reference 14), No-U-Turn sampler, target acceptance 0.9. The primary model used four chains with 4,000 warm-up and 4,000 retained draws per chain (seed 20261013); it had 0 divergences, maximum rank-normalised R-hat 1.0024 for scalar parameters and 1.0015 for municipal effects, and minimum bulk effective sample size 1,449 and 7,235, respectively. Period-specific fits used seeds 20261014 and 20261015. Settings and diagnostics of the other fits are given below.

| Fit | Warm-up | Retained draws | Divergences | Max R-hat, scalar parameters | Min bulk ESS, scalar parameters | Max R-hat, municipal effects | Min bulk ESS, municipal effects |
|---|---|---|---|---|---|---|---|
| Structure, exchangeable | 2,000 | 1,000 | 0 | 1.007 | 655 | 1.009 | 572 |
| Structure and composition, exchangeable | 4,000 | 4,000 | 0 | 1.001 | 2,026 | 1.002 | 1,446 |
| Structure, restricted spatial | 4,000 | 8,000 | 0 | 1.002 | 759 | 1.001 | 8,145 |
| Structure and composition, restricted spatial (primary) | 4,000 | 4,000 | 0 | 1.002 | 1,449 | 1.002 | 7,235 |
| Primary plus health-region effect | 4,000 | 4,000 | 0 | 1.004 | 722 | 1.002 | 5,385 |
| Structure, unrestricted spatial | 4,000 | 4,000 | 0 | 1.002 | 1,042 | 1.002 | 7,064 |
| Structure and composition, unrestricted spatial | 4,000 | 8,000 | 0 | 1.009 | 615 | 1.004 | 1,713 |
| Unrestricted plus health-region effect | 4,000 | 8,000 | 0 | 1.016 | 354 | 1.005 | 1,125 |
| Exchangeable plus access indicators | 4,000 | 8,000 | 0 | 1.001 | 11,725 | 1.001 | 8,506 |
| Exchangeable plus pregnancy and birth characteristics | 4,000 | 8,000 | 0 | 1.001 | 5,419 | 1.002 | 3,317 |
| Exchangeable plus registration indicators | 2,000 | 1,000 | 0 | 1.004 | 1,139 | 1.010 | 798 |
| Primary specification, 2014-2019 | 4,000 | 4,000 | 0 | 1.010 | 510 | 1.003 | 3,253 |
| Primary specification, 2020-2024 | 4,000 | 4,000 | 0 | 1.011 | 529 | 1.003 | 1,931 |
| Scenario fit, extreme: schooling q = 0, race/skin colour q = 0 | 4,000 | 4,000 | 0 | 1.001 | 1,538 | 1.002 | 11,956 |
| Scenario fit, extreme: schooling q = 0, race/skin colour q = 1 | 4,000 | 4,000 | 0 | 1.004 | 1,349 | 1.002 | 9,340 |
| Scenario fit, extreme: schooling q = 1, race/skin colour q = 0 | 4,000 | 4,000 | 0 | 1.002 | 1,286 | 1.002 | 10,742 |
| Scenario fit, extreme: schooling q = 1, race/skin colour q = 1 | 4,000 | 4,000 | 0 | 1.004 | 1,328 | 1.002 | 7,495 |
| Scenario fit, intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | 4,000 | 4,000 | 0 | 1.002 | 1,602 | 1.002 | 7,604 |
| Scenario fit, intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | 4,000 | 4,000 | 0 | 1.004 | 1,422 | 1.002 | 8,529 |
| Scenario fit, intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | 4,000 | 4,000 | 0 | 1.001 | 1,540 | 1.002 | 9,807 |
| Scenario fit, intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | 4,000 | 4,000 | 0 | 1.002 | 1,477 | 1.002 | 8,505 |

Simulation fits. In the 40 datasets with residual heterogeneity, the maximum R-hat of the municipal effects was 1.009 and the minimum bulk ESS was 876; hyperparameter R-hat exceeded 1.01 in 33 fits and hyperparameter ESS was below 400 in 38. In the 12 null datasets (no residual heterogeneity), no municipality was flagged and the posterior median realised standard deviation of the municipal effect was 0.003 to 0.016, but mixing of the municipal effects was weak in all 12 fits (R-hat 1.03 to 1.23; bulk ESS 12 to 85) even after the sampling was extended to 4,000 warm-up and 4,000 retained draws per chain, a pattern expected when the scale of the municipal effect is close to zero. The 12 null fits therefore support the absence of flags only qualitatively. The maximum exceedance probability in the null fits was not stored.

Numerical classification uncertainty. The rule is fixed at P(RR > 1.10) ≥ 0.95. Exceedance probabilities are estimated from correlated draws; Monte Carlo standard errors were computed from 20 batch means per chain of the primary posterior (16,000 draws). Eleven municipalities lay within 1.96 Monte Carlo standard errors of the 0.95 threshold, and the median standard error of the exceedance probabilities was 0.003. These errors quantify computational uncertainty, not posterior uncertainty, and municipalities near the threshold can change classification in an independent run.

Distance to the nearest municipality with neonatal intensive care beds was computed from the establishments registered in the last year each municipality appears in CNES; presence of beds uses any year of 2014-2024. Heterogeneity is reported as the standard deviation of $u$ across municipalities within a posterior draw (median over draws), because under a scaled ICAR the BYM2 scale parameter $\sigma_u$ is not the realised spread of the field.

## eMethods 3. Maternal age and missing-input sensitivity

The 10-60-year admissible range was applied to the known-age denominator and to both age-category numerators. Aggregation of the individual-record files was reconciled with the municipal-year counts of births and known-age denominators (61,435 municipality-year rows; zero absolute discrepancy). In the analysed set, applying the range to the numerators removed 18 records from the younger-than-20 numerator and 372 from the age-35-or-older numerator, affecting 227 municipalities. Births, deaths, known-age denominators, non-age inputs and the adjacency graph were unchanged. No known-value denominator was zero in the full period or in either separate period, so no imputation of these proportions was performed.

Missingness of the maternal fields, nationally and by municipality:

| Variable | Missing or unusable records, national | National percentage | Municipalities with more than 5.0% missing | Maximum municipal percentage | Flagged municipalities with more than 5.0% missing | Maximum width of the extreme-allocation input range among flagged municipalities (percentage points) |
|---|---|---|---|---|---|---|
| Maternal age | 771 | 0.00 | 0 | 0.6 | 0 | 0.6 |
| Maternal schooling | 368,213 | 1.21 | 343 | 65.4 | 14 | 19.7 |
| Maternal race/skin colour | 1,101,219 | 3.62 | 1,035 | 90.6 | 26 | 21.6 |
| Antenatal visits | 160,625 | 0.53 | 13 | 14.6 | 0 | 4.4 |
| Gestational age | 477,241 | 1.57 | 469 | 43.2 | 17 | 18.8 |
| Birthweight | 8,603 | 0.03 | 2 | 14.1 | 1 | 14.1 |

For each schooling or race/skin colour category, let N denote registered births, K the known-value count, C the known category count and M = N − K. Extreme scenarios used (C + qM)/N with q equal to zero or one, in all four combinations across the two variables. Intermediate scenarios used q = ωp/(1 − p + ωp), where p = C/K and ω was 0.5 or 2, in all four combinations; with ω = 1 the input equals C/K exactly. These are deterministic input scenarios; expected category allocations need not be integers. They make no claim that missingness follows any of these mechanisms, and their results do not bound all possible posteriors.

Every scenario recomputed the across-municipality means and population standard deviations, the standardised design matrix and its orthogonal projection, and refitted the primary model; births, deaths, MHDI, population and graph were held fixed. Coefficients were also expressed per one standard deviation of each covariate in the corrected primary input, using exp(β × SD_reference / SD_scenario). Selected sums use each scenario's own selected set.

Extreme scenarios and intermediate scenarios are input stress tests, not mathematical bounds on coefficients, risk estimates, classifications or count departures. In the two extreme scenarios that allocated all missing race/skin colour records outside the Black, brown or Indigenous category (q = 0), the national net count departure interval excluded zero.

| Scenario | Flagged above expectation | Retained from primary | Removed | Added | Spearman correlation of RR with primary | SD of municipal effect, median (95% CrI) | National net count departure, median (95% CrI) | Selected count departure, median (95% CrI) |
|---|---|---|---|---|---|---|---|---|
| Extreme: schooling q = 0, race/skin colour q = 0 | 112 | 106 | 4 | 6 | 0.996 | 0.158 (95% CrI 0.151, 0.164) | 1,814 (95% CrI 83, 3,528) | 5,898 (95% CrI 5,552, 6,247) |
| Extreme: schooling q = 0, race/skin colour q = 1 | 114 | 110 | 0 | 4 | 0.999 | 0.158 (95% CrI 0.151, 0.164) | 1,681 (95% CrI -37, 3,400) | 6,139 (95% CrI 5,770, 6,520) |
| Extreme: schooling q = 1, race/skin colour q = 0 | 107 | 104 | 6 | 3 | 0.996 | 0.158 (95% CrI 0.151, 0.164) | 1,871 (95% CrI 146, 3,601) | 5,871 (95% CrI 5,535, 6,220) |
| Extreme: schooling q = 1, race/skin colour q = 1 | 110 | 109 | 1 | 1 | 0.999 | 0.158 (95% CrI 0.151, 0.165) | 1,715 (95% CrI -25, 3,461) | 6,134 (95% CrI 5,763, 6,507) |
| Intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | 110 | 107 | 3 | 3 | 1.000 | 0.157 (95% CrI 0.151, 0.164) | 1,684 (95% CrI -49, 3,408) | 5,849 (95% CrI 5,502, 6,197) |
| Intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | 110 | 109 | 1 | 1 | 1.000 | 0.157 (95% CrI 0.151, 0.164) | 1,681 (95% CrI -66, 3,416) | 6,108 (95% CrI 5,732, 6,478) |
| Intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | 108 | 106 | 4 | 2 | 1.000 | 0.157 (95% CrI 0.151, 0.164) | 1,700 (95% CrI -61, 3,434) | 5,845 (95% CrI 5,502, 6,195) |
| Intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | 110 | 109 | 1 | 1 | 1.000 | 0.157 (95% CrI 0.151, 0.164) | 1,698 (95% CrI -28, 3,434) | 6,114 (95% CrI 5,737, 6,490) |

| Scenario | Covariate | Mortality rate ratio per SD (95% CrI), per scenario SD | Per SD of the corrected primary input (95% CrI) |
|---|---|---|---|
| Extreme: schooling q = 0, race/skin colour q = 0 | MHDI income subindex | 0.977 (95% CrI 0.964, 0.990) | 0.977 (95% CrI 0.964, 0.990) |
| Extreme: schooling q = 0, race/skin colour q = 0 | MHDI education subindex | 1.010 (95% CrI 0.995, 1.025) | 1.010 (95% CrI 0.995, 1.025) |
| Extreme: schooling q = 0, race/skin colour q = 0 | Logarithm of resident population | 1.005 (95% CrI 0.997, 1.012) | 1.005 (95% CrI 0.997, 1.012) |
| Extreme: schooling q = 0, race/skin colour q = 0 | Mothers younger than 20 years | 1.060 (95% CrI 1.043, 1.077) | 1.060 (95% CrI 1.043, 1.077) |
| Extreme: schooling q = 0, race/skin colour q = 0 | Mothers aged 35 years or older | 1.015 (95% CrI 1.001, 1.028) | 1.015 (95% CrI 1.001, 1.028) |
| Extreme: schooling q = 0, race/skin colour q = 0 | Mothers with 7 years of schooling or less | 1.024 (95% CrI 1.012, 1.037) | 1.025 (95% CrI 1.012, 1.038) |
| Extreme: schooling q = 0, race/skin colour q = 0 | Black, brown or Indigenous mothers | 1.068 (95% CrI 1.057, 1.079) | 1.071 (95% CrI 1.059, 1.083) |
| Extreme: schooling q = 0, race/skin colour q = 1 | MHDI income subindex | 0.982 (95% CrI 0.969, 0.995) | 0.982 (95% CrI 0.969, 0.995) |
| Extreme: schooling q = 0, race/skin colour q = 1 | MHDI education subindex | 1.009 (95% CrI 0.994, 1.024) | 1.009 (95% CrI 0.994, 1.024) |
| Extreme: schooling q = 0, race/skin colour q = 1 | Logarithm of resident population | 1.004 (95% CrI 0.997, 1.011) | 1.004 (95% CrI 0.997, 1.011) |
| Extreme: schooling q = 0, race/skin colour q = 1 | Mothers younger than 20 years | 1.063 (95% CrI 1.046, 1.080) | 1.063 (95% CrI 1.046, 1.080) |
| Extreme: schooling q = 0, race/skin colour q = 1 | Mothers aged 35 years or older | 1.015 (95% CrI 1.001, 1.028) | 1.015 (95% CrI 1.001, 1.028) |
| Extreme: schooling q = 0, race/skin colour q = 1 | Mothers with 7 years of schooling or less | 1.023 (95% CrI 1.010, 1.035) | 1.023 (95% CrI 1.011, 1.036) |
| Extreme: schooling q = 0, race/skin colour q = 1 | Black, brown or Indigenous mothers | 1.070 (95% CrI 1.058, 1.082) | 1.070 (95% CrI 1.059, 1.083) |
| Extreme: schooling q = 1, race/skin colour q = 0 | MHDI income subindex | 0.979 (95% CrI 0.966, 0.992) | 0.979 (95% CrI 0.966, 0.992) |
| Extreme: schooling q = 1, race/skin colour q = 0 | MHDI education subindex | 1.009 (95% CrI 0.994, 1.024) | 1.009 (95% CrI 0.994, 1.024) |
| Extreme: schooling q = 1, race/skin colour q = 0 | Logarithm of resident population | 1.004 (95% CrI 0.997, 1.011) | 1.004 (95% CrI 0.997, 1.011) |
| Extreme: schooling q = 1, race/skin colour q = 0 | Mothers younger than 20 years | 1.062 (95% CrI 1.046, 1.079) | 1.062 (95% CrI 1.046, 1.079) |
| Extreme: schooling q = 1, race/skin colour q = 0 | Mothers aged 35 years or older | 1.016 (95% CrI 1.003, 1.029) | 1.016 (95% CrI 1.003, 1.029) |
| Extreme: schooling q = 1, race/skin colour q = 0 | Mothers with 7 years of schooling or less | 1.024 (95% CrI 1.012, 1.035) | 1.023 (95% CrI 1.012, 1.033) |
| Extreme: schooling q = 1, race/skin colour q = 0 | Black, brown or Indigenous mothers | 1.068 (95% CrI 1.057, 1.080) | 1.071 (95% CrI 1.060, 1.083) |
| Extreme: schooling q = 1, race/skin colour q = 1 | MHDI income subindex | 0.983 (95% CrI 0.970, 0.996) | 0.983 (95% CrI 0.970, 0.996) |
| Extreme: schooling q = 1, race/skin colour q = 1 | MHDI education subindex | 1.005 (95% CrI 0.990, 1.020) | 1.005 (95% CrI 0.990, 1.020) |
| Extreme: schooling q = 1, race/skin colour q = 1 | Logarithm of resident population | 1.004 (95% CrI 0.997, 1.012) | 1.004 (95% CrI 0.997, 1.012) |
| Extreme: schooling q = 1, race/skin colour q = 1 | Mothers younger than 20 years | 1.067 (95% CrI 1.051, 1.084) | 1.067 (95% CrI 1.051, 1.084) |
| Extreme: schooling q = 1, race/skin colour q = 1 | Mothers aged 35 years or older | 1.016 (95% CrI 1.003, 1.030) | 1.016 (95% CrI 1.003, 1.030) |
| Extreme: schooling q = 1, race/skin colour q = 1 | Mothers with 7 years of schooling or less | 1.018 (95% CrI 1.006, 1.029) | 1.017 (95% CrI 1.006, 1.028) |
| Extreme: schooling q = 1, race/skin colour q = 1 | Black, brown or Indigenous mothers | 1.069 (95% CrI 1.057, 1.081) | 1.069 (95% CrI 1.058, 1.081) |
| Intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | MHDI income subindex | 0.983 (95% CrI 0.970, 0.996) | 0.983 (95% CrI 0.970, 0.996) |
| Intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | MHDI education subindex | 1.010 (95% CrI 0.995, 1.026) | 1.010 (95% CrI 0.995, 1.026) |
| Intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | Logarithm of resident population | 1.003 (95% CrI 0.996, 1.010) | 1.003 (95% CrI 0.996, 1.010) |
| Intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | Mothers younger than 20 years | 1.060 (95% CrI 1.044, 1.077) | 1.060 (95% CrI 1.044, 1.077) |
| Intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | Mothers aged 35 years or older | 1.014 (95% CrI 1.000, 1.027) | 1.014 (95% CrI 1.000, 1.027) |
| Intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | Mothers with 7 years of schooling or less | 1.025 (95% CrI 1.013, 1.038) | 1.026 (95% CrI 1.013, 1.038) |
| Intermediate: odds multiplier 0.5 (schooling), 0.5 (race/skin colour) | Black, brown or Indigenous mothers | 1.072 (95% CrI 1.060, 1.084) | 1.072 (95% CrI 1.060, 1.084) |
| Intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | MHDI income subindex | 0.982 (95% CrI 0.969, 0.996) | 0.982 (95% CrI 0.969, 0.996) |
| Intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | MHDI education subindex | 1.010 (95% CrI 0.995, 1.026) | 1.010 (95% CrI 0.995, 1.026) |
| Intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | Logarithm of resident population | 1.003 (95% CrI 0.996, 1.011) | 1.003 (95% CrI 0.996, 1.011) |
| Intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | Mothers younger than 20 years | 1.061 (95% CrI 1.044, 1.078) | 1.061 (95% CrI 1.044, 1.078) |
| Intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | Mothers aged 35 years or older | 1.014 (95% CrI 1.000, 1.027) | 1.014 (95% CrI 1.000, 1.027) |
| Intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | Mothers with 7 years of schooling or less | 1.025 (95% CrI 1.013, 1.037) | 1.025 (95% CrI 1.013, 1.038) |
| Intermediate: odds multiplier 0.5 (schooling), 2 (race/skin colour) | Black, brown or Indigenous mothers | 1.071 (95% CrI 1.059, 1.083) | 1.071 (95% CrI 1.059, 1.083) |
| Intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | MHDI income subindex | 0.983 (95% CrI 0.970, 0.996) | 0.983 (95% CrI 0.970, 0.996) |
| Intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | MHDI education subindex | 1.011 (95% CrI 0.996, 1.026) | 1.011 (95% CrI 0.996, 1.026) |
| Intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | Logarithm of resident population | 1.003 (95% CrI 0.996, 1.010) | 1.003 (95% CrI 0.996, 1.010) |
| Intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | Mothers younger than 20 years | 1.060 (95% CrI 1.044, 1.077) | 1.060 (95% CrI 1.044, 1.077) |
| Intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | Mothers aged 35 years or older | 1.014 (95% CrI 1.000, 1.027) | 1.014 (95% CrI 1.000, 1.027) |
| Intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | Mothers with 7 years of schooling or less | 1.026 (95% CrI 1.013, 1.039) | 1.026 (95% CrI 1.013, 1.038) |
| Intermediate: odds multiplier 2 (schooling), 0.5 (race/skin colour) | Black, brown or Indigenous mothers | 1.072 (95% CrI 1.060, 1.084) | 1.072 (95% CrI 1.060, 1.084) |
| Intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | MHDI income subindex | 0.983 (95% CrI 0.970, 0.996) | 0.983 (95% CrI 0.970, 0.996) |
| Intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | MHDI education subindex | 1.010 (95% CrI 0.995, 1.026) | 1.010 (95% CrI 0.995, 1.026) |
| Intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | Logarithm of resident population | 1.003 (95% CrI 0.996, 1.011) | 1.003 (95% CrI 0.996, 1.011) |
| Intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | Mothers younger than 20 years | 1.062 (95% CrI 1.045, 1.079) | 1.062 (95% CrI 1.045, 1.079) |
| Intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | Mothers aged 35 years or older | 1.014 (95% CrI 1.001, 1.027) | 1.014 (95% CrI 1.001, 1.027) |
| Intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | Mothers with 7 years of schooling or less | 1.025 (95% CrI 1.013, 1.038) | 1.025 (95% CrI 1.013, 1.037) |
| Intermediate: odds multiplier 2 (schooling), 2 (race/skin colour) | Black, brown or Indigenous mothers | 1.071 (95% CrI 1.059, 1.083) | 1.071 (95% CrI 1.059, 1.083) |

## eTable 1. Alternative specifications

The BYM2 mixing parameter and the health-region scale are defined before projection and are not fractions of the realised variance. WAIC differences are descriptive; standard errors were not computed, and differences of a few units are not interpreted.

| Model | SD of residual log RR (95% CrI) | BYM2 mixing parameter (95% CrI) | Health-region SD (95% CrI) | Difference in expected log predictive density (WAIC) | Above expectation | Below expectation | Spearman correlation with primary |
|---|---|---|---|---|---|---|---|
| Structure, exchangeable | 0.169 (95% CrI 0.162, 0.176) | NA | NA | -231 | 74 | 32 | 0.76 |
| Structure and composition, exchangeable | 0.157 (95% CrI 0.150, 0.164) | NA | NA | -188 | 57 | 28 | 0.78 |
| Structure, restricted spatial | 0.168 (95% CrI 0.162, 0.175) | 0.96 (95% CrI 0.93, 0.97) | NA | -15 | 151 | 58 | 0.89 |
| Structure and composition, restricted spatial (primary) | 0.157 (95% CrI 0.151, 0.164) | 0.95 (95% CrI 0.93, 0.97) | NA | +0 | 110 | 50 | 1.00 |
| Primary plus health-region effect | 0.156 (95% CrI 0.150, 0.163) | 0.92 (95% CrI 0.87, 0.95) | 0.068 (95% CrI 0.055, 0.080) | +14 | 128 | 46 | 0.97 |
| Structure, unrestricted spatial | 0.176 (95% CrI 0.168, 0.184) | 0.96 (95% CrI 0.93, 0.98) | NA | -16 | 205 | 50 | 0.80 |
| Structure and composition, unrestricted spatial | 0.159 (95% CrI 0.152, 0.166) | 0.95 (95% CrI 0.93, 0.97) | NA | +0 | 108 | 47 | 0.99 |
| Unrestricted plus health-region effect | 0.158 (95% CrI 0.151, 0.165) | 0.92 (95% CrI 0.87, 0.95) | 0.068 (95% CrI 0.055, 0.081) | +13 | 119 | 50 | 0.96 |
| Exchangeable plus access indicators | 0.156 (95% CrI 0.149, 0.163) | NA | NA | -187 | 55 | 25 | 0.77 |
| Exchangeable plus pregnancy and birth characteristics | 0.148 (95% CrI 0.141, 0.155) | NA | NA | -119 | 45 | 18 | 0.76 |
| Exchangeable plus registration indicators | 0.156 (95% CrI 0.150, 0.163) | NA | NA | -185 | 55 | 26 | 0.78 |

## eTable 2. Examples among municipalities classified above expectation

The twelve largest municipal count-departure estimates among the municipalities classified above expectation (P(RR > 1.10) ≥ 0.95) in the primary model. Expected deaths are those expected from structure and maternal composition with the municipal effect set to zero. These examples are conditional on the classification.

| Municipality (state) | Neonatal deaths | Expected deaths | Residual RR (95% CrI) | Selected count departure (95% CrI) | Hospital births (%) |
|---|---|---|---|---|---|
| Salvador (BA) | 3,938 | 3,016 | 1.30 (95% CrI 1.25, 1.36) | 918 (95% CrI 777, 1,062) | 97.9 |
| São Luís (MA) | 1,811 | 1,418 | 1.27 (95% CrI 1.21, 1.33) | 379 (95% CrI 293, 467) | 99.4 |
| Aracaju (SE) | 1,140 | 806 | 1.40 (95% CrI 1.32, 1.49) | 324 (95% CrI 258, 391) | 99.4 |
| Macapá (AP) | 1,225 | 887 | 1.35 (95% CrI 1.27, 1.43) | 309 (95% CrI 242, 377) | 98.1 |
| Belém (PA) | 2,041 | 1,769 | 1.14 (95% CrI 1.09, 1.20) | 256 (95% CrI 166, 349) | 99.6 |
| Santarém (PA) | 930 | 680 | 1.35 (95% CrI 1.27, 1.44) | 240 (95% CrI 183, 300) | 93.9 |
| Campos dos Goytacazes (RJ) | 905 | 665 | 1.33 (95% CrI 1.24, 1.42) | 217 (95% CrI 162, 276) | 99.6 |
| Itabuna (BA) | 438 | 263 | 1.64 (95% CrI 1.50, 1.78) | 167 (95% CrI 131, 206) | 99.4 |
| Petrolina (PE) | 789 | 629 | 1.26 (95% CrI 1.18, 1.34) | 162 (95% CrI 111, 216) | 97.1 |
| Feira de Santana (BA) | 1,083 | 916 | 1.17 (95% CrI 1.10, 1.24) | 157 (95% CrI 93, 221) | 99.5 |
| Ilhéus (BA) | 388 | 243 | 1.56 (95% CrI 1.42, 1.71) | 137 (95% CrI 102, 173) | 98.7 |
| Guarujá (SP) | 506 | 372 | 1.33 (95% CrI 1.22, 1.45) | 123 (95% CrI 83, 167) | 99.4 |

## eTable 3. E-values

Computed from unrounded estimates with the formula and functions of the R package EValue (version 4.1.4). These summaries concern municipal covariate contrasts, not residual classifications. The coefficient intervals are conditional on restricted spatial regression and may understate uncertainty.

| Covariate | Mortality rate ratio per SD (95% CrI) | E-value | E-value at the interval limit closer to 1 |
|---|---|---|---|
| MHDI income subindex | 0.9827 (95% CrI 0.9695, 0.9959) | 1.152 | 1.068 |
| MHDI education subindex | 1.0106 (95% CrI 0.9955, 1.0262) | 1.114 | 1.000 |
| Logarithm of resident population | 1.0030 (95% CrI 0.9957, 1.0104) | 1.058 | 1.000 |
| Mothers younger than 20 years | 1.0608 (95% CrI 1.0441, 1.0777) | 1.315 | 1.259 |
| Mothers aged 35 years or older | 1.0137 (95% CrI 1.0003, 1.0274) | 1.131 | 1.019 |
| Mothers with 7 years of schooling or less | 1.0256 (95% CrI 1.0131, 1.0382) | 1.188 | 1.128 |
| Black, brown or Indigenous mothers | 1.0716 (95% CrI 1.0596, 1.0834) | 1.348 | 1.311 |
