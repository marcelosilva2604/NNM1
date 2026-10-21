# Supplementary material

**Neonatal mortality beyond socioeconomic structure and maternal composition across Brazilian municipalities, 2014-2024: a Bayesian spatial analysis**

## S1. Adjacency graph and BYM2 scaling

Queen contiguity was computed on the IBGE municipal mesh (46 invalid geometries repaired before use) for the 5,565 analyzed municipalities, giving 16,129 undirected neighbor pairs (mean 5.8 neighbors; minimum 1, maximum 23). Four municipalities had no contiguous neighbor (offshore or coastal polygons that do not touch) and the graph had several components; each isolated component was joined to its nearest municipality by centroid distance (Brazil Polyconic projection, EPSG:5880) until the graph was connected, as an intrinsic conditional autoregressive prior is improper on every disconnected component. Nazária (PI), emancipated in 2008 and absent from the mesh, was placed at the median coordinates of its eight CNES establishments and joined to its three nearest neighbors. The artificial joins numbered 10, with a median length of 19.1 km; the longest was 374.8 km. Component-specific constraints and scaling are an alternative that was not fitted. The BYM2 scaling factor, the geometric mean of the diagonal of the generalized inverse of the graph Laplacian (Riebler et al., 2016), was 3.58.

## S2. Model specification, priors, and diagnostics

For municipality $i$ with $b_i$ live births, standardized covariates $x_i$, and $y_i$ neonatal deaths:

$$y_i \sim \mathrm{Poisson}\big(b_i \exp(\beta_0 + x_i'\beta + u_i)\big)$$

$$u^{*} = \sigma_u \big( (\rho/s)^{1/2} \phi + (1-\rho)^{1/2} \theta \big)$$

with $\phi$ an intrinsic conditional autoregressive field on $W$ (ICAR) and $\theta_i \sim N(0, 1)$ exchangeable, and

$$u = \big(I - X_a (X_a' X_a)^{-1} X_a'\big)\, u^{*}, \qquad X_a = [\mathbf{1}, X]$$

with $s$ the scaling factor and $W$ the adjacency matrix. The projection in the last line is restricted spatial regression (Hodges & Reich, 2010; Hughes & Haran, 2013): the random effect is confined to the orthogonal complement of the covariate space, so that $u$ is orthogonal to the covariate space and is what the covariates do not predict; centering $u$ on the log scale does not center $\exp(u)$, so the national net excess need not be zero. Without the projection, the smooth spatial field absorbed the covariates (income IRR per SD 0.955 with structure only and 0.976 with composition, against 0.913 and 0.983 under the projection; ρ 0.95) while producing a similar residual map (Spearman 0.99 with the restricted fit; 106 versus 112 municipalities above expectation) and a similar fit (Table S1, Δelpd −2); the projection changes the interpretation of the coefficients more than the map. Khan & Calder (2022) show that restricted methods can understate the uncertainty of the covariate coefficients; the coefficients are not the primary estimand here.

Priors: $\beta_0 \sim N(\log(\bar y / \bar b), 1)$, $\beta_j \sim N(0, 0.5)$ on the standardized scale, $\sigma_u \sim$ half-normal(0.5), $\rho \sim$ Beta(1, 1), and in the health-region model $\gamma_r = \sigma_{hr} z_r$ with $z_r \sim N(0,1)$ and $\sigma_{hr} \sim$ half-normal(0.5), added to $u^{*}$ before the projection. The exchangeable models replace the convolution by $\sigma_u \theta$.

Sampling: PyMC 6.2 (Abril-Pla et al., 2023), No-U-Turn sampler, target acceptance 0.9, four chains of 2,000 warm-up and 1,000 retained draws (4,000 draws). Primary model: 0 divergent transitions; maximum split $\hat R$ 1.008 (hyperparameters 1.007); minimum bulk effective sample size 346 for hyperparameters and 1,807 for the municipality effects. All eleven full-period fits had 0 divergences and maximum $\hat R$ ≤ 1.055 (the unrestricted BYM2 with composition, whose covariate coefficients mix slowly for the reason given above; its residual map is similar to that of the restricted fit, Spearman 0.99; maximum $\hat R$ for the municipality effects 1.02).

Distance to the nearest municipality with neonatal intensive care beds was computed from the establishments registered in the last year each municipality appears in CNES; presence of beds uses any year of 2014-2024. Heterogeneity is reported as the standard deviation of $u$ across municipalities within a posterior draw (median over draws), because under a scaled ICAR the BYM2 scale parameter $\sigma_u$ (0.37 in the primary model) is not the realized spread of the field.

## S3. Tables

**Table S1.** All fitted models: residual heterogeneity, spatial share, health-region standard deviation, fit, and classification counts. Δelpd is the difference in expected log predictive density (WAIC) from the primary model; positive values favor the row model. Differences are shown without standard errors, and differences of a few units are not interpreted. The spatial share ρ and the health-region scale are parameters of the prior defined before projection, not shares of the realized variance of the projected field. NA, not applicable.

| Model | Covariates | Prior on u | SD of log RR (95% CrI) | ρ | Health-region effect, SD on the log scale | Δelpd vs primary | Above | Below | Spearman with primary |
|---|---|---|---|---|---|---|---|---|---|
| Structure, exchangeable | structure | exchangeable | 0.169 (0.162-0.176) | NA | NA | −229 | 74 | 32 | 0.76 |
| Structure + composition, exchangeable | structure + composition | exchangeable | 0.157 (0.150-0.164) | NA | NA | −187 | 59 | 26 | 0.78 |
| Structure, spatial restricted | structure | BYM2, restricted | 0.169 (0.162-0.175) | 0.96 (0.93-0.97) | NA | −15 | 153 | 57 | 0.89 |
| **Primary** | structure + composition | BYM2, restricted | 0.157 (0.151-0.164) | 0.95 (0.93-0.97) | NA | 0 | 112 | 48 | 1.00 |
| Primary + health region | structure + composition | BYM2, restricted, + health region | 0.156 (0.150-0.163) | 0.92 (0.87-0.95) | 0.068 (0.055-0.080) | +13 | 128 | 49 | 0.97 |
| Structure, spatial unrestricted | structure | BYM2, unrestricted | 0.176 (0.168-0.184) | 0.96 (0.93-0.98) | NA | −14 | 209 | 51 | 0.80 |
| Structure + composition, spatial unrestricted | structure + composition | BYM2, unrestricted | 0.158 (0.152-0.166) | 0.95 (0.92-0.97) | NA | −2 | 106 | 47 | 0.99 |
| Unrestricted + health region | structure + composition | BYM2, unrestricted, + health region | 0.158 (0.151-0.165) | 0.92 (0.87-0.95) | 0.067 (0.055-0.080) | +13 | 119 | 53 | 0.96 |
| + access to care | structure + composition + access | exchangeable | 0.156 (0.149-0.162) | NA | NA | −185 | 57 | 26 | 0.77 |
| + intermediates | + antenatal visits, preterm, low birth weight | exchangeable | 0.148 (0.141-0.155) | NA | NA | −115 | 46 | 19 | 0.76 |
| + registration indicators | + hospital births, missing schooling, missing race/color | exchangeable | 0.156 (0.150-0.163) | NA | NA | −186 | 56 | 26 | 0.78 |

**Table S2.** Sensitivity covariates (exchangeable models with structure and composition; incidence-rate ratio per standard deviation of the covariate, with 95% credible intervals). The structure and composition coefficients of these models differ from the primary model's by at most 0.02. Distance to intensive care and the presence of beds in the municipality both have ratios above 1, which a simple access explanation would not predict; referral patterns and place of residence may contribute.

| Model | Covariate (SD) | IRR (95% CrI) |
|---|---|---|
| Access | Distance to the nearest municipality with NICU beds, per 100 km (0.61) | 1.010 (1.002-1.018) |
| Access | NICU beds in the municipality, any year (0.25) | 1.013 (1.005-1.021) |
| Access | Births outside the municipality of residence, share (0.34) | 1.001 (0.990-1.012) |
| Access | Births in a municipality with NICU beds that year, share (0.34) | 0.987 (0.976-0.997) |
| Intermediates | Fewer than 7 antenatal visits, share (0.12) | 1.017 (1.007-1.027) |
| Intermediates | Gestational age under 37 weeks, share (0.022) | 1.005 (0.997-1.014) |
| Intermediates | Birth weight under 2,500 g, share (0.016) | 1.077 (1.065-1.089) |
| Registration | Hospital births, share (0.052) | 0.991 (0.984-0.997) |
| Registration | Births with missing maternal schooling, share (0.032) | 1.001 (0.994-1.009) |
| Registration | Births with missing maternal race/color, share (0.076) | 0.991 (0.983-0.999) |

**Table S3.** Minimum residual relative risk detectable with 80% power (one-sided α = 0.05) from sampling variation alone (normal approximation; the exact one-sided Poisson power at 50 expected deaths and a relative risk of 1.35 is about 73%), by expected deaths over the decade, and the number of municipalities at or above each expected count. Adding the overdispersion implied by the residual heterogeneity (variance E + (exp(σ²) − 1)E²) gives 1.53 at 50 expected deaths, 1.43 at 200, and 1.40 at 1,000, approaching 1.40 for large counts.

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

**Table S4.** Registration indicators in municipalities classified above expectation. Medians among the 112 above expectation and among the other 5,453 municipalities, and the 12 of the 112 with less than 90% of births in hospital (national share 98.4%). These indicators describe record quality and do not measure under-registration of births or deaths.

| | Above expectation (n = 112) | Other municipalities (n = 5,453) |
|---|---|---|
| Hospital births, median (%) | 98.7 | 99.4 |
| Births with missing maternal schooling, median (%) | 0.93 | 0.59 |
| Births with missing maternal race/color, median (%) | 1.99 | 1.49 |

| Municipality (state) | Neonatal deaths | Relative risk (posterior median) | Hospital births (%) | Missing schooling (%) | Missing race/color (%) |
|---|---|---|---|---|---|
| Caracaraí (RR) | 156 | 2.25 | 83.2 | 2.2 | 0.9 |
| Barcelos (AM) | 162 | 1.91 | 59.5 | 5.9 | 1.1 |
| Amajari (RR) | 98 | 1.63 | 50.8 | 7.7 | 0.2 |
| Alto Alegre (RR) | 133 | 1.43 | 34.2 | 9.1 | 0.1 |
| Eirunepé (AM) | 162 | 1.42 | 84.9 | 1.9 | 0.2 |
| Envira (AM) | 78 | 1.33 | 88.2 | 0.0 | 0.1 |
| Aveiro (PA) | 33 | 1.32 | 61.0 | 2.2 | 1.1 |
| Uruçuca (BA) | 42 | 1.31 | 84.7 | 3.6 | 4.3 |
| Paratinga (BA) | 66 | 1.31 | 87.3 | 3.5 | 13.9 |
| Nhamundá (AM) | 55 | 1.30 | 86.5 | 19.5 | 4.5 |
| São Gabriel da Cachoeira (AM) | 199 | 1.29 | 54.8 | 0.9 | 0.3 |
| Óbidos (PA) | 150 | 1.25 | 89.7 | 0.9 | 0.6 |

The Spearman correlation of the municipal relative risk with the share of hospital births was -0.06 and with the share of births with missing race/color -0.01.

## S4. Posterior predictive check (primary model)

From 1,000 posterior draws, counts were simulated for every municipality. 1.2% of municipalities fell outside their 95% posterior predictive interval (a low value is expected for a model with one random effect per municipality, so this check has limited power to detect misfit); the share of municipalities with zero deaths was 1.6% observed versus 1.8% predicted; the largest count was 12,293 observed (São Paulo) versus 12,294 predicted (median); the variance-to-mean ratio of the counts was 1,195 observed versus 1,195 predicted.

## S5. Sum of positive excesses under noise

The sum over all municipalities of the positive part of observed minus expected deaths, computed at the posterior median expectations, was 20,279. The model-based sum of the positive part of $E_i(RR_i - 1)$, which shrinks small municipalities toward their expectation, was 15,789 (15,039-16,573) under the primary model. Simulating Poisson counts around the fitted expectations (across posterior draws of the expectations, whereas the observed sum above uses their median) with no residual heterogeneity at all, the same sum of positive parts was 10,923 (10,401-11,462), and 4,539 (4,053-5,075) among municipalities with at least 50 deaths. The positive part of noise is positive, which is why this quantity is not reported as a burden.

## S6. Figure

**Figure S1.** The thirty municipalities with the largest posterior-median excess deaths, with the residual relative risk (median and 95% credible interval). The vertical lines mark 1.00 and 1.10. Classification as above expectation uses the posterior probability P(RR > 1.10) ≥ 0.95, not the interval shown; 35 of the 112 classified municipalities have a 95% interval whose lower bound is at or below 1.10.

## S7. Supplementary data

`municipal_estimates.csv`: for each of the 5,565 municipalities, deaths, births, expected deaths (median and 95% CrI), residual relative risk (median and 95% CrI), posterior probabilities that the relative risk exceeds 1, 1.10, and 1.25, or falls below 0.90, excess deaths (median and 95% CrI), classification, and the three registration indicators.

## S8. Persistence between periods

**Table S5.** Persistence of the residual relative risk between periods. The primary model was fitted separately to each period (period-specific deaths, births, population, and maternal composition; same graph, priors, and sampler settings; the HDI, dated 2010, is the same in both periods).

| | 2014-2019 | 2020-2024 |
|---|---|---|
| Neonatal deaths in the analysis set | 152,418 | 107,547 |
| Live births | 17,568,151 | 12,887,307 |
| SD of log relative risk (95% CrI) | 0.169 (0.161-0.178) | 0.165 (0.156-0.176) |
| Municipalities above expectation | 73 | 80 |
| Municipalities below expectation | 47 | 9 |
| Divergent transitions; maximum R-hat of u | 0; 1.008 | 0; 1.010 |
| Net excess deaths, Brazil (95% CrI) | 602 (−717 to 1,857) | 1,298 (378 to 2,155) |
| Net excess deaths, North (95% CrI) | 465 (128 to 805) | 401 (109 to 673) |
| Net excess deaths, Northeast (95% CrI) | 1,651 (1,210 to 2,130) | 975 (635 to 1,329) |
| Net excess deaths, Central-West (95% CrI) | 60 (−174 to 296) | −31 (−227 to 169) |
| Net excess deaths, Southeast (95% CrI) | −1,263 (−2,067 to −506) | 312 (−232 to 813) |
| Net excess deaths, South (95% CrI) | −314 (−566 to −59) | −359 (−551 to −150) |

Above expectation in both periods: 35; in 2014-2019 only: 38; in 2020-2024 only: 45. Benchmarks for the overlap: 1.0 expected if flags were assigned at random among all municipalities (hypergeometric, p = 3.3 × 10^-48^), and 3.0 expected (maximum 10 in 10,000 permutations; none reached 35, so the Monte Carlo estimate is p = 0.0001) if the 2020-2024 flags were assigned at random within strata of pooled expected deaths (below 20, 20-49, 50-99, 100-199, 200-499, 500 or more). Neither benchmark preserves spatial clustering. Of the 73 above expectation in 2014-2019, 92% had a median relative risk above 1.10 and 99% above 1.00 in 2020-2024; the median relative risk of these municipalities was 1.37 in 2014-2019 and 1.25 in 2020-2024, as expected when municipalities are selected on an estimate. For 31 municipalities, the joint probability of exceeding 1.10 in both periods was at least 0.95 (assuming independent fits on disjoint data). Of the 112 above expectation in the pooled analysis, 35 were above in both periods and 23 in neither. Spearman correlation of the median relative risks between periods: 0.51 (all municipalities), 0.61 (the 783 with at least 25 deaths in each period), 0.58 (the 978 with at least 20 expected deaths in each period); correlation with the pooled estimate 0.88 and 0.80. Deaths are assigned by year of death and births by year of birth, so a death in January 2020 of a child born in December 2019 falls in the second period's numerator and the first period's denominator; the effect was not quantified. The second period includes the COVID-19 pandemic.

## S9. Simulation-based evaluation of the classification rule

Forty datasets were generated with the observed births and covariates, one posterior draw of the primary model as the known truth (intercept, coefficients, and municipal effects) and Poisson counts (39 distinct truths: two replicates drew the same posterior draw); 12 further datasets had no residual (all municipal effects zero). Each was refitted with the primary model (four chains of 2,000 warm-up and 1,000 retained draws; no divergent transitions in the 40 truth datasets) and compared with the truth. The truth is conditional on a generator calibrated to the fitted model and the model is correctly specified, so the calibration in Table S7 is partly expected by construction; these are the results of a simulation-based evaluation under favorable conditions, not a measure of performance on the real data.

In the 40 truth datasets the rule flagged 114.5 municipalities on average (range 86-150) of 1,475 with a true relative risk above 1.10; 1.7% of flagged municipalities (80 of 4,580) had a true relative risk of 1.10 or less (by dataset, 0%-4.8%; Monte Carlo standard errors for the pooled sensitivities were not computed). Overall sensitivity was 7.6% for a true relative risk above 1.10 and 20.4% above 1.25. In the 12 datasets with no residual the rule flagged no municipality above or below expectation. The 95% credible interval covered the true relative risk in 95.0% of municipalities. The net excess (true mean 1,672) was underestimated by 247 on average (15%), with the interval covering the truth in 95% of datasets; the excess among flagged municipalities (true mean 6,065, estimated mean 5,963) was underestimated by 102 (1.7%), with the interval covering the truth in 80% of datasets, below the nominal 95%.

**Table S6.** Sensitivity of the rule and coverage of the 95% credible interval, by true expected deaths of the municipality (pooled over truth datasets).

| Expected deaths | Municipalities with true RR > 1.10 | Flagged | Sensitivity | True RR > 1.25 | Flagged | Sensitivity | 95% interval coverage |
|---|---|---|---|---|---|---|---|
| <20 | 36,576 | 801 | 2.2% | 10,274 | 670 | 6.5% | 95.0% |
| 20-49 | 13,363 | 1,437 | 10.8% | 4,374 | 1,217 | 27.8% | 94.8% |
| 50-199 | 7,323 | 1,385 | 18.9% | 1,955 | 1,041 | 53.2% | 95.1% |
| >=200 | 1,722 | 877 | 50.9% | 609 | 588 | 96.6% | 95.4% |

**Table S7.** Coarse calibration of the exceedance probability P(RR > 1.10): share of municipalities with a true relative risk above 1.10, by bin of estimated probability (pooled over truth datasets; mean predicted probability within bins not computed).

| Estimated P(RR > 1.10) | Municipality-datasets | Share truly above 1.10 |
|---|---|---|
| 0.00 to 0.05 | 45,420 | 1.7% |
| 0.05 to 0.50 | 139,058 | 22.2% |
| 0.50 to 0.80 | 25,885 | 62.4% |
| 0.80 to 0.95 | 7,657 | 87.2% |
| 0.95 to 1.00 (flagged) | 4,580 | 98.3% |

## S10. Comparison with other approaches

**Table S8.** Municipalities flagged above expectation (P(RR > 1.10) ≥ 0.95) by the primary model (restricted spatial regression with structure and maternal composition), by two structure-only models (non-spatial with an exchangeable effect; spatial with the same restricted BYM2 prior), and by crude neonatal mortality rate.

| Comparison | Number |
|---|---|
| Flagged by the primary model | 112 |
| Flagged by the non-spatial structure-only model | 74 |
| Flagged by both | 58 |
| Primary model only | 54 (45 in the Northeast) |
| Non-spatial structure-only model only | 16 |
| Spearman correlation of the municipal relative risks, primary versus non-spatial structure-only | 0.76 |
| Flagged by the spatial structure-only model | 153 |
| Flagged by both the primary and the spatial structure-only model | 93 |
| Primary model only (added by maternal composition) | 19 |
| Spatial structure-only model only (removed by maternal composition) | 60 |
| Flagged by the primary model with at least 50 deaths | 67 |
| ... of which among the 112 highest crude rates (municipalities with at least 50 deaths) | 46 |
| ... of which among the 67 highest crude rates (like-for-like) | 38 |
| Expected number of flagged municipalities with a true relative risk of 1.10 or less (model-based, sum of 1 − P(RR > 1.10)) | 1.8 (1.6%) |
