# Analysis plan for the evaluation of the classification rule (written before the simulation was run, 23 November 2026)

English translation of the plan as it was written, with the amendment of 24 November 2026. Two phrases that referred to files outside this repository were generalized.

Question: what are the operating characteristics (sensitivity, proportion of false flags, coverage, calibration) of the rule P(RR > tau) >= c in a BYM2 model with restricted spatial regression under sparse counts?

## Generator (known)
- Primary fitted model (`primary_fit/model/idata_bym2_rsr_full.nc`), one posterior extraction per replicate: b0, beta, u. Observed births and covariates.
- Heterogeneity scenarios: true u = k * extracted u, with k in {0, 0.5, 1.0, 1.5, 2.0} (k = 0 is the null; k = 1 repeats the fitted heterogeneity). This preserves the spatial structure and changes the scale of the heterogeneity.
- 100 replicates per scenario (500 fits). Seeds: base 20261123 + 1000 * scenario index + replicate. Poisson counts with mean E * RR.
- Each replicate refits the primary model (4 chains, 1000 retained draws, 2000 tuning draws, target acceptance 0.9, as in the earlier simulation script, which is not included).

## Rules evaluated (no new fit; computed from the same posteriors)
- Family: P(RR > tau) >= c, with tau in {1.05, 1.10, 1.20} and c in {0.80, 0.90, 0.95, 0.99}. True target: true RR > tau (same tau).
- Comparators (target true RR > 1.10):
  - C1: posterior median RR > 1.10.
  - C2: one-sided exact Poisson tail test of H0: mean = 1.10 * E_hat, p < 0.05 without adjustment (E_hat = births * exp(b0 + X beta), posterior mean).
  - C3: the same test with Benjamini-Hochberg adjustment at 5% across the 5,565 municipalities.
- Primary rule of the evaluation: tau = 1.10, c = 0.95 (the rule used in the registry application).

## Metrics, by scenario, rule, and stratum of true E (<20, 20-49, 50-199, >=200)
- Sensitivity = flagged with true target / municipalities with true target (pooled across replicates and as the mean of replicate-level values with Monte Carlo standard error = SD / square root of 100).
- Proportion of false flags = flagged without true target / flagged.
- Null: number flagged per replicate and replicates with any flag.
- Coverage of the 95% interval (all, and among flagged only); calibration of P(RR > 1.10) in bins (0-0.05, 0.05-0.5, 0.5-0.8, 0.8-0.95, 0.95-1).
- Convergence: divergences and maximum R-hat per replicate.

## Limits declared in advance
- The generator is calibrated on the fitted model and the model is correctly specified: this shows performance under favorable conditions, not in the registry data.
- Heterogeneity scaled by k is a scenario, not a truth about the world; priors and generator do not vary.
- No simulated result is converted into a count of correct classifications in the registry data.
- Sensitivity to the neighborhood matrix is done separately, on the registry data (step 2), not in the simulation.

## Amendment of 24 November 2026 (made after seeing 6 null replicates with the original sampling; the other scenarios did not influence the amendment)
- In the first 6 null replicates (k = 0, 1000 retained and 2000 tuning draws), maximum R-hat was 1.67 to 2.22 (no municipality flagged), that is, weak mixing of the municipal effects when the generator has no heterogeneity.
- Decision: null scenario with 4000 tuning and 4000 retained draws (same seeds), limited to 40 replicates (not 100) because of cost; the 6 initial JSON files are kept in `simulation/outputs/grid_null_initial/` and replicates 0 to 5 were redone. A seventh initial null replicate (R-hat 1.74) was identified later and treated in the same way.
- R-hat is recorded separately for the hyperparameters (b0, beta, sigma_u, rho) and for the municipal effects u. Replicates run before the amendment do not have the field; their overall R-hat is in the JSON.
- Any replicate with hyperparameter R-hat above 1.10 is reported as such and not discarded.
