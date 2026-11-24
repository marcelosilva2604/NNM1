**eTable 1.** Share of municipality-replicates with true RR > 1.10 by bin of estimated P(RR > 1.10), with the number of municipality-replicates in each bin.

| Estimated probability | k = 0.5 | k = 1.0 | k = 1.5 | k = 2.0 |
|---|---|---|---|---|
| 0.00-0.05 | 2.3% (n = 238,387) | 1.8% (n = 112,814) | 1.4% (n = 110,869) | 1.2% (n = 120,843) |
| 0.05-0.50 | 14.8% (n = 302,434) | 22.2% (n = 349,029) | 23.8% (n = 293,001) | 24.0% (n = 251,209) |
| 0.50-0.80 | 62.7% (n = 13,030) | 62.6% (n = 64,818) | 63.2% (n = 87,894) | 64.0% (n = 90,372) |
| 0.80-0.95 | 88.0% (n = 2,137) | 87.2% (n = 18,578) | 87.7% (n = 33,989) | 87.9% (n = 42,058) |
| 0.95-1.00 | 99.4% (n = 512) | 98.4% (n = 11,261) | 98.7% (n = 30,747) | 98.9% (n = 52,018) |

**eTable 2.** Rules and comparators by k, shown as sensitivity / false-flag proportion / flagged municipality-replicates (pooled; target true RR above the rule's own tau, which is 1.10 except in the first and sixth rows). Under k = 0, C2 flagged 3,159 municipality-replicates in 40 replicates; the other rules flagged none.

| Rule | k = 0.5 | k = 1.0 | k = 1.5 | k = 2.0 |
|---|---|---|---|---|
| P(RR > 1.05) >= 0.95 | 2.0% / 2.3% / 2,909 | 9.7% / 1.6% / 20,448 | 18.6% / 1.3% / 43,202 | 26.9% / 1.2% / 65,495 |
| P(RR > 1.10) >= 0.99 | 0.3% / 0.0% / 193 | 3.7% / 0.1% / 5,398 | 9.3% / 0.2% / 17,366 | 16.0% / 0.2% / 33,269 |
| P(RR > 1.10) >= 0.95 (primary) | 0.8% / 0.6% / 512 | 7.5% / 1.6% / 11,261 | 16.2% / 1.3% / 30,747 | 24.7% / 1.1% / 52,018 |
| P(RR > 1.10) >= 0.90 | 1.6% / 3.3% / 1,008 | 11.3% / 3.7% / 17,291 | 22.1% / 3.0% / 42,572 | 31.6% / 2.6% / 67,571 |
| P(RR > 1.10) >= 0.80 | 3.9% / 9.8% / 2,649 | 18.5% / 8.6% / 29,839 | 32.2% / 7.1% / 64,736 | 42.5% / 6.0% / 94,076 |
| P(RR > 1.20) >= 0.95 | 0.1% / 0.0% / 9 | 5.2% / 1.2% / 3,535 | 13.0% / 1.4% / 15,644 | 21.6% / 1.2% / 33,294 |
| C1: posterior median RR > 1.10 | 17.3% / 32.6% / 15,670 | 46.1% / 28.3% / 94,613 | 61.8% / 24.2% / 152,571 | 70.3% / 20.7% / 184,401 |
| C2: Poisson tail test, p < 0.05 | 7.6% / 61.3% / 11,900 | 14.0% / 19.1% / 25,398 | 21.8% / 7.7% / 44,245 | 29.4% / 4.4% / 63,990 |
| C3: Poisson tail test, BH 5% | 0.1% / 0.0% / 33 | 1.9% / 0.6% / 2,772 | 6.3% / 0.4% / 11,752 | 13.3% / 0.3% / 27,758 |

**eTable 3.** Pooled coverage of the 95% credible interval for the residual RR, by generating expected deaths, and among flagged municipality-replicates under the primary rule.

| Stratum or subset | k = 0 | k = 0.5 | k = 1.0 | k = 1.5 | k = 2.0 |
|---|---|---|---|---|---|
| <20 | 100.0% | 94.6% | 94.9% | 95.0% | 95.1% |
| 20-49 | 100.0% | 94.6% | 95.0% | 95.0% | 95.0% |
| 50-199 | 100.0% | 94.7% | 95.1% | 95.1% | 95.1% |
| >=200 | 100.0% | 94.8% | 94.4% | 95.0% | 94.9% |
| Among flagged, pooled | NA | 95.1% | 93.1% | 94.5% | 94.8% |
| Among flagged, mean of replicate-level values | NA | 94.5% | 93.0% | 94.5% | 94.8% |

**eTable 4.** Convergence diagnostics by scenario. The hyperparameter and municipal-effect R-hat were recorded only for replicates run after the plan amendment.

| k | Replicates | With divergences | Divergences | Hyperparameter R-hat > 1.10 | Max hyperparameter R-hat | Municipal-effect R-hat > 1.10 | Max municipal-effect R-hat | Overall R-hat > 1.10 | Max overall R-hat |
|---|---|---|---|---|---|---|---|---|---|
| 0.0 | 40 | 0 | 0 | 2 of 40 | 1.134 | 24 of 40 | 1.247 | 40 of 40 | 2.374 |
| 0.5 | 100 | 2 | 15 | 3 of 93 | 1.125 | 0 of 93 | 1.057 | 12 of 100 | 1.165 |
| 1.0 | 100 | 1 | 1 | 0 of 93 | 1.031 | 0 of 93 | 1.010 | 0 of 100 | 1.031 |
| 1.5 | 100 | 0 | 0 | 0 of 94 | 1.050 | 0 of 94 | 1.011 | 0 of 100 | 1.050 |
| 2.0 | 100 | 0 | 0 | 0 of 94 | 1.094 | 0 of 94 | 1.011 | 0 of 100 | 1.094 |

**eTable 5.** Standard deviation of the generating log RR across municipalities (mean and range over replicates) and share of municipality-replicates with true RR > 1.10.

| k | SD of generating log RR | True RR > 1.10 |
|---|---|---|
| 0.5 | 0.078 (0.075-0.082) | 10.9% |
| 1.0 | 0.157 (0.147-0.166) | 26.5% |
| 1.5 | 0.236 (0.224-0.253) | 33.6% |
| 2.0 | 0.315 (0.300-0.328) | 37.4% |
