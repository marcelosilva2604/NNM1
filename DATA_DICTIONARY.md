# Data dictionary / Dicionário de dados

`v2/supplement/municipal_estimates.csv` (supplement section S7): one row per municipality, 5,565 rows,
24 columns. / uma linha por município, 5.565 linhas, 24 colunas.

CrI = credible interval. / ICr = intervalo de credibilidade.

| Column / Coluna | English | Português |
|---|---|---|
| `CODMUNRES` | Six-digit IBGE code of the municipality of residence | Código IBGE de 6 dígitos do município de residência |
| `municipio` | Municipality name | Nome do município |
| `state` | State abbreviation (UF) | Sigla da UF |
| `REG` | Macro-region (North, Northeast, Central-West, Southeast, South) | Grande região (North, Northeast, Central-West, Southeast, South) |
| `deaths` | Neonatal deaths (0 to 27 days), 2014-2024 pooled | Óbitos neonatais (0 a 27 dias), 2014-2024 agrupados |
| `births` | Live births, 2014-2024 pooled | Nascidos vivos, 2014-2024 agrupados |
| `nmr` | Crude neonatal mortality rate per 1,000 live births | Taxa bruta de mortalidade neonatal por 1.000 nascidos vivos |
| `expected_median`, `expected_lo95`, `expected_hi95` | Expected deaths given municipal structure and maternal composition: posterior median and 95% CrI limits | Óbitos esperados dada a estrutura municipal e a composição materna: mediana a posteriori e limites do ICr95% |
| `rr_median`, `rr_lo95`, `rr_hi95` | Residual relative risk, exp(u): posterior median and 95% CrI limits | Risco relativo residual, exp(u): mediana a posteriori e limites do ICr95% |
| `p_rr_gt1` | Posterior probability that the residual RR exceeds 1 | Probabilidade a posteriori de o RR residual superar 1 |
| `p_rr_gt1_10` | Posterior probability that the residual RR exceeds 1.10 | Probabilidade a posteriori de o RR residual superar 1,10 |
| `p_rr_gt1_25` | Posterior probability that the residual RR exceeds 1.25 | Probabilidade a posteriori de o RR residual superar 1,25 |
| `p_rr_lt0_90` | Posterior probability that the residual RR is below 0.90 | Probabilidade a posteriori de o RR residual ficar abaixo de 0,90 |
| `excess_median`, `excess_lo95`, `excess_hi95` | Excess deaths, E x (RR - 1): posterior median and 95% CrI limits (can be negative) | Óbitos em excesso, E x (RR - 1): mediana a posteriori e limites do ICr95% (podem ser negativos) |
| `sh_hospital` | Share of hospital births among births with known place of birth (0 to 1) | Proporção de partos hospitalares entre os nascimentos com local conhecido (0 a 1) |
| `sh_miss_esc` | Share of births with missing maternal schooling (0 to 1) | Proporção de nascimentos com escolaridade materna ausente (0 a 1) |
| `sh_miss_race` | Share of births with missing maternal race/color (0 to 1) | Proporção de nascimentos com raça/cor materna ausente (0 a 1) |
| `class` | `above`: P(RR > 1.10) at least 0.95; `below`: P(RR < 0.90) at least 0.95; `undetermined`: otherwise | `above`: P(RR > 1,10) de pelo menos 0,95; `below`: P(RR < 0,90) de pelo menos 0,95; `undetermined`: nos demais casos |

Notes / Notas:

- The three registration indicators describe record quality; they do not measure under-registration of births or deaths. / Os três indicadores de registro descrevem a qualidade do registro; não medem o sub-registro de nascimentos ou óbitos.
- The classification is a screening aid, not a count of avoidable deaths. / A classificação é um apoio à triagem, não uma contagem de óbitos evitáveis.
- Column definitions follow `v2/src/05_results.py` and `v2/src/03_build_dataset.py`. / As definições das colunas seguem `v2/src/05_results.py` e `v2/src/03_build_dataset.py`.
