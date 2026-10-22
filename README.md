# NNM1: municipalities with neonatal mortality above that expected from socioeconomic structure and maternal composition, Brazil, 2014-2024

[English](#english) | [Português](#português)

---

## English

Analysis code, derived municipal results, and supplementary material for the manuscript
**"Municipalities with neonatal mortality above that expected from socioeconomic structure and
maternal composition: an ecological study with Bayesian spatial analysis, Brazil, 2014-2024"**
(Portuguese title: *Municípios com mortalidade neonatal acima do esperado pela estrutura
socioeconômica e pela composição materna: estudo ecológico com análise espacial bayesiana, Brasil,
2014-2024*), prepared for submission to *Epidemiologia e Serviços de Saúde:
Revista do Sistema Único de Saúde do Brasil* (Epidemiol Serv Saude).

This README is bilingual. Code, code comments, and file names are in English. Everything else
(READMEs, data dictionary, supplementary material) is in English and Portuguese.

### Why this study

Neonatal deaths account for a large share of deaths in children under five, and ending preventable
newborn deaths is a target of Sustainable Development Goal 3.2. Brazil reduced infant mortality for
two decades under a universal health system, but the decline slowed after 2014, and neonatal
mortality still follows the municipal Human Development Index and is higher in the North and
Northeast. That gradient combines two components: the risk profile of the women who give birth in
each place, and the outcome of those births. A model-based expectation for each municipality, with
its uncertainty, describes how far observed mortality departs from what the specified population
profile predicts. The residual has limits: it contains unmeasured population risk, model
misspecification, and registration error, in addition to any deficiency in care, and it is
detectable only where expected counts are large enough (most Brazilian municipalities expect fewer
than 15 neonatal deaths in a decade).

**Objective:** to identify the Brazilian municipalities with neonatal mortality above that expected
from socioeconomic structure and maternal composition, with uncertainty, and to verify whether the
excess persists between 2014-2019 and 2020-2024.

### Relevance for the journal

*Epidemiologia e Serviços de Saúde* is the journal of Brazil's Unified Health System (SUS). The
study uses SUS information systems (SIM, Sinasc, CNES) at the municipal level, and the results are
framed as a surveillance aid: they suggest where confidential investigation of neonatal deaths,
within infant and fetal death surveillance, may concentrate, and where checking the quality of
registration is a reasonable first step. They are not an estimate of avoidable deaths.

### Design and data

- Ecological study of 5,565 municipalities (municipality of maternal residence as the unit),
  2014-2024: 259,965 neonatal deaths (0 to 27 days) and 30,455,458 live births.
- Mortality System (SIM) and Live Birth System (Sinasc), Ministry of Health; CNES (neonatal intensive
  care beds); IBGE population; Atlas of Human Development in Brazil (HDI, 2010).
- Poisson model with a BYM2 spatial random effect under restricted spatial regression (the random
  effect is orthogonal to the covariates), adjusted for municipal structure and for maternal age,
  schooling, and race/color from every birth record. Estimands: residual relative risk (RR),
  exceedance probabilities, and excess deaths, with 95% credible intervals.
- Registration indicators, period-specific fits, alternative specifications, and a simulation-based
  evaluation of the classification rule are in the supplementary material.

### Main findings reported in the manuscript

- Residual heterogeneity: SD of log RR 0.157 (95% CrI 0.151-0.164) in the primary model.
- Net excess over expectation: 1,669 deaths (95% CrI -137 to 3,341). The North and Northeast
  remained above expectation and the South below.
- 112 municipalities were classified above expectation (posterior probability of RR > 1.10 of at
  least 0.95) and 48 below (probability of RR < 0.90 of at least 0.95).
- Maternal composition removed 60 of the 153 municipalities flagged by a spatial model with
  structure only and added 19.
- 35 municipalities were above expectation in both 2014-2019 and 2020-2024.
- In a simulation-based evaluation under a correctly specified model, 1.7% of flagged municipalities
  had a true RR at or below 1.10, and 7.6% of the municipalities with a true RR above 1.10 were
  flagged. This is a simulation under favorable conditions, not a validation on real data.

The classification is a screening aid. A residual can reflect unmeasured risk, misspecification, or
registration error, and the analysis does not exclude registration error.

### Repository map

```
README.md                      this file (English and Portuguese)
DATA_DICTIONARY.md             columns of municipal_estimates.csv (English and Portuguese)
v2/                            current analysis: scripts, municipal aggregates, outputs, supplement
  src/                         numbered scripts 01-10, in logical order
  data/                        municipal aggregates and the adjacency graph (no individual records)
  outputs/                     municipal posterior summaries, tables, figures, run summaries
  supplement/                  full supplementary material (English and Portuguese), Figure S1,
                               municipal_estimates.csv
notebooks/, outputs/           data-preparation notebooks and their results (see below)
requirements.txt               dependencies of the data-preparation notebooks
```

`notebooks/01_cohort_indicators.ipynb` builds the neonatal cohort and the municipal base that
`v2/src/03_build_dataset.py` reads. `notebooks/02_excess_mortality.ipynb` is an auxiliary non-spatial
analysis (standardized mortality ratio, including an avoidable-cause decomposition that is not part
of this manuscript) whose output is read by step 03 for comparison columns (`v1_*`). Neither is the
primary analysis reported in the manuscript.

### Reproducing

Inputs are listed as the scripts read them; the table does not claim that every step runs from a
clean clone without the files marked "not distributed".

| Step | Reads | Distributed here? | Writes |
|---|---|---|---|
| 01 maternal composition | yearly SINASC files (individual level), `cnes_muni_ano.csv` and `nascidos_muni_ano.csv` | no (DATASUS and processed tables; set `NNM1_SINASC`) | `v2/data/maternal_composition*.csv` |
| 02 adjacency graph | `analysis_muni.csv`, IBGE municipal mesh (`br_municipios.geojson`) | table yes; mesh no | `adjacency_pairs.csv`, `municipalities.csv` |
| 03 analysis dataset | processed municipal tables (DATASUS, IBGE, Atlas, CNES), health-region and immediate-region tables (set `NNM1_REF`), notebook outputs | no | `v2/data/analysis_muni.csv` |
| 04 model fitting | `analysis_muni.csv`, `adjacency_pairs.csv`, `bym2_scaling.json` | yes | `v2/model/*.nc` (about 2.7 GB, not versioned), `muni_*.csv` |
| 05 results | model output, `analysis_muni.csv`, municipal and state meshes (`br_municipios.geojson`, `br_states.geojson`) | table yes; meshes no; model output from 04 | `v2/outputs/` tables, figures, `NUMBERS.json` |
| 06 manuscript tables | `v2/outputs/tables/` | yes | `manuscript_tables.md` |
| 07 period persistence | annual deaths table (`deaths_muni_year.csv`, distributed; the death cohort only if it is missing), processed births and population tables | deaths table yes; processed tables no | `period_estimates.csv`, `period_persistence.json` |
| 08 simulation | `model/idata_bym2_rsr_full.nc` (from 04) | from 04 | `v2/outputs/simulation/`, calibration tables |
| 09 comparison with routine outputs | model summaries (04), `municipal_estimates.csv` (05) | from 04 and 05 | `added_value.json` |
| 10 Portuguese figures | `v2/supplement/municipal_estimates.csv`, IBGE municipal and state meshes | CSV yes; meshes no | `v2/outputs/figures_pt/` |

Environment: see `v2/requirements.txt` (PyMC 6.2, geopandas, libpysal; steps 02-09) and polars (step 01).
The provided results in `v2/outputs/` are those generated by the authors' run.
`v2/outputs/tables/manuscript_tables.md` was written for an earlier table numbering and includes the
exploratory cause analysis; the final tables are in the manuscript and in the supplementary material.
The figures in `outputs/figures/` (repository root) belong to the auxiliary notebooks. The figures of the
manuscript are those in `v2/outputs/figures_pt/` (Portuguese, with scale bar and north arrow; script
`v2/src/10_figures_pt.py`); `v2/outputs/figures/` holds the earlier English versions, without scale bar and
north arrow. Supplementary Figure S1 is in `v2/supplement/`.
The steps are numbered in logical order; step 02 reads the table written by step 03, so to start at step 04
use the distributed `analysis_muni.csv` and `adjacency_pairs.csv`.

### Data sources (all public)

- SIM and Sinasc, Ministry of Health, [DATASUS](https://datasus.saude.gov.br/transferencia-de-arquivos/)
- CNES (neonatal intensive care beds), DATASUS
- Municipal population, [IBGE/SIDRA](https://sidra.ibge.gov.br/)
- Municipal Human Development Index (2010), [Atlas Brasil](https://www.atlasbrasil.org.br/)

Individual-level files are not distributed here. The repository holds municipal aggregates only.

### Supplementary material

`v2/supplement/`: `Supplementary_Material_EN.pdf` and `Supplementary_Material_PT.pdf` (sections
S1-S10, Tables S1-S8, A4; Markdown sources `supplement_en.md` and `supplement_pt.md`), `Figure_S1.pdf`
(English) and `Figure_S1_PT.pdf` (Portuguese), and `municipal_estimates.csv` (S7). Column definitions are in `DATA_DICTIONARY.md`.

---

## Português

Código de análise, resultados municipais derivados e material suplementar do manuscrito
**"Municípios com mortalidade neonatal acima do esperado pela estrutura socioeconômica e pela
composição materna: estudo ecológico com análise espacial bayesiana, Brasil, 2014-2024"** (título em
inglês: *Municipalities with neonatal mortality above that expected from socioeconomic structure and
maternal composition: an ecological study with Bayesian spatial analysis, Brazil, 2014-2024*),
preparado para submissão à revista *Epidemiologia e Serviços de Saúde:
Revista do Sistema Único de Saúde do Brasil* (Epidemiol Serv Saude).

Este README é bilíngue. O código, os comentários do código e os nomes de arquivo estão em inglês.
O restante (READMEs, dicionário de dados, material suplementar) está em inglês e em português.

### Justificativa

Os óbitos neonatais respondem por grande parcela das mortes em menores de cinco anos, e o fim das
mortes evitáveis de recém-nascidos é meta do Objetivo de Desenvolvimento Sustentável 3.2. O Brasil,
sob um sistema universal de saúde, reduziu a mortalidade infantil por duas décadas, mas a queda
desacelerou após 2014, e a mortalidade neonatal ainda acompanha o Índice de Desenvolvimento Humano
Municipal e é mais alta no Norte e no Nordeste. Esse gradiente combina dois componentes: o perfil de
risco das mulheres que dão à luz em cada lugar e o desfecho desses nascimentos. Uma expectativa
explícita, baseada em modelo, para cada município, com sua incerteza, descreve o quanto a
mortalidade observada se afasta do que o perfil populacional especificado prediz. O resíduo tem
limites: contém risco populacional não medido, erro de especificação do modelo e erro de registro,
além de qualquer deficiência do cuidado, e só é detectável onde as contagens esperadas são grandes
o bastante (na maioria dos municípios brasileiros esperam-se menos de 15 óbitos neonatais em uma
década).

**Objetivo:** identificar os municípios brasileiros com mortalidade neonatal acima do esperado pela
estrutura socioeconômica e pela composição materna, com incerteza, e verificar se o excesso persiste
entre 2014-2019 e 2020-2024.

### Pertinência para a revista

A *Epidemiologia e Serviços de Saúde* é a revista do Sistema Único de Saúde (SUS) do Brasil. O
estudo usa sistemas de informação do SUS (SIM, Sinasc, CNES) em nível municipal, e os resultados
são apresentados como apoio à vigilância: sugerem onde a investigação confidencial dos óbitos
neonatais, no âmbito da vigilância do óbito infantil e fetal, pode se concentrar, e onde verificar
a qualidade do registro é um primeiro passo razoável. Não são uma estimativa de óbitos evitáveis.

### Desenho e dados

- Estudo ecológico de 5.565 municípios (município de residência da mãe como unidade), 2014-2024:
  259.965 óbitos neonatais (0 a 27 dias) e 30.455.458 nascidos vivos.
- Sistema de Informações sobre Mortalidade (SIM) e Sistema de Informações sobre Nascidos Vivos
  (Sinasc), Ministério da Saúde; CNES (leitos de terapia intensiva neonatal); população do IBGE;
  Atlas do Desenvolvimento Humano no Brasil (IDHM, 2010).
- Modelo de Poisson com efeito aleatório espacial BYM2 sob regressão espacial restrita (o efeito
  aleatório é ortogonal às covariáveis), ajustado pela estrutura municipal e pela idade, escolaridade
  e raça/cor das mães, obtidas de todos os registros de nascimento. Estimandos: risco relativo (RR)
  residual, probabilidades de excedência e óbitos em excesso, com intervalos de credibilidade de 95%.
- Indicadores de registro, ajustes por período, especificações alternativas e uma avaliação por
  simulação da regra de classificação estão no material suplementar.

### Principais resultados descritos no manuscrito

- Heterogeneidade residual: DP do log do RR de 0,157 (ICr95% 0,151; 0,164) no modelo primário.
- Excesso líquido sobre o esperado: 1.669 óbitos (ICr95% -137; 3.341). Nordeste e Norte
  permaneceram acima do esperado, e o Sul, abaixo.
- 112 municípios foram classificados acima do esperado (probabilidade a posteriori de RR > 1,10 de
  pelo menos 0,95) e 48 abaixo (probabilidade de RR < 0,90 de pelo menos 0,95).
- A composição materna removeu 60 dos 153 municípios sinalizados por um modelo espacial apenas com
  estrutura e incluiu 19.
- 35 municípios estiveram acima do esperado em 2014-2019 e em 2020-2024.
- Em uma avaliação por simulação com o modelo corretamente especificado, 1,7% dos municípios
  sinalizados tinham RR verdadeiro igual ou inferior a 1,10, e 7,6% dos municípios com RR verdadeiro
  acima de 1,10 foram sinalizados. É uma simulação em condições favoráveis, não uma validação em
  dados reais.

A classificação é um apoio à triagem. O resíduo pode refletir risco não medido, erro de especificação
ou erro de registro, e a análise não exclui erro de registro.

### Mapa do repositório

```
README.md                      este arquivo (inglês e português)
DATA_DICTIONARY.md             colunas de municipal_estimates.csv (inglês e português)
v2/                            análise atual: scripts, agregados municipais, saídas, suplemento
  src/                         scripts numerados de 01 a 10, em ordem lógica
  data/                        agregados municipais e grafo de adjacências (sem registros individuais)
  outputs/                     resumos a posteriori por município, tabelas, figuras, resumos das execuções
  supplement/                  material suplementar completo (inglês e português), Figura S1,
                               municipal_estimates.csv
notebooks/, outputs/           notebooks de preparação dos dados e seus resultados (ver abaixo)
requirements.txt               dependências dos notebooks de preparação dos dados
```

O `notebooks/01_cohort_indicators.ipynb` constrói a coorte neonatal e a base municipal que o
`v2/src/03_build_dataset.py` lê. O `notebooks/02_excess_mortality.ipynb` é uma análise não espacial
auxiliar (razão de mortalidade padronizada, com uma decomposição por causas evitáveis que não faz
parte deste manuscrito), cuja saída é lida pelo passo 03 para colunas de comparação (`v1_*`).
Nenhum dos dois é a análise primária descrita no manuscrito.

### Reprodução

As entradas estão listadas como os scripts as leem; a tabela não afirma que todo passo roda a partir
de um clone limpo sem os arquivos marcados como "não distribuído".

| Passo | Lê | Distribuído aqui? | Escreve |
|---|---|---|---|
| 01 composição materna | arquivos anuais do Sinasc (nível individual), `cnes_muni_ano.csv` e `nascidos_muni_ano.csv` | não (DATASUS e tabelas processadas; definir `NNM1_SINASC`) | `v2/data/maternal_composition*.csv` |
| 02 grafo de adjacências | `analysis_muni.csv`, malha municipal do IBGE (`br_municipios.geojson`) | tabela sim; malha não | `adjacency_pairs.csv`, `municipalities.csv` |
| 03 base de análise | tabelas municipais processadas (DATASUS, IBGE, Atlas, CNES), tabelas de região de saúde e região imediata (definir `NNM1_REF`), saídas dos notebooks | não | `v2/data/analysis_muni.csv` |
| 04 ajuste dos modelos | `analysis_muni.csv`, `adjacency_pairs.csv`, `bym2_scaling.json` | sim | `v2/model/*.nc` (cerca de 2,7 GB, não versionado), `muni_*.csv` |
| 05 resultados | saída do modelo, `analysis_muni.csv`, malhas municipal e estadual (`br_municipios.geojson`, `br_states.geojson`) | tabela sim; malhas não; saída do modelo vem do 04 | tabelas, figuras e `NUMBERS.json` em `v2/outputs/` |
| 06 tabelas do manuscrito | `v2/outputs/tables/` | sim | `manuscript_tables.md` |
| 07 persistência entre períodos | tabela anual de óbitos (`deaths_muni_year.csv`, distribuída; a coorte de óbitos só se ela faltar), tabelas processadas de nascimentos e população | tabela de óbitos sim; tabelas processadas não | `period_estimates.csv`, `period_persistence.json` |
| 08 simulação | `model/idata_bym2_rsr_full.nc` (do 04) | vem do 04 | `v2/outputs/simulation/`, tabelas de calibração |
| 09 comparação com saídas de rotina | resumos dos modelos (04), `municipal_estimates.csv` (05) | vêm do 04 e do 05 | `added_value.json` |
| 10 figuras em português | `v2/supplement/municipal_estimates.csv`, malhas municipal e estadual do IBGE | CSV sim; malhas não | `v2/outputs/figures_pt/` |

Ambiente: ver `v2/requirements.txt` (PyMC 6.2, geopandas, libpysal; passos 02 a 09) e polars (passo 01).
Os resultados fornecidos em `v2/outputs/` são os gerados pela execução dos autores.
`v2/outputs/tables/manuscript_tables.md` foi escrito para uma numeração anterior das tabelas e inclui a
análise exploratória por causa; as tabelas finais estão no manuscrito e no material suplementar.
As figuras de `outputs/figures/` (raiz do repositório) pertencem aos notebooks auxiliares. As figuras do
manuscrito são as de `v2/outputs/figures_pt/` (em português, com barra de escala e seta do Norte; script
`v2/src/10_figures_pt.py`); `v2/outputs/figures/` traz as versões anteriores em inglês, sem barra de escala
e seta do Norte. A Figura S1 do suplemento está em `v2/supplement/`.
Os passos estão numerados em ordem lógica; o passo 02 lê a tabela escrita pelo passo 03, de modo que, para
começar no passo 04, use o `analysis_muni.csv` e o `adjacency_pairs.csv` distribuídos.

### Fontes de dados (todas públicas)

- SIM e Sinasc, Ministério da Saúde, [DATASUS](https://datasus.saude.gov.br/transferencia-de-arquivos/)
- CNES (leitos de terapia intensiva neonatal), DATASUS
- População municipal, [IBGE/SIDRA](https://sidra.ibge.gov.br/)
- Índice de Desenvolvimento Humano Municipal (2010), [Atlas Brasil](https://www.atlasbrasil.org.br/)

Os arquivos de nível individual não são distribuídos aqui. O repositório contém apenas agregados
municipais.

### Material suplementar

`v2/supplement/`: `Supplementary_Material_EN.pdf` e `Supplementary_Material_PT.pdf` (seções S1 a S10,
Tabelas S1 a S8, A4; fontes em Markdown `supplement_en.md` e `supplement_pt.md`), `Figure_S1.pdf`
(inglês) e `Figure_S1_PT.pdf` (português), e `municipal_estimates.csv` (S7). As definições das colunas estão em `DATA_DICTIONARY.md`.
