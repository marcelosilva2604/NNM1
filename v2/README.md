# Analysis code / Código de análise (v2)

**English.** Bayesian spatial analysis of neonatal mortality beyond socioeconomic structure and
maternal composition across Brazilian municipalities, 2014-2024 (BYM2 with restricted spatial
regression, Poisson counts, 5,565 municipalities). This folder is the primary analysis of the
manuscript prepared for submission to *Epidemiologia e Serviços de Saúde*. The repository README
(`../README.md`) gives the rationale, design, main findings, and the reproduction table.

**Português.** Análise espacial bayesiana da mortalidade neonatal além da estrutura socioeconômica e
da composição materna nos municípios brasileiros, 2014-2024 (BYM2 com regressão espacial restrita,
contagens de Poisson, 5.565 municípios). Esta pasta é a análise primária do manuscrito preparado
para submissão à *Epidemiologia e Serviços de Saúde*. O README do repositório (`../README.md`)
traz a justificativa, o desenho, os principais resultados e a tabela de reprodução.

```
src/        numbered scripts, in order of execution / scripts numerados, em ordem de execução
              01 maternal composition and place of birth (from SINASC) / composição materna e local do parto
              02 adjacency graph / grafo de adjacências
              03 analysis dataset / base de análise
              04 model fitting (PyMC) / ajuste dos modelos
              05 results and figures / resultados e figuras
              06 manuscript tables / tabelas do manuscrito
              07 persistence between 2014-2019 and 2020-2024 / persistência entre os períodos
              08 simulation-based calibration of the classification rule / calibração por simulação da regra de classificação
              09 comparison with routine outputs / comparação com saídas de rotina
data/       municipal aggregates and the adjacency graph (no individual records) /
              agregados municipais e grafo de adjacências (sem registros individuais)
outputs/    municipal posterior summaries, tables, figures (PNG) and run summaries /
              resumos a posteriori por município, tabelas, figuras (PNG) e resumos das execuções
supplement/ full supplementary material in English and Portuguese (Markdown and PDF), Figure S1
              (English and Portuguese), and municipal_estimates.csv (S7) / material suplementar completo
              em inglês e português (Markdown e PDF), Figura S1 (inglês e português) e municipal_estimates.csv (S7)
```

`outputs/tables/cause_negative_control.csv` is the output of an exploratory analysis that is not
reported in the manuscript. / `outputs/tables/cause_negative_control.csv` é a saída de uma análise
exploratória que não consta no manuscrito.

Posterior samples (`model/*.nc`, about 2.7 GB) and the raw DATASUS files are not versioned. Raw
sources: SIM, SINASC and CNES (DATASUS), IBGE and Atlas Brasil. The simulation study evaluates the
method under a correctly specified model and is not a test on the real data. / As amostras a
posteriori (`model/*.nc`, cerca de 2,7 GB) e os arquivos brutos do DATASUS não são versionados. O
estudo de simulação avalia o método com o modelo corretamente especificado e não é um teste nos
dados reais.

Environment / Ambiente: see `requirements.txt` / ver `requirements.txt`.
