# Material suplementar

**Mortalidade neonatal além da estrutura socioeconômica e da composição materna: estudo ecológico com análise espacial bayesiana nos municípios brasileiros, 2014-2024**

*As Tabelas S1 e S6 também acompanham o manuscrito, como Tabela suplementar 1 e Tabela suplementar 2.*

## S1. Grafo de adjacências e escalonamento do BYM2

A contiguidade do tipo rainha foi calculada sobre a malha municipal do IBGE (46 geometrias inválidas corrigidas antes do uso) para os 5.565 municípios analisados, o que resultou em 16.129 pares de vizinhos não direcionados (média de 5,8 vizinhos; mínimo 1, máximo 23). Quatro municípios não tinham vizinho contíguo (polígonos oceânicos ou costeiros que não se tocam) e o grafo tinha vários componentes; cada componente isolado foi ligado ao município mais próximo pela distância entre centroides (projeção Policônica do Brasil, EPSG:5880) até o grafo ficar conexo, pois a distribuição a priori autorregressiva condicional intrínseca é imprópria em cada componente desconexo. Nazária (PI), emancipado em 2008 e ausente da malha, foi posicionado nas coordenadas medianas de seus oito estabelecimentos do CNES e ligado aos seus três vizinhos mais próximos. Foram 10 ligações artificiais, com comprimento mediano de 19,1 km; a maior foi de 374,8 km. Restrições e escalonamento específicos por componente são uma alternativa que não foi ajustada. O fator de escala do BYM2, a média geométrica da diagonal da inversa generalizada do laplaciano do grafo (Riebler et al., 2016), foi 3,58.

## S2. Especificação do modelo, distribuições a priori e diagnósticos

Para o município $i$ com $b_i$ nascidos vivos, covariáveis padronizadas $x_i$ e $y_i$ óbitos neonatais:

$$y_i \sim \mathrm{Poisson}\big(b_i \exp(\beta_0 + x_i'\beta + u_i)\big)$$

$$u^{*} = \sigma_u \big( (\rho/s)^{1/2} \phi + (1-\rho)^{1/2} \theta \big)$$

em que $\phi$ é um campo autorregressivo condicional intrínseco em $W$ (ICAR) e $\theta_i \sim N(0, 1)$ é permutável, e

$$u = \big(I - X_a (X_a' X_a)^{-1} X_a'\big)\, u^{*}, \qquad X_a = [\mathbf{1}, X]$$

sendo $s$ o fator de escala e $W$ a matriz de adjacências. A projeção na última linha é a regressão espacial restrita (Hodges & Reich, 2010; Hughes & Haran, 2013): o efeito aleatório fica confinado ao complemento ortogonal do espaço das covariáveis, de modo que $u$ é ortogonal a esse espaço e corresponde ao que as covariáveis não predizem; centrar $u$ na escala logarítmica não centra $\exp(u)$, de modo que o excesso líquido nacional não precisa ser zero. Sem a projeção, o campo espacial suave absorveu as covariáveis (RTI da renda por desvio-padrão 0,955 apenas com a estrutura e 0,976 com a composição, contra 0,913 e 0,983 sob a projeção; ρ 0,95), ao mesmo tempo que produziu um mapa de resíduos semelhante (Spearman 0,99 com o ajuste restrito; 106 contra 112 municípios acima do esperado) e um ajuste semelhante (Tabela S1, Δelpd −2); a projeção altera a interpretação dos coeficientes mais do que o mapa. Khan & Calder (2022) mostram que os métodos restritos podem subestimar a incerteza dos coeficientes das covariáveis; os coeficientes não são o estimando principal aqui.

Distribuições a priori: $\beta_0 \sim N(\log(\bar y / \bar b), 1)$, $\beta_j \sim N(0, 0{,}5)$ na escala padronizada, $\sigma_u \sim$ meia-normal(0,5), $\rho \sim$ Beta(1, 1) e, no modelo com região de saúde, $\gamma_r = \sigma_{hr} z_r$ com $z_r \sim N(0,1)$ e $\sigma_{hr} \sim$ meia-normal(0,5), somado a $u^{*}$ antes da projeção. Os modelos permutáveis substituem a convolução por $\sigma_u \theta$.

Amostragem: PyMC 6.2 (Abril-Pla et al., 2023), amostrador No-U-Turn, taxa de aceitação-alvo de 0,9, quatro cadeias com 2.000 iterações de aquecimento e 1.000 amostras retidas (4.000 amostras). Modelo primário: 0 transições divergentes; $\hat R$ dividido máximo de 1,008 (hiperparâmetros 1,007); tamanho efetivo de amostra (bulk ESS) mínimo de 346 para os hiperparâmetros e 1.807 para os efeitos municipais. Todos os onze ajustes do período completo tiveram 0 divergências e $\hat R$ máximo ≤ 1,055 (o BYM2 irrestrito com composição, cujos coeficientes das covariáveis misturam lentamente pelo motivo acima; seu mapa de resíduos é semelhante ao do ajuste restrito, Spearman 0,99; $\hat R$ máximo dos efeitos municipais 1,02).

A distância ao município mais próximo com leitos de terapia intensiva neonatal foi calculada a partir dos estabelecimentos registrados no último ano em que cada município aparece no CNES; a presença de leitos usa qualquer ano de 2014-2024. A heterogeneidade é informada como o desvio-padrão de $u$ entre os municípios dentro de uma amostra a posteriori (mediana das amostras), porque, sob um ICAR escalonado, o parâmetro de escala $\sigma_u$ do BYM2 (0,37 no modelo primário) não é a dispersão realizada do campo.

## S3. Tabelas

**Tabela S1.** Todos os modelos ajustados: heterogeneidade residual, parcela espacial, desvio-padrão da região de saúde, ajuste e contagens de classificação. Δelpd é a diferença na densidade preditiva logarítmica esperada (WAIC) em relação ao modelo primário; valores positivos favorecem o modelo da linha. As diferenças são mostradas sem erros-padrão, e diferenças de poucas unidades não são interpretadas. A parcela espacial ρ e a escala da região de saúde são parâmetros da distribuição a priori definida antes da projeção, não parcelas da variância realizada do campo projetado. NA, não aplicável.

| Modelo | Covariáveis | Distribuição a priori de u | DP do log RR (ICr95%) | ρ | Efeito de região de saúde, DP na escala logarítmica | Δelpd em relação ao primário | Acima | Abaixo | Spearman com o primário |
|---|---|---|---|---|---|---|---|---|---|
| Estrutura, permutável | estrutura | permutável | 0,169 (0,162; 0,176) | NA | NA | −229 | 74 | 32 | 0,76 |
| Estrutura + composição, permutável | estrutura + composição | permutável | 0,157 (0,150; 0,164) | NA | NA | −187 | 59 | 26 | 0,78 |
| Estrutura, espacial restrito | estrutura | BYM2, restrito | 0,169 (0,162; 0,175) | 0,96 (0,93; 0,97) | NA | −15 | 153 | 57 | 0,89 |
| **Primário** | estrutura + composição | BYM2, restrito | 0,157 (0,151; 0,164) | 0,95 (0,93; 0,97) | NA | 0 | 112 | 48 | 1,00 |
| Primário + região de saúde | estrutura + composição | BYM2, restrito, + região de saúde | 0,156 (0,150; 0,163) | 0,92 (0,87; 0,95) | 0,068 (0,055; 0,080) | +13 | 128 | 49 | 0,97 |
| Estrutura, espacial irrestrito | estrutura | BYM2, irrestrito | 0,176 (0,168; 0,184) | 0,96 (0,93; 0,98) | NA | −14 | 209 | 51 | 0,80 |
| Estrutura + composição, espacial irrestrito | estrutura + composição | BYM2, irrestrito | 0,158 (0,152; 0,166) | 0,95 (0,92; 0,97) | NA | −2 | 106 | 47 | 0,99 |
| Irrestrito + região de saúde | estrutura + composição | BYM2, irrestrito, + região de saúde | 0,158 (0,151; 0,165) | 0,92 (0,87; 0,95) | 0,067 (0,055; 0,080) | +13 | 119 | 53 | 0,96 |
| + acesso à assistência | estrutura + composição + acesso | permutável | 0,156 (0,149; 0,162) | NA | NA | −185 | 57 | 26 | 0,77 |
| + intermediários | + consultas de pré-natal, prematuridade, baixo peso ao nascer | permutável | 0,148 (0,141; 0,155) | NA | NA | −115 | 46 | 19 | 0,76 |
| + indicadores de registro | + partos hospitalares, escolaridade ausente, raça/cor ausente | permutável | 0,156 (0,150; 0,163) | NA | NA | −186 | 56 | 26 | 0,78 |

**Tabela S2.** Covariáveis de sensibilidade (modelos permutáveis com estrutura e composição; razão de taxas de incidência (RTI) por desvio-padrão da covariável, com intervalos de credibilidade de 95%). Os coeficientes de estrutura e composição desses modelos diferem dos do modelo primário em no máximo 0,02. A distância até a terapia intensiva e a presença de leitos no município têm razões acima de 1, o que uma explicação simples por acesso não prediria; padrões de referenciamento e o local de residência podem contribuir.

| Modelo | Covariável (DP) | RTI (ICr95%) |
|---|---|---|
| Acesso | Distância ao município mais próximo com leitos de UTI neonatal, por 100 km (0,61) | 1,010 (1,002; 1,018) |
| Acesso | Leitos de UTI neonatal no município, qualquer ano (0,25) | 1,013 (1,005; 1,021) |
| Acesso | Partos fora do município de residência, proporção (0,34) | 1,001 (0,990; 1,012) |
| Acesso | Partos em município com leitos de UTI neonatal no ano, proporção (0,34) | 0,987 (0,976; 0,997) |
| Intermediários | Menos de 7 consultas de pré-natal, proporção (0,12) | 1,017 (1,007; 1,027) |
| Intermediários | Idade gestacional inferior a 37 semanas, proporção (0,022) | 1,005 (0,997; 1,014) |
| Intermediários | Peso ao nascer inferior a 2.500 g, proporção (0,016) | 1,077 (1,065; 1,089) |
| Registro | Partos hospitalares, proporção (0,052) | 0,991 (0,984; 0,997) |
| Registro | Nascimentos com escolaridade materna ausente, proporção (0,032) | 1,001 (0,994; 1,009) |
| Registro | Nascimentos com raça/cor materna ausente, proporção (0,076) | 0,991 (0,983; 0,999) |

**Tabela S3.** Risco relativo residual mínimo detectável com poder de 80% (α unilateral = 0,05) apenas pela variação amostral (aproximação normal; o poder unilateral exato de Poisson com 50 óbitos esperados e risco relativo de 1,35 é de cerca de 73%), segundo os óbitos esperados na década, e número de municípios com contagem esperada igual ou superior. Acrescentar a superdispersão implicada pela heterogeneidade residual (variância E + (exp(σ²) − 1)E²) resulta em 1,53 com 50 óbitos esperados, 1,43 com 200 e 1,40 com 1.000, aproximando-se de 1,40 para contagens grandes.

| Óbitos esperados | RR mínimo detectável | Municípios com pelo menos essa expectativa | Parcela de todos os óbitos neonatais nesses municípios |
|---|---|---|---|
| 10 | 1,79 | 3.265 | 95% |
| 20 | 1,56 | 2.153 | 89% |
| 50 | 1,35 | 942 | 74% |
| 100 | 1,25 | 448 | 61% |
| 200 | 1,18 | 196 | 48% |
| 500 | 1,11 | 63 | 32% |
| 1.000 | 1,08 | 21 | 21% |
| 2.000 | 1,06 | 7 | 13% |
| 4.000 | 1,04 | 2 | 7% |

**Tabela S4.** Indicadores de registro nos municípios classificados acima do esperado. Medianas entre os 112 acima do esperado e entre os outros 5.453 municípios, e os 12 dos 112 com menos de 90% dos partos em hospital (proporção nacional 98,4%). Esses indicadores descrevem a qualidade do registro e não medem o sub-registro de nascimentos ou óbitos.

| | Acima do esperado (n = 112) | Demais municípios (n = 5.453) |
|---|---|---|
| Partos hospitalares, mediana (%) | 98,7 | 99,4 |
| Nascimentos com escolaridade materna ausente, mediana (%) | 0,93 | 0,59 |
| Nascimentos com raça/cor materna ausente, mediana (%) | 1,99 | 1,49 |

| Município (UF) | Óbitos neonatais | Risco relativo (mediana a posteriori) | Partos hospitalares (%) | Escolaridade ausente (%) | Raça/cor ausente (%) |
|---|---|---|---|---|---|
| Caracaraí (RR) | 156 | 2,25 | 83,2 | 2,2 | 0,9 |
| Barcelos (AM) | 162 | 1,91 | 59,5 | 5,9 | 1,1 |
| Amajari (RR) | 98 | 1,63 | 50,8 | 7,7 | 0,2 |
| Alto Alegre (RR) | 133 | 1,43 | 34,2 | 9,1 | 0,1 |
| Eirunepé (AM) | 162 | 1,42 | 84,9 | 1,9 | 0,2 |
| Envira (AM) | 78 | 1,33 | 88,2 | 0,0 | 0,1 |
| Aveiro (PA) | 33 | 1,32 | 61,0 | 2,2 | 1,1 |
| Uruçuca (BA) | 42 | 1,31 | 84,7 | 3,6 | 4,3 |
| Paratinga (BA) | 66 | 1,31 | 87,3 | 3,5 | 13,9 |
| Nhamundá (AM) | 55 | 1,30 | 86,5 | 19,5 | 4,5 |
| São Gabriel da Cachoeira (AM) | 199 | 1,29 | 54,8 | 0,9 | 0,3 |
| Óbidos (PA) | 150 | 1,25 | 89,7 | 0,9 | 0,6 |

A correlação de Spearman do risco relativo municipal com a proporção de partos hospitalares foi −0,06 e com a proporção de nascimentos com raça/cor ausente, −0,01.

## S4. Verificação preditiva a posteriori (modelo primário)

A partir de 1.000 amostras a posteriori, foram simuladas contagens para todos os municípios. Em 1,2% dos municípios, o valor observado ficou fora do intervalo preditivo a posteriori de 95% (valor baixo é esperado para um modelo com um efeito aleatório por município, de modo que esta verificação tem poder limitado para detectar falta de ajuste); a proporção de municípios sem óbitos foi de 1,6% observada contra 1,8% predita; a maior contagem foi de 12.293 observada (São Paulo) contra 12.294 predita (mediana); a razão variância/média das contagens foi de 1.195 observada contra 1.195 predita.

## S5. Soma dos excessos positivos sob ruído

A soma, sobre todos os municípios, da parte positiva da diferença entre óbitos observados e esperados, calculada nas expectativas da mediana a posteriori, foi 20.279. A soma, baseada no modelo, da parte positiva de $E_i(RR_i - 1)$, que contrai os municípios pequenos em direção à sua expectativa, foi 15.789 (15.039; 16.573) no modelo primário. Simulando contagens de Poisson em torno das expectativas ajustadas (entre amostras a posteriori das expectativas, ao passo que a soma observada acima usa sua mediana), sem nenhuma heterogeneidade residual, a mesma soma das partes positivas foi 10.923 (10.401; 11.462) e 4.539 (4.053; 5.075) entre os municípios com pelo menos 50 óbitos. A parte positiva do ruído é positiva, e é por isso que essa quantidade não é apresentada como carga.

## S6. Figura

**Figura S1.** Os trinta municípios com maior excesso de óbitos pela mediana a posteriori, com o risco relativo residual (mediana e intervalo de credibilidade de 95%). As linhas verticais marcam 1,00 e 1,10. A classificação como acima do esperado usa a probabilidade a posteriori P(RR>1,10)≥0,95, não o intervalo mostrado; 35 dos 112 municípios classificados têm intervalo de 95% cujo limite inferior é igual ou inferior a 1,10.

## S7. Dados suplementares

`municipal_estimates.csv`: para cada um dos 5.565 municípios, óbitos, nascimentos, óbitos esperados (mediana e ICr95%), risco relativo residual (mediana e ICr95%), probabilidades a posteriori de o risco relativo superar 1, 1,10 e 1,25, ou ficar abaixo de 0,90, óbitos em excesso (mediana e ICr95%), classificação e os três indicadores de registro.

## S8. Persistência entre períodos

**Tabela S5.** Persistência do risco relativo residual entre períodos. O modelo primário foi ajustado separadamente a cada período (óbitos, nascimentos, população e composição materna específicos do período; mesmo grafo, mesmas distribuições a priori e mesmas configurações do amostrador; o IDHM, de 2010, é o mesmo nos dois períodos).

| | 2014-2019 | 2020-2024 |
|---|---|---|
| Óbitos neonatais no conjunto de análise | 152.418 | 107.547 |
| Nascidos vivos | 17.568.151 | 12.887.307 |
| DP do log do risco relativo (ICr95%) | 0,169 (0,161; 0,178) | 0,165 (0,156; 0,176) |
| Municípios acima do esperado | 73 | 80 |
| Municípios abaixo do esperado | 47 | 9 |
| Transições divergentes; R-chapéu máximo de u | 0; 1,008 | 0; 1,010 |
| Óbitos em excesso líquido, Brasil (ICr95%) | 602 (−717; 1.857) | 1.298 (378; 2.155) |
| Óbitos em excesso líquido, Norte (ICr95%) | 465 (128; 805) | 401 (109; 673) |
| Óbitos em excesso líquido, Nordeste (ICr95%) | 1.651 (1.210; 2.130) | 975 (635; 1.329) |
| Óbitos em excesso líquido, Centro-Oeste (ICr95%) | 60 (−174; 296) | −31 (−227; 169) |
| Óbitos em excesso líquido, Sudeste (ICr95%) | −1.263 (−2.067; −506) | 312 (−232; 813) |
| Óbitos em excesso líquido, Sul (ICr95%) | −314 (−566; −59) | −359 (−551; −150) |

Acima do esperado nos dois períodos: 35; apenas em 2014-2019: 38; apenas em 2020-2024: 45. Referências para a sobreposição: 1,0 esperado se as sinalizações fossem atribuídas ao acaso entre todos os municípios (hipergeométrica, p = 3,3 × 10^-48^), e 3,0 esperados (máximo de 10 em 10.000 permutações; nenhuma alcançou 35, de modo que a estimativa de Monte Carlo é p = 0,0001) se as sinalizações de 2020-2024 fossem atribuídas ao acaso dentro de estratos de óbitos esperados agrupados (abaixo de 20, 20-49, 50-99, 100-199, 200-499, 500 ou mais). Nenhuma das duas referências preserva a aglomeração espacial. Dos 73 acima do esperado em 2014-2019, 91,8% tinham risco relativo mediano acima de 1,10 e 99% acima de 1,00 em 2020-2024; o risco relativo mediano desses municípios foi 1,37 em 2014-2019 e 1,25 em 2020-2024, como esperado quando os municípios são selecionados com base em uma estimativa. Para 31 municípios, a probabilidade conjunta de superar 1,10 nos dois períodos foi de pelo menos 0,95 (supondo ajustes independentes em dados disjuntos). Dos 112 acima do esperado na análise agrupada, 35 estavam acima nos dois períodos e 23 em nenhum. Correlação de Spearman dos riscos relativos medianos entre os períodos: 0,51 (todos os municípios), 0,61 (os 783 com pelo menos 25 óbitos em cada período), 0,58 (os 978 com pelo menos 20 óbitos esperados em cada período); correlação com a estimativa agrupada 0,88 e 0,80. Os óbitos são atribuídos pelo ano do óbito e os nascimentos pelo ano do nascimento, de modo que um óbito em janeiro de 2020 de uma criança nascida em dezembro de 2019 cai no numerador do segundo período e no denominador do primeiro; o efeito não foi quantificado. O segundo período inclui a pandemia de covid-19.

## S9. Avaliação por simulação da regra de classificação

Foram gerados quarenta conjuntos de dados com os nascimentos e as covariáveis observados, uma amostra a posteriori do modelo primário como valores verdadeiros simulados, conhecidos por construção (intercepto, coeficientes e efeitos municipais) e contagens de Poisson (39 conjuntos distintos de valores verdadeiros: duas réplicas sortearam a mesma amostra a posteriori); outros 12 conjuntos não tinham resíduo (todos os efeitos municipais iguais a zero). Cada um foi reajustado com o modelo primário (quatro cadeias com 2.000 iterações de aquecimento e 1.000 amostras retidas; sem transições divergentes nos 40 conjuntos com valores verdadeiros simulados) e comparado com os valores verdadeiros simulados. Esses valores são condicionais a um gerador calibrado no modelo ajustado e o modelo está corretamente especificado, de modo que a calibração da Tabela S7 é em parte esperada por construção; estes são resultados de uma avaliação por simulação em condições favoráveis, não uma medida de desempenho nos dados reais.

Nos 40 conjuntos com valores verdadeiros simulados, a regra sinalizou em média 114,5 municípios (intervalo 86-150) dos 1.475 com risco relativo verdadeiro acima de 1,10; 1,7% dos municípios sinalizados (80 de 4.580) tinham risco relativo verdadeiro igual ou inferior a 1,10 (por conjunto, 0%-4,8%; os erros-padrão de Monte Carlo das sensibilidades agrupadas não foram calculados). A sensibilidade global foi de 7,6% para risco relativo verdadeiro acima de 1,10 e de 20,4% acima de 1,25. Nos 12 conjuntos sem resíduo, a regra não sinalizou nenhum município acima ou abaixo do esperado. O intervalo de credibilidade de 95% cobriu o risco relativo verdadeiro em 95,0% dos municípios. O excesso líquido (média verdadeira 1.672) foi subestimado em 247 em média (14,8%), com o intervalo cobrindo o valor verdadeiro em 95% dos conjuntos; o excesso entre os municípios sinalizados (média verdadeira 6.065, média estimada 5.963) foi subestimado em 102 (1,7%), com o intervalo cobrindo o valor verdadeiro em 80% dos conjuntos, abaixo dos 95% nominais.

**Tabela S6.** Sensibilidade da regra e cobertura do intervalo de credibilidade de 95%, segundo os óbitos esperados verdadeiros do município (agrupado entre os conjuntos com valores verdadeiros simulados).

| Óbitos esperados | Municípios com RR verdadeiro > 1,10 | Sinalizados | Sensibilidade | RR verdadeiro > 1,25 | Sinalizados | Sensibilidade | Cobertura do intervalo de 95% |
|---|---|---|---|---|---|---|---|
| <20 | 36.576 | 801 | 2,2% | 10.274 | 670 | 6,5% | 95,0% |
| 20-49 | 13.363 | 1.437 | 10,8% | 4.374 | 1.217 | 27,8% | 94,8% |
| 50-199 | 7.323 | 1.385 | 18,9% | 1.955 | 1.041 | 53,2% | 95,1% |
| ≥200 | 1.722 | 877 | 50,9% | 609 | 588 | 96,6% | 95,4% |

**Tabela S7.** Calibração grosseira da probabilidade de excedência P(RR>1,10): proporção de municípios com risco relativo verdadeiro acima de 1,10, por faixa de probabilidade estimada (agrupado entre os conjuntos com valores verdadeiros simulados; a probabilidade média predita dentro das faixas não foi calculada).

| P(RR>1,10) estimada | Municípios-conjuntos | Proporção realmente acima de 1,10 |
|---|---|---|
| 0,00 a 0,05 | 45.420 | 1,7% |
| 0,05 a 0,50 | 139.058 | 22,2% |
| 0,50 a 0,80 | 25.885 | 62,4% |
| 0,80 a 0,95 | 7.657 | 87,2% |
| 0,95 a 1,00 (sinalizados) | 4.580 | 98,3% |

## S10. Comparação com outras abordagens

**Tabela S8.** Municípios sinalizados acima do esperado (P(RR>1,10)≥0,95) pelo modelo primário (regressão espacial restrita com estrutura e composição materna), por dois modelos apenas com estrutura (não espacial com efeito permutável; espacial com a mesma distribuição a priori BYM2 restrita) e pela taxa bruta de mortalidade neonatal.

| Comparação | Número |
|---|---|
| Sinalizados pelo modelo primário | 112 |
| Sinalizados pelo modelo não espacial apenas com estrutura | 74 |
| Sinalizados pelos dois | 58 |
| Somente pelo modelo primário | 54 (45 no Nordeste) |
| Somente pelo modelo não espacial apenas com estrutura | 16 |
| Correlação de Spearman dos riscos relativos municipais, primário versus não espacial apenas com estrutura | 0,76 |
| Sinalizados pelo modelo espacial apenas com estrutura | 153 |
| Sinalizados pelo modelo primário e pelo espacial apenas com estrutura | 93 |
| Somente pelo modelo primário (acrescentados pela composição materna) | 19 |
| Somente pelo modelo espacial apenas com estrutura (removidos pela composição materna) | 60 |
| Sinalizados pelo modelo primário com pelo menos 50 óbitos | 67 |
| ... dos quais entre as 112 maiores taxas brutas (municípios com pelo menos 50 óbitos) | 46 |
| ... dos quais entre as 67 maiores taxas brutas (comparação equivalente) | 38 |
| Número esperado de municípios sinalizados com risco relativo verdadeiro igual ou inferior a 1,10 (baseado no modelo, soma de 1 − P(RR>1,10)) | 1,8 (1,6%) |
