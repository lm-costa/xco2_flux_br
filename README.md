# OCO-2 v11 MIP: XCO2 × fluxos gradeados (Brasil)

Pipeline em R que estima o fluxo de CO2 do Brasil (2015–2023) a partir do XCO2 do OCO-2, **sem cálculos de transporte atmosférico**, e o compara com o `net_flux` gradeado do OCO-2 MIP v11.

**Ideia.** O XCO2 simulado pelo experimento do MIP que assimila só dados *in situ* (`IS`) serve de fundo. A diferença `delta = XCO2 observado (OCO-2) − XCO2 simulado (IS)` carrega o sinal de trocas próximas à superfície que a rede *in situ* não restringe, e é convertida em fluxo por um balanço de massa na coluna (fração `f_a` da coluna e tempo de residência `τ`). O fluxo estimado é agregado nas células de 1° do MIP e comparado, por bioma, ano e mês, com o fluxo do experimento `LNLGIS` (que também assimila OCO-2), com desmatamento (PRODES), área queimada (AQ1KM) e clima (NASA POWER).

**Status.** Trabalho exploratório. Os resultados mostram concordância razoável no balanço anual (r = 0,61 para o Brasil) e fraca na escala mensal, com viés positivo; a abordagem não substitui as inversões. Veja as limitações no manuscrito e na seção 8.

**Índice:** 1. Visão geral · 2. Arquivos `.rds` · 3. Experimentos · 4. Submissions · 5. Scripts · 6. Incerteza · 7. Ressalvas · 8. Análises adicionais · 9. Estrutura e como rodar · 10. Fontes de dados.



---

## 1. Visão geral do processamento

```
NetCDF do MIP (por submission × experimento)
   │  mip_extractor()
   ▼
.rds de XCO2 por sondagem (Brasil)                     [1 linha = 1 sondagem]
   │  build_ensemble()
   ▼
ensemble__<EXP>.rds   (obs + média/dp do ensemble)     [1 linha = 1 sondagem]
   │  delta = xco2 − xco2_ens_mean   (experimento IS)
   ▼
delta → fluxo pontual por sondagem                      (balanço de massa: f_a = 0,2; τ = 2,5 dias)
   │  agg_xco2_flux()   (célula do net_flux × mês)
   ▼
teste_comp.rds   (estimado + net_flux, 1 linha = célula × mês)
   │  compare_flux() / map_flux_annual()
   ▼
métricas, balanços (Tg C), séries, mapas

NetCDF de fluxos gradeados do MIP
   │  flux_extractor()
   ▼
EnsMean__<EXP>.rds, EnsStd__<EXP>.rds                   [1 linha = célula × mês]
```

### Etapas

1. **Extração do XCO2 (`mip_extractor`).** Lê os NetCDF do MIP de cada submission × experimento e recorta o Brasil.
2. **Ensemble (`build_ensemble`).** Combina as submissions de um experimento e gera, por sondagem, a média (`xco2_ens_mean`) e o desvio padrão (`xco2_ens_sd`) do XCO2 simulado, junto com o XCO2 observado (`xco2`).
3. **Extração dos fluxos (`flux_extractor`).** Lê os fluxos gradeados (1°, `net_flux`) e salva a média do ensemble (`EnsMean`) e o desvio entre membros (`EnsStd`) por célula × mês.
4. **Delta e fluxo pontual (`fco2_from_delta`).** `delta = xco2 − xco2_ens_mean`, com o simulado do experimento **IS**. O delta vira fluxo por sondagem por um balanço de massa na coluna, sem vento:
   `n_dry = (Ps/g − TCWV) / M_ar` (ar seco por m², com a pressão de superfície e o vapor d'água da coluna de cada sondagem), `dM = f_a · delta · 1e-6 · n_dry · M_C` e `fco2 = dM / τ` (g C m⁻² d⁻¹). No `main.R` usa-se `f_a = 0,2` (`flux_fraction`) e `τ = 2,5 dias` (`tau_days`). O fluxo escala com `f_a/τ`: os dois parâmetros não são separáveis pelos dados e controlam a amplitude, não a correlação com o MIP.
5. **Agregação (`agg_xco2_flux`).** Cada sondagem vai para a célula do fluxo mais próxima (centro da célula) e para o mês do seu `datetime`. O método padrão é `dois_estagios`: média por célula × dia e depois média dos dias do mês, de modo que cada dia pesa igual. `n_days` é o número de dias independentes; `<var>_mean`, `<var>_sd` e `<var>_sem` são as estatísticas por célula × mês.
6. **Comparação (`compare_flux`, `map_flux_annual`).** Estimado × `net_flux` no Brasil e por bioma: correlação, RMSE, viés, balanço mensal/anual/total e regressão célula a célula. Os mapas anuais ficam em `map_flux_annual()`.
7. **Método beta (legado: `beta_xco2_flux`, `compare_beta`, `balance_beta`).** Tendência do XCO2 por célula e ano convertida em fluxo, reamostrada para 1° e comparada com o estimado e o `net_flux`. É a abordagem da primeira versão do trabalho, mantida só como referência (o código fica comentado no `main.R`). A análise com SIF foi retirada do projeto.

### Convenções e unidades

- Grade de fluxo: 1°, centros em x.5. A área da célula é fixa em 110 × 110 km.
- Balanços: Tg C = fluxo × área × dias, somado nas células pareadas.
- Delta e beta são calculados com produtos de XCO2 diferentes (MIP v11 vs. Lite FP v11.1 no repositório do beta). Isso é uma ressalva da comparação.

---

## 2. Nomes dos arquivos `.rds`

Padrões observados nos scripts:

| Arquivo | Conteúdo | Granularidade |
|---|---|---|
| `<Submission>__<EXP>.rds` (ex.: `data/xco2_pre_processed/Ames__IS.rds`) | XCO2 observado + XCO2 simulado nesse `<EXP>` | 1 linha por sondagem |
| `ensemble__<EXP>.rds` (ex.: `data/xco2_pre_processed/ensemble/ensemble__IS.rds`) | XCO2 observado + média e desvio do ensemble simulado do experimento `<EXP>` | 1 linha por sondagem |
| `<Submission>__<EXP>.rds` (ex.: `data/Fluxes_MIP_completo/Ames__IS.rds`) | `net_flux`  do `<Submission>` para o experimento `<EXP>` | 1 linha por célula × mês |
| `EnsMean__<EXP>.rds` (ex.: `data/Fluxes_MIP_completo/EnsMean__LNLGIS.rds`) | `net_flux` médio do ensemble para o experimento `<EXP>` | 1 linha por célula × mês |
| `EnsStd__<EXP>.rds` (ex.: `data/Fluxes_MIP_completo/EnsStd__LNLGIS.rds`) | desvio padrão entre membros do `net_flux` do experimento `<EXP>` (candidato a `ref_sigma`) | 1 linha por célula × mês |
| `comparison.rds` | saída de `agg_xco2_flux()` da versão anterior (τ = 30 d, sem `f_a`) | 1 linha por célula × mês |
| `comparison_new.rds` | saída de `agg_xco2_flux()` com `f_a = 0,2` e `τ = 2,5 d` (base das análises atuais) | 1 linha por célula × mês |
| `data/fire/aq1km_cell_month.rds` | área queimada do AQ1KM agregada em células de 1° (listas `monthly`, `annual`, `land`, `total_biome`) | célula × mês / ano |
| `data/prodes/prodes_cell_year.rds` | desmatamento do PRODES agregado em células de 1° (listas `cell`, `total`) | célula × ano |
| `data/nasa_power/power_cell_month.rds` | NASA POWER mensal por célula (prec, qg, rh, temp, ws) | célula × mês |
| `output/*.rds`, `output/*.csv` | resultados das análises (métricas, balanços, comparações). `metrics_temporal.rds` e `metrics_spatial.rds` não são versionados (> 100 MB) e se regeneram no `main.R` | variável |

Regra geral: `<TIPO>__<EXPERIMENTO>.rds`, com **dois sublinhados** entre o tipo (`ensemble`, `EnsMean`, `EnsStd`) e o experimento.

**TODO:** padrão dos arquivos individuais por submission × experimento (antes do ensemble), por exemplo `<SUBMISSION>__<EXP>.rds`. Preencha com o nome real.

---

## 3. Experimentos (`exp`)

Definição do protocolo OCO-2 MIP (confirme contra a documentação do MIP que você está usando):

| Código | Dados assimilados |
|---|---|
| `IS` | apenas medidas in situ (superfície/aeronave); **não assimila OCO-2**. É a referência de XCO2 simulado usada no delta (`obs − IS`). |
| `LNLG` | OCO-2 sobre terra: nadir + glint |
| `LNLGIS` | `LNLG` + in situ. É o experimento do `net_flux` usado como referência neste projeto |
| `LNLGOGIS` | `LNLG` + glint sobre o oceano + in situ |

Resumo da lógica: `IS` não vê o OCO-2, então `obs − IS` mede o que o satélite "vê" além do que os dados in situ já explicam. É esse sinal que vira fluxo no método próprio.



---

## 4. Submissions

Cada submission é o resultado de um grupo de inversão participante do OCO-2 v11 MIP. Cada uma roda os experimentos da seção 3 (resultados de 2015 a 2024) e é identificada pelo nome do modelo. A lista abaixo segue a tabela oficial do MIP (https://gml.noaa.gov/ccgg/OCO2_v11mip/index.php), atualizada em 2026-08-17. O ensemble combina as submissions de cada experimento.

| Submission | Instituição | Modelo de transporte | Meteorologia | Método de inversão |
|---|---|---|---|---|
| Ames | NASA Ames Research Center | GEOS-Chem | MERRA-2 | 4D-Var |
| CAMS | LSCE (França) | LMDz-Dispersion | ERA5 (via LMDz) | Variacional |
| CMS-Flux | NASA JPL | GEOS-Chem | MERRA-2 | 4D-Var |
| CMS-MFlux | NASA JPL | GEOS-Chem | MERRA-2 | Analítico / multirresolução |
| CT | Univ. do Colorado e NOAA GML | TM5-zoom | ERA5 | EnKF |
| CTE-FMI | Finnish Meteorological Institute | TM5-MP | ERA5 | EnKF |
| GCAS | Chinese Academy of Sciences | GEOS-Chem | MERRA-2 | NLS-4DVar |
| GONGGA | Nanjing University | MOZART-4 | GEOS5 | EnKF |
| JHU | Johns Hopkins University | GEOS-Chem | MERRA-2 | Geoestatístico |
| MAGI | NASA JPL | GEOS-Chem | MERRA-2 | Analítico |
| NISMON-CO2 | NIES | NICAM-TM | JRA-3Q | 4D-Var |
| NTFVAR | NIES Satellite Observation Center | NIES-TM e FLEXPART | ERA5 e JRA-55 | Variacional |
| PCTM | Colorado State University | PCTM | MERRA-2 | 4D-Var |
| TM5-4DVAR | Univ. de Maryland e NASA GMAO | TM5 | ERA5 | 4D-Var |
| WOMBAT-GC | Univ. of Western Australia e Colorado State Univ. | GEOS-Chem | MERRA-2 | Bayesiano hierárquico |
| WOMBAT-TM5 | Univ. of Western Australia e Univ. do Colorado | TM5 | ERA5 | Bayesiano hierárquico |

Observações do site que afetam o uso dos dados:

- **GONGGA:** os cosamples não passaram nas verificações, com viés negativo em SPO e ajuste pior da taxa de crescimento global.
- **GCAS:** o submission v2 ainda não foi checado quanto às correções de problemas anteriores.
- **MAGI:** usa o OCO-2 v1r2.
- **NISMON-CO2:** o LNLG usado é a versão 2.
- Os demais passam nas quatro verificações de cosample (ObsPack, TCCON, OCO-2 e OCO-3).
- Cada submission tem número de versão e revisão no site; use a versão que está nos seus arquivos.


---

## 5. Scripts

Os scripts ficam em `r/` e são carregados pelo `main.R` (`purrr::map(list.files('r'), source)`).

| Script | Funções principais |
|---|---|
| `setup_pack.R` | `setup_packages()`: instala/carrega os pacotes |
| `mip_extractor.R`, `mip_ensemble.R`, `flux_extractor.R` | `mip_extractor()`, `mip_ensemble()` (ensemble do XCO2 por sondagem), `flux_extractor()` (fluxos gradeados) |
| `fco2_from_delta.R` | `fco2_from_delta()`: delta → fluxo (balanço de massa com `f_a` e `τ`) e propagação de incerteza |
| `aggregation_xco2.R` | `agg_xco2_flux()`: agrega por célula × mês e junta com o fluxo |
| `compare_flux.R` | `compare_flux()`, `make_flux_plots()`, `save_flux_plots()`: métricas, balanços (mensal/anual/total), séries e dispersão |
| `map_flux_annual.R`, `map_theme.R` | `map_flux_annual()`, `save_flux_maps()`: mapas anuais de estimado, referência, viés, r, RMSE, inclinação e incerteza |
| `plot_lat_profile.R` | `plot_lat_profile()`: perfil latitudinal do fluxo, uma linha por ano |
| `compare_climate.R` | `compare_climate()`, `save_climate_plots()`, `plot_climate_maps()`: fluxos × anomalias do NASA POWER; mapas das médias climáticas |
| `disturbance.R` | `aq1km_download()`, `aq1km_aggregate()`, `prodes_aggregate()`, `map_disturbance_annual()`, `save_disturbance_maps()`, `compare_disturbance()`: fogo e desmatamento |
| `compare_models.R` | `compare_models()`, `plot_models()`, `plot_biome_scatter()`, `compare_models_monthly()`, `plot_models_flux()`, `plot_seasonal_cycle()`: comparação modelo a modelo do MIP |
| `compare_increment.R` | `mip_increment()`, `compare_increment()`, `plot_increment()`, `increment_disturbance()`, `plot_increment_disturbance()`: comparação com o incremento do OCO-2 (LNLGIS − IS) |
| `beta_fco2.R`, `compare_beta.R` | método beta (legado; ver seção 1) |

Ordem de uso típica:

```r
df_n  <- fco2_from_delta(df_n, delta = "delta", sigma_delta = "sigma", method = "massa",
                         tau_days = 2.5, flux_fraction = 0.2)
comp  <- agg_xco2_flux(df_n, df_flux, vars = "fco2", method = "dois_estagios", sigma = "sigma_f",
                       out_file = "data/comparison_new.rds")
res   <- compare_flux(comp, biomes = biomas, years = 2015:2023, ref_sigma = "net_flux_sd")
maps  <- map_flux_annual(comp, biomes = biomas, years = 2015:2023)
```

---

## 6. Incerteza do fluxo estimado a partir do delta

A incerteza é propagada em três níveis: sondagem → célula × mês → balanço (Tg C).

**1) Por sondagem (`sigma` → `sigma_f`).** O delta e a incerteza total do delta são calculados por sondagem:

```r
df_n <- df |>
  mutate(
    delta = xco2 - xco2_ens_mean,
    sigma = sqrt((xco2_ens_se^2) + (uncertanty^2) + (model_error^2))
  )
```

`sigma` é 1 desvio padrão do delta (ppm), somando em quadratura três termos, tratados como independentes entre si:

| Termo | Significado |
|---|---|
| `xco2_ens_se` | erro padrão da média do ensemble de XCO2 simulado (dispersão entre submissions) |
| `uncertanty` | incerteza da medida de XCO2 do OCO-2 (coluna `uncertanty` da base) |
| `model_error` | erro de modelo (transporte e representação) |

`fco2_from_delta()` converte `delta` em `fco2` (g C m⁻² d⁻¹) e leva `sigma` para a mesma unidade, gerando `sigma_f` (1 desvio padrão do fluxo da sondagem). Se a conversão é linear no delta, `sigma_f` é `sigma` multiplicado pelo mesmo fator que leva `delta` a `fco2`.

**2) Célula × mês (`agg_xco2_flux(..., sigma = "sigma_f")`).** Depende do método de agregação:

- `dois_estagios` (padrão): dentro de um dia os erros são tratados como **totalmente correlacionados**, então σ_dia = média dos `sigma_f` das sondagens do dia. Entre dias os erros são **independentes**:
  `σ_mês = sqrt(Σ σ_dia²) / n_days` → coluna `sigma_f_cell`.
- `direto`: média simples das sondagens, com três colunas conforme a correlação ρ entre sondagens: `sigma_f_cell_indep` (ρ = 0, `σ/√n`), `sigma_f_cell_corr` (ρ = 1) e, se `rho` for informado, `sigma_f_cell_rho` (`σ·√(1/n + (1 − 1/n)·ρ)`).
- Checagem empírica: `fco2_sem = fco2_sd / √n_days` (dispersão real entre dias). Se `fco2_sem` for bem maior que `sigma_f_cell`, a incerteza propagada está subestimada.

**3) Balanço (`compare_flux`).** O balanço de cada célula × mês é `fluxo × área × dias do mês`. A incerteza dele sai em duas versões, que são limites:

| Versão | Hipótese entre células | Cálculo |
|---|---|---|
| `*_sigma_indep_TgC` | células independentes (limite inferior) | `sqrt(Σ (σ·w)²)` |
| `*_sigma_corr_TgC` | células totalmente correlacionadas (limite superior) | `Σ σ·w` |

Do mensal para o anual e o total: a versão `indep` soma em quadratura (`sqrt(Σ σ²)`) e a `corr` soma direto. No gráfico do balanço (`sigma_bars`), a barra grossa é a `indep` e a haste fina é a `corr`. Com erros correlacionados a haste costuma ser muito maior que o próprio balanço; `sigma_bars = "indep"` deixa o gráfico legível, mas mostra só o limite inferior.

**Referência (`net_flux`).** A incerteza do MIP entra por `ref_sigma` (coluna de desvio padrão por célula × mês, em geral a do arquivo `EnsStd__<EXP>.rds`; `ref_factor` é aplicado como no `net_flux`). Sem `ref_sigma`, a incerteza da referência não entra. A incerteza da **diferença** (estimado − referência) é `sqrt(σ_est² + σ_ref²)`, nas duas versões.

**O que a incerteza não cobre**

- Fração da coluna `f_a` (0,2) e tempo de residência `τ` (2,5 dias): são hipóteses, não erros aleatórios, e não entram na incerteza propagada (a menos que se informe `sigma_flux_fraction` e `sigma_tau_days`). A magnitude do fluxo escala com `f_a/τ`.
- Erro sistemático do XCO2 simulado do experimento IS (viés do modelo de transporte), que entra no delta como se fosse sinal.
- Unidade do `net_flux` (`ref_factor`).
- O método beta (legado) não tem propagação de incerteza; só há o p-valor da regressão, que fica inflado por autocorrelação.

---

## 7. Ressalvas abertas

- A magnitude do fluxo estimado depende de `f_a` e `τ` (constantes no espaço e no tempo, sem calibração).
- O transporte atmosférico não é modelado: o ΔXCO2 de uma célula reflete também os fluxos a montante. A concordância melhora em médias anuais e regionais.
- Sinal e ruído: com os parâmetros usados, um fluxo de 0,3 g C m⁻² d⁻¹ equivale a cerca de 0,9 ppm de ΔXCO2, parecido com a incerteza de uma sondagem (≈ 1 ppm). Um deslocamento do fundo (≈ +0,27 ppm) é indistinguível de uma fonte líquida persistente.
- Amostragem: só há células × meses com sondagens válidas (cerca de 33% da área e do tempo do Brasil, 22% na Amazônia); os balanços anuais se referem às células pareadas, não aos biomas inteiros.
- O fundo (IS) e o fluxo de referência (LNLGIS) vêm do mesmo MIP, então a comparação não é independente. O fluxo do `LNLGIS` contém também a parte já restrita pelos dados *in situ*; ver o teste com `LNLGIS − IS` na seção 8.
- Não há validação com observações independentes (perfis aéreos, torres de fluxo) neste repositório.
- Poucos anos por célula (5–9): correlações anuais com fogo, desmatamento e clima são fracas e com IC largos.

---

## 8. Análises adicionais

Todas rodam a partir de `data/comparison_new.rds` e estão nas seções do `main.R`:

| Seção do `main.R` | O que faz | Saídas |
|---|---|---|
| 5 | `compare_flux()` e `map_flux_annual()`: estimado × `net_flux` do LNLGIS por bioma, ano e mês; mapas de viés, r e RMSE | `output/metrics_temporal.rds`, `output/metrics_spatial.rds`, `figs/` |
| 6 e 6b | clima (NASA POWER): correlação de anomalias anuais com o fluxo estimado e o do MIP; mapas das médias climáticas | `output/compare_climate_*.rds`, `figs_climate/` |
| 7 | fogo (AQ1KM) e desmatamento (PRODES) agregados em 1°, mapas anuais e associação (Spearman, IC por bootstrap) com os fluxos | `data/fire/`, `data/prodes/`, `figs_disturbance/`, `output/compare_*.rds` |
| 8 e 8b | comparação modelo a modelo do MIP (balanço anual, médias mensais e ciclo sazonal) | `output/compare_models*.rds`, `figs_manuscript/` |
| 9 e 9b | **exploratório:** comparação com o incremento do OCO-2, `LNLGIS − IS` (grupo a grupo, 16 grupos), e sua associação com fogo e desmatamento | `output/compare_increment.rds`, `output/increment_disturbance.rds`, `figs_increment/` |

As figuras do artigo e do material suplementar estão em `figs_manuscript/`, `figs_climate/`, `figs_disturbance/` e `figs_increment/`. O teste com o incremento (seções 9 e 9b) ainda é uma hipótese: depende de como o incremento é definido e não está no texto principal do artigo.

---

## 9. Estrutura do repositório e como rodar

```
main.R            roteiro completo (seções 0 a 9b)
r/                funções (ver seção 5)
data/             dados processados (rds) usados nas análises
data-raw/         dados brutos (não versionados; ver seção 10)
output/           resultados das análises
figs*/            figuras (figs, figs_climate, figs_disturbance, figs_increment, figs_manuscript)
```

1. Abra `xco2_flux_br.Rproj` e rode `main.R` por seções (as etapas pesadas — extração do MIP, download do AQ1KM e do PRODES, download do NASA POWER — estão comentadas e já têm os resultados salvos em `data/`).
2. Os pacotes são instalados e carregados por `setup_packages()` (por exemplo: dplyr, ggplot2, sf, geobr, ncdf4, readxl, vroom, nasapower, ggspatial, patchwork).
3. Para repetir tudo desde o início, baixe os dados brutos para `data-raw/` e descomente as etapas de extração (`mip_extractor()`, `flux_extractor()`, `aq1km_download()`, `prodes_aggregate()`, `compare_climate(download = TRUE)`).
4. Os arquivos `output/metrics_temporal.rds` e `output/metrics_spatial.rds` são gerados pela seção 5 e estão no `.gitignore` (acima do limite de 100 MB por arquivo do GitHub).

---

## 10. Fontes de dados

| Dado | Fonte |
|---|---|
| XCO2 do OCO-2 (L2 Lite FP, v11.2r) | NASA GES DISC: https://disc.gsfc.nasa.gov/datasets/OCO2_L2_Lite_FP_11.2r/summary |
| OCO-2 MIP v11 (XCO2 co-amostrado e fluxos gradeados) | NOAA GML: https://gml.noaa.gov/ccgg/OCO2_v11mip/ (DOI 10.15138/V5FS-0616) |
| Desmatamento (PRODES) | INPE TerraBrasilis: https://terrabrasilis.dpi.inpe.br/en/download-files/ |
| Área queimada (AQ1KM, Coleção 2) | INPE Queimadas: https://dataserver-coids.inpe.br/queimadas/queimadas/area_queimada/colecao2/tif/ |
| Clima (precipitação, umidade, temperatura, radiação, vento) | NASA POWER: https://power.larc.nasa.gov/ (pacote `nasapower`) |
| Limites dos biomas (IBGE, 2019) | pacote `geobr`: https://github.com/ipeaGIT/geobr |
