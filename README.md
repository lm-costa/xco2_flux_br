# OCO-2 v11 MIP: XCO2 × fluxos gradeados (Brasil)

Pipeline em R que compara o fluxo de CO2 estimado a partir do XCO2 do OCO-2 MIP com o `net_flux` gradeado do próprio MIP.



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
delta → fluxo pontual por sondagem                      (inventário de coluna, τ = 1 dia)
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
4. **Delta e fluxo pontual.** `delta = xco2 − xco2_ens_mean`, com o simulado do experimento **IS**. O delta é convertido em fluxo por sondagem como inventário de coluna, sem vento, com tempo de residência τ = 1 dia (default). A magnitude depende de τ.
5. **Agregação (`agg_xco2_flux`).** Cada sondagem vai para a célula do fluxo mais próxima (centro da célula) e para o mês do seu `datetime`. O método padrão é `dois_estagios`: média por célula × dia e depois média dos dias do mês, de modo que cada dia pesa igual. `n_days` é o número de dias independentes; `<var>_mean`, `<var>_sd` e `<var>_sem` são as estatísticas por célula × mês.
6. **Comparação (`compare_flux`, `map_flux_annual`).** Estimado × `net_flux` no Brasil e por bioma: correlação, RMSE, viés, balanço mensal/anual/total e regressão célula a célula. Os mapas anuais ficam em `map_flux_annual()`.
7. **Método beta (`beta_xco2_flux`, `compare_beta`, `balance_beta`).** Tendência do XCO2 por célula e ano convertida em fluxo, reamostrada para 1° e comparada com o estimado e o `net_flux`. Só as células presentes nas duas bases são mantidas (co-locadas).

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
| `comparison.rds` | saída de `agg_xco2_flux()`: estimado + `net_flux` | 1 linha por célula × mês |

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

| Script | Função principal |
|---|---|
| `agg_xco2_flux.R` | `agg_xco2_flux()`: agrega por célula × mês e junta com o fluxo |
| `compare_flux.R` | `compare_flux()`, `make_flux_plots()`, `save_flux_plots()`: métricas, balanços (mensal/anual/total), séries e scatter |
| `map_flux_annual.R` | `map_flux_annual()`, `save_flux_maps()`: mapas anuais de estimado, referência, viés, r, RMSE e inclinação |
| `beta_xco2_flux.R` | `beta_xco2_flux()`, `pair_beta_flux()`: método beta (tendência do XCO2) |
| `compare_beta.R` | `compare_beta()`, `balance_beta()`, `save_beta_plots()`, `save_beta_maps()`, `save_beta_balance()`: beta × estimado × `net_flux` |

Ordem de uso típica:

```r
comp <- agg_xco2_flux(xco2, flux, method = "dois_estagios", out_file = "teste_comp.rds")
res  <- compare_flux(comp, biomes = biomas)
maps <- map_flux_annual(comp, biomes = biomas)
rb   <- compare_beta(dfall, comp, biomes = biomas, unit = "physical", maps = TRUE)
bal  <- balance_beta(rb)
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

- Tempo de residência τ (1 dia): é uma hipótese, não um erro aleatório. A magnitude do fluxo escala com ela.
- Erro sistemático do XCO2 simulado do experimento IS (viés do modelo de transporte), que entra no delta como se fosse sinal.
- Unidade do `net_flux` (`ref_factor`).
- O método beta não tem propagação de incerteza; só há o p-valor da regressão, que fica inflado por autocorrelação.

---

## 7. Ressalvas abertas

- Magnitude do fluxo estimado depende de τ (1 dia).
- Poucos anos por célula (5–9): regressões fracas e p-valores inflados por autocorrelação.
- Beta ignora o transporte.
