############################################
#   0. loading functions and libraries     #
############################################

biomas   <- geobr::read_biomes(year = 2019)


purrr::map(list.files('r',full.names = T),source)
setup_packages()

biomas <- biomas |>
  filter(
    name_biome!='Sistema Costeiro'
  )

south_file <- list.files('data-raw/South_America/',pattern = 'shp',full.names = T)
south_america <- sf::read_sf(south_file[1])

############################################
#       1. Pre - processesing              #
############################################

# YOU DON'T NEED TO RUN THIS SECTION EVERYTIME
# ALSO, ALL PRE PROCESSED DATA ARE AVAILABLE AT DATA FOLDER

## XCO2 submodel and expirement extraction considering a geometry

mip_extractor(
  nc_file = 'data-raw/OCO2.nc', # raw data from <https://gml.noaa.gov/ccgg/OCO2_v11mip/release_downloads.php>
  out_dir = 'data/xco2_pre_processed', #output folder
  geometry = geobr::read_country(year = 2025), # geometry used

)

## ensemble of XCO2 data by experiments
mip_ensemble(
  in_dir = 'data/xco2_pre_processed',
  experiments = 'IS'
)

## Flux submodel and expirement extraction considering a geometry

# df_flux <- flux_extractor(
#   in_dir      = 'data-raw/OCO2_v11MIP_gridded_fluxes_all_20260729.v2r3',
#   out_dir     = "data/Fluxes_MIP",
#   geometry    = geobr::read_country(year = 2025),
#   submissions = "EnsMean",
#   experiments = "LNLGIS"
# )
#
# df_flux_std <- flux_extractor('data-raw/OCO2_v11MIP_gridded_fluxes_all_20260729.v2r3',
#                               "data/fluxes_MIP_std",
#                               geometry = geobr::read_country(year = 2025),
#                               submissions = "EnsStd",
#                               experiments = "LNLGIS")

df_flux <- list.files('data/Fluxes_MIP',full.names = T)
df_flux_std <- list.files('data/Fluxes_MIP_std',full.names = T)

df_flux <- readr::read_rds(df_flux) # or you can insert the path directly. e.g: 'data/Fluxes_MIP_completo/<Submission>__<Exp>.rds

df_std <- readr::read_rds(df_flux_std) |> # or you can insert the path directly. e.g: 'data/Fluxes_MIP_completo/<Submission>__<Exp>.rds
  dplyr::transmute(lon = round(lon, 4), lat = round(lat, 4),
                   year, month, net_flux_sd = net_flux)

#####

df <- readr::read_rds(
  list.files('data/xco2_pre_processed/ensemble',
             full.names = T)
)

## general view
df |>
  sample_n(1000) |>
  ggplot(
    aes(x=lon,y=lat,col=xco2)
  )+
  geom_point()

############################################
#       2. Delta and total uncertanty      #
############################################
df_n <- df |>
  #rowwise() |>
  mutate(
    delta = xco2 - xco2_ens_mean,
    sigma = sqrt((xco2_ens_se^2)+(uncertanty^2)+(model_error^2))
  )

df_n |>
  dplyr::mutate(
    date = lubridate::date(datetime),
    year = lubridate::year(date),
    month = lubridate::month(date)
  ) |>
  dplyr::filter(year %in% 2015:2023) |>
  dplyr::mutate(
    date = lubridate::make_date(year,month,'15')
  ) |>
  dplyr::group_by(date) |>
  dplyr::summarise(xco2_mean=mean(delta)) |>
  ggplot2::ggplot(ggplot2::aes(x=date,y=xco2_mean))+
  ggplot2::geom_point(shape=21,color="black",fill="gray") +
  ggplot2::geom_line(color="red")+
  ggplot2::geom_smooth(method = "lm") +
  #ggplot2::ylim(390,420)+
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        ggplot2::aes(label = paste(..eq.label..,
                                                   ..rr.label..,
                                                   #..p.value.label..,
                                                   sep = "*`,`~")),
                        label.y = 0.01,
                        parse = TRUE
  )+
  #ggplot2::facet_wrap(~year,scales ='free_x')+
  ggplot2::theme_bw()+
  ggplot2::labs(x='',y=expression(Delta~'Xco'[2]~' (ppm)'),fill='' )

############################################
#       3. Fluxes estimative               #
############################################

df_n <- fco2_from_delta(
  df_n,
  delta='delta',
  sigma_delta = 'sigma',
  method = 'massa',
  tau_days = 2.5,
  flux_fraction = .2

)


############################################
# 4. Aggregation of OCO-2 to MIP 1deg      #
############################################

df_comp <- agg_xco2_flux(
  df_n,
  df_flux,
  vars = 'fco2',
  method = 'dois_estagios',
  sigma = 'sigma_f',
  min_days = 1,
  join='left',
  out_file = 'data/comparison_new.rds',

)

df_comp <- dplyr::left_join(
  dplyr::mutate(df_comp |>
                  filter(
                    fco2_mean > -15

                    ), lon = round(lon, 4), lat = round(lat, 4)),
  df_std, by = c("lon", "lat", "year", "month")
  )



############################################
# 5. Compare fluxes estimated with MIP     #
############################################

### temporal

res <- compare_flux(
  df_comp |> filter(
    fco2_mean > -15
  ),
  flux_all   = df_flux,
  biomes     = biomas,
  min_days   = 1,            # exige >= 3 dias com sondagem na célula x mês
  years      = 2015:2023,    # período em que existe net_flux
  ref_factor = 1,            # ajuste se o net_flux não estiver em g C m-2 d-1,
  ref_sigma = "net_flux_sd",
  sigma_bars = "indep"
)




res$metrics          # Brasil + biomas: bias, rmse, r, slope...
res$balance_total    # balanço acumulado (Tg C): estimado x net_flux
res$balance_year     # por ano
res$flux_monthly |> view()     # fluxo médio mensal (g C m-2 d-1) por grupo
res$plots$scatter; res$plots$balance; res$plots$monthly;res$plots$balance_year   # por bioma, em inglês

res |> readr::write_rds('output/metrics_temporal.rds')

save_flux_plots(res, dir = "figs/")                       # salva PNGs

#### Spatial


out <- map_flux_annual(df_comp, biomes = biomas,
                       south_america = south_america,   # seu objeto sf
                       map_theme = map_theme,         # seu tema
                       est_sigma = 'sigma_f_cell',
                       ref_sigma = 'net_flux_sd',
                       years = 2015:2023,
                       min_months = 3,   # meses pareados mínimos na célula x ano
                       min_years  = 3,)   # anos mínimos para a regressão da célula

out$summary         # Brasil: medianas de r, slope, rmse, bias; fração de células com p < 0.05
out$summary_biome   # o mesmo por bioma
out$regression      # métricas célula a célula (lon, lat, slope, r, p, rmse, bias...)
out$maps$flux_estimated; out$maps$flux_reference;
out$maps$sigma_estimated;out$maps$sigma_reference;
out$maps$bias; out$maps$r; out$maps$slope;out$maps$rmse

out |> readr::write_rds('output/metrics_spatial.rds')

save_flux_maps(out, "figs/")

#### Latitudinal distribution (one line per year)

lat <- plot_lat_profile(out$annual,
                        bin = 2,          # largura da faixa de latitude (graus)
                        min_cells = 5,    # mínimo de células por faixa x ano
                        min_months = 3,   # mínimo de meses pareados na célula x ano
                        flip = TRUE)      # latitude no eixo y
lat$plot                                  # estimado, MIP e diferença; linha preta = média dos anos
lat$data                                  # média por faixa de latitude x ano x série
ggplot2::ggsave("figs/flux_latitudinal_profile.png", lat$plot,
                width = 12, height = 6, dpi = 300, bg = "white")


############################################################
# 5.1 Compare fluxes estimated with MIP and forme beta     #
############################################################

# dfall   <- readxl::read_xlsx("data/beta_significant.xlsx")   # ou o dfall do seu loop
# res <-  compare_beta(dfall, df_comp, biomes = biomas, maps = TRUE,
#                      map_args = list(south_america = south_america, map_theme = map_theme),
#                      unit='repo')
# res$metrics     # Brasil e biomas, três pares
# res$plots       # scatter_brazil, scatter_beta_ref, scatter_beta_est, scatter_est_ref
#
#
# res$maps
#
# bal <- balance_beta(res)            # aceita a saída da compare_beta (ou o data frame res$pairs)
# bal$year                            # Tg C por grupo × ano (beta, estimado, net_flux)
# bal$total                           # Tg C por grupo (soma de todas as células e anos)
# bal$plots$balance_year              # 3 barras por ano, uma facet por Brasil/bioma
# bal$plots$balance_total             # 3 barras por grupo
# save_beta_balance(bal, "figs_beta/")
#
# save_beta_plots(res, "figs_beta/")
# save_beta_maps(res, "figs_beta_maps/")


############################################################
# 6.  Compare fluxes estimated and from MIP with climate   #
############################################################

res_climatica <- compare_climate(df_comp,
                       biomes = biomas,
                       time_scale = "annual",
                       data_dir = "data-raw/nasa_power",
                       download = FALSE,
                       out_dir = "data/nasa_power",   # salva o POWER mensal por célula, antes do match com df_comp
                       map_args = list(
                         south_america = south_america,
                         map_theme = map_theme
                         )
                       )

res_climatica$summary                # mediana de r por bioma, variável e fonte (estimado x MIP)
res_climatica$pooled                 # r, rho, inclinação do scatter geral por bioma
res_climatica$plots$prec             # scatter: fluxo x anomalia da precipitação, por bioma
res_climatica$maps$temp              # mapa de r (estimado e MIP lado a lado) com a temperatura
save_climate_plots(res_climatica, "figs_climate/")


############################################################
# 7.  Fire (INPE AQ1KM) and deforestation (PRODES)         #
############################################################

# 7.1 Área queimada mensal do INPE (AQ1KM, Coleção 2, 1 km, MODIS): download dos GeoTIFFs mensais.
#     YOU DON'T NEED TO RUN THIS EVERY TIME: os arquivos ficam em data-raw/aq1km (~90 MB para 2015-2023)
# aq1km_download(dir = "data-raw/aq1km", years = 2015:2023)

# 7.2 Agregação para a grade de 1° do fluxo (só pixels dentro dos biomas). `annual` usa a união dos meses
fire <- aq1km_aggregate("data-raw/aq1km",
                        years    = 2015:2023,
                        biomes   = biomas,
                        out_dir  = "data/fire")
# fire <- readRDS("data/fire/aq1km_cell_month.rds")
fire$annual         # lon, lat, year, burned_km2 (área queimada ao menos uma vez no ano)
fire$total_biome    # km2 queimados por bioma x ano (conferir contra outras fontes)

# 7.3 Desmatamento anual do PRODES (TerraBrasilis, WFS): baixa os polígonos em páginas e já agrega para 1°
#     (os polígonos não são guardados; cada bioma x ano fica em data/prodes). `year` = ano PRODES (ago-jul)
#     YOU DON'T NEED TO RUN THIS EVERY TIME: ~1 GB de transferência, retoma de onde parou
# prodes <- prodes_aggregate(years = 2015:2023, out_dir = "data/prodes")
prodes <- readRDS("data/prodes/prodes_cell_year.rds")
prodes$cell         # biome, lon, lat, year, def_km2, n_poly
prodes$total        # km2 desmatados por bioma x ano (conferir contra os números oficiais do PRODES)

# Biomas x anos com dado no PRODES: Caatinga, Mata Atlântica, Pampa e Pantanal só têm série a partir de 2016
# (em 2015 vêm vazios, o Pampa com 1 polígono), então 2015 conta como "sem dado" neles, e não como zero
prodes_cov <- prodes$total |>
  dplyr::filter(n_poly > 0,
                !(biome %in% c("Caatinga", "Mata Atlântica", "Pampa", "Pantanal") & year == 2015)) |>
  dplyr::select(biome, year)

# 7.4 Mapas anuais de área queimada e de desmatamento nas células de 1° (escala em Mha; 1 Mha = 10 mil km2)
m_fire <- map_disturbance_annual(fire$annual,
                                 value = "burned_km2", label = "Burned area",
                                 years = 2015:2023,
                                 cells = fire$land,              # células sem registro de fogo entram como 0
                                 biomes = biomas,
                                 south_america = south_america,
                                 map_theme = map_theme)
m_fire$maps$annual     # um painel por ano, mesma escala (0 até o percentil 99, em Mha)
m_fire$maps$total      # acumulado 2015-2023
save_disturbance_maps(m_fire, "figs_disturbance/", prefix = "fire")

m_def <- map_disturbance_annual(prodes$cell,
                                value = "def_km2", label = "Deforestation",
                                years = 2015:2023,
                                cells = fire$land,              # mesma grade de células (dentro dos biomas)
                                coverage = prodes_cov,   # biomas x anos com dado no PRODES
                                biomes = biomas,
                                south_america = south_america,
                                map_theme = map_theme)
m_def$maps$annual
save_disturbance_maps(m_def, "figs_disturbance/", prefix = "deforestation")

# acumulado só de 2016-2023, para todos os biomas terem os mesmos anos
m_def_total <- map_disturbance_annual(prodes$cell,
                                      value = "def_km2", label = "Deforestation",
                                      years = 2016:2023, cells = fire$land, coverage = prodes_cov,
                                      biomes = biomas, south_america = south_america, map_theme = map_theme)
m_def_total$maps$total
ggplot2::ggsave("figs_disturbance/deforestation_total_2016_2023.png", m_def_total$maps$total,
                width = 8.4, height = 7.2, dpi = 300, bg = "white")

