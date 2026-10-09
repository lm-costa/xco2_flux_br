purrr::map(list.files('r',full.names = T),source)
south_file <- list.files('data-raw/South_America/',pattern = 'shp',full.names = T)
south_america <- sf::read_sf(south_file[1])



mip_extractor(
  nc_file = 'data-raw/OCO2.nc',
  out_dir = 'data/xco2_pre_processed',
  geometry = geobr::read_country(year = 2025),

)

mip_ensemble(
  'data/xco2_pre_processed',
  'IS'
)


#####
library(tidyverse)

df <- readr::read_rds(
  list.files('data/xco2_pre_processed/ensemble',
             full.names = T)
)


df |>
  sample_n(1000) |>
  ggplot(
    aes(x=lon,y=lat,col=xco2)
  )+
  geom_point()


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

#####

df_n <- fco2_from_delta(
  df_n,
  delta='delta',
  sigma_delta = 'sigma',
  method = 'massa',
  tau_days = 30
)

df_n

#####

df_flux <- flux_extractor(
  in_dir      = 'data-raw/OCO2_v11MIP_gridded_fluxes_all_20260729.v2r3',
  out_dir     = "data/Fluxes_MIP",
  geometry    = geobr::read_country(year = 2025),
  submissions = "EnsMean",
  experiments = "LNLGIS"
)

df_flux_std <- flux_extractor('data-raw/OCO2_v11MIP_gridded_fluxes_all_20260729.v2r3',
                         "data/fluxes_MIP_std",
                         geometry = geobr::read_country(year = 2025),
                         submissions = "EnsStd",
                         experiments = "LNLGIS")

df_flux <- readr::read_rds(df_flux)

df_std <- readRDS(df_flux_std) |>
  dplyr::transmute(lon = round(lon, 4), lat = round(lat, 4),
                   year, month, net_flux_sd = net_flux)

df_flux |>
  filter(month==1) |>
  ggplot(
    aes(x=lon,y=lat,fill=net_flux)
  )+
  geom_raster()




####

df_comp <- agg_xco2_flux(
  df_n,
  df_flux,
  vars = 'fco2',
  method = 'dois_estagios',
  sigma = 'sigma_f',
  min_days = 1,
  join='left',
  out_file = 'data/comparison.rds',

)

df_comp <- dplyr::left_join(
  dplyr::mutate(df_comp |>
                  filter(
                    fco2_mean > -15

                    ), lon = round(lon, 4), lat = round(lat, 4)),
  df_std, by = c("lon", "lat", "year", "month")
  )


#####
biomas   <- geobr::read_biomes(year = 2019)
biomas <- biomas |>
  filter(
    name_biome!='Sistema Costeiro'
  )


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
res$flux_monthly     # fluxo médio mensal (g C m-2 d-1) por grupo
res$plots$scatter; res$plots$balance; res$plots$monthly;res$plots$balance_year   # por bioma, em inglês


save_flux_plots(res, dir = "figs/")                       # salva PNGs

####


out <- map_flux_annual(df_comp, biomes = biomas,
                       south_america = south_america,   # seu objeto sf
                       map_theme = map_theme,         # seu tema
                       years = 2015:2023,
                       min_months = 3,   # meses pareados mínimos na célula x ano
                       min_years  = 3)   # anos mínimos para a regressão da célula

out$summary         # Brasil: medianas de r, slope, rmse, bias; fração de células com p < 0.05
out$summary_biome   # o mesmo por bioma
out$regression      # métricas célula a célula (lon, lat, slope, r, p, rmse, bias...)
out$maps$flux_estimated; out$maps$flux_reference; out$maps$bias; out$maps$r; out$maps$slope
out$maps$rmse

save_flux_maps(out, "figs/")


##### BETA

dfall   <- readxl::read_xlsx("data/beta_significant.xlsx")   # ou o dfall do seu loop
res <-  compare_beta(dfall, df_comp, biomes = biomas, maps = TRUE,
                     map_args = list(south_america = south_america, map_theme = map_theme),
                     unit='repo')
res$metrics     # Brasil e biomas, três pares
res$plots       # scatter_brazil, scatter_beta_ref, scatter_beta_est, scatter_est_ref


res$maps

bal <- balance_beta(res)            # aceita a saída da compare_beta (ou o data frame res$pairs)
bal$year                            # Tg C por grupo × ano (beta, estimado, net_flux)
bal$total                           # Tg C por grupo (soma de todas as células e anos)
bal$plots$balance_year              # 3 barras por ano, uma facet por Brasil/bioma
bal$plots$balance_total             # 3 barras por grupo
save_beta_balance(bal, "figs_beta/")

save_beta_plots(res, "figs_beta/")
save_beta_maps(res, "figs_beta_maps/")
