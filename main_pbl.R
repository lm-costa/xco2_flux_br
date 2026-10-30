############################################################
#  Delta com f_a por bioma (profundidade da PBL)           #
############################################################
# Refaz as secoes 3 a 5 e 8 do main.R com o fator f_a da Eq. 3 dependente do bioma, a pedido de D. Baker:
#   f_a = 0.05 na Amazonia e 0.18 nos demais biomas (fracao da coluna de 10 km ocupada pela PBL, da profundidade da PBL do ECMWF
#   por estacao), multiplicado por 2 (o sinal esta concentrado na metade de baixo da coluna), e tau = 2 dias (como no texto do artigo).
#   flux_fraction = 0.10 (Amazonia) e 0.36 (demais), equivalente a f_a = 1 com esses fatores por bioma.
# As figuras vao para figs_pbl/ e figs_manuscript_pbl/, sem sobrescrever as do main.R. Rodar a partir da raiz do projeto.

purrr::map(list.files('r', full.names = TRUE), source)
setup_packages()
library(dplyr)

fa_amazon <- 0.05; fa_other <- 0.18; mult <- 2; tau <- 2

biomas <- geobr::read_biomes(year = 2019) |> filter(name_biome != 'Sistema Costeiro')
south_america <- sf::read_sf(list.files('data-raw/South_America/', pattern = 'shp', full.names = TRUE)[1])
df_flux <- readr::read_rds(list.files('data/Fluxes_MIP', full.names = TRUE))
df_std  <- readr::read_rds(list.files('data/Fluxes_MIP_std', full.names = TRUE)) |>
  dplyr::transmute(lon = round(lon, 4), lat = round(lat, 4), year, month, net_flux_sd = net_flux)
df <- readr::read_rds(list.files('data/xco2_pre_processed/ensemble', full.names = TRUE))

# bioma de cada sondagem (fora dos biomas: usa o valor dos demais)
biome_of <- function(xy) {
  s2 <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(s2))
  b <- sf::st_make_valid(sf::st_transform(biomas, 4326))
  hit <- sf::st_intersects(sf::st_as_sf(xy[, c('lon', 'lat')], coords = c('lon', 'lat'), crs = 4326), b)
  idx <- vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
  as.character(sf::st_drop_geometry(b)$name_biome)[idx]
}

df_n <- df |>
  mutate(delta = xco2 - xco2_ens_mean,
         sigma = sqrt((xco2_ens_se^2) + (uncertanty^2) + (model_error^2)),
         biome = biome_of(pick(lon, lat)),
         flux_fraction = mult * ifelse(biome %in% 'Amazônia', fa_amazon, fa_other))
table(df_n$biome, useNA = 'ifany')

df_n <- fco2_from_delta(df_n, delta = 'delta', sigma_delta = 'sigma', method = 'massa',
                        tau_days = tau, flux_fraction = df_n$flux_fraction)

df_comp <- agg_xco2_flux(df_n, df_flux, vars = 'fco2', method = 'dois_estagios', sigma = 'sigma_f',
                         min_days = 1, join = 'left', out_file = 'data/comparison_pbl.rds')
df_comp <- df_comp |>
  filter(fco2_mean > -15) |>
  mutate(lon = round(lon, 4), lat = round(lat, 4)) |>
  left_join(df_std, by = c('lon', 'lat', 'year', 'month'))

############################################################
# Comparacao com o MIP (LNLGIS) e figuras                  #
############################################################
res <- compare_flux(df_comp, flux_all = df_flux, biomes = biomas, min_days = 1, years = 2015:2023,
                    ref_factor = 1, ref_sigma = 'net_flux_sd', sigma_bars = 'indep')
res$metrics; res$balance_total; res$balance_year
saveRDS(res, 'output/metrics_temporal_pbl.rds')
dir.create('figs_pbl', showWarnings = FALSE); dir.create('figs_manuscript_pbl', showWarnings = FALSE)
save_flux_plots(res, dir = 'figs_pbl/')                           # Fig. 4, Fig. 5 e Supp. Fig. 4

out <- map_flux_annual(df_comp, biomes = biomas, south_america = south_america, map_theme = map_theme,
                       est_sigma = 'sigma_f_cell', ref_sigma = 'net_flux_sd', years = 2015:2023, min_months = 3, min_years = 3)
saveRDS(out, 'output/metrics_spatial_pbl.rds')
save_flux_maps(out, 'figs_pbl/')                                  # Fig. 7, Fig. 8 e Supp. Fig. 5

lat <- plot_lat_profile(out$annual, bin = 2, min_cells = 5, min_months = 3, flip = TRUE)   # Fig. 6
ggplot2::ggsave('figs_pbl/flux_latitudinal_profile.png', lat$plot, width = 12, height = 6, dpi = 300, bg = 'white')

mm <- compare_models(res$cells, 'data/Fluxes_MIP_completo', experiment = 'LNLGIS')            # Supp. Fig. 3
saveRDS(mm, 'output/compare_models_pbl.rds')
plot_models(mm, 'figs_manuscript_pbl/supp_model_by_model.png')
plot_biome_scatter(res$balance_year, 'figs_manuscript_pbl/supp_biome_scatter.png')            # Supp. Fig. 2

mmm <- compare_models_monthly(res$cells, 'data/Fluxes_MIP_completo', experiment = 'LNLGIS')
met <- mmm$metrics
plot_models_flux(mmm$monthly, met[met$escala == 'mensal' & met$nivel == 'biomas (pooled)', ], 'figs_manuscript_pbl/supp_model_monthly.png')
plot_models_flux(mmm$seasonal, met[met$escala == 'sazonal' & met$nivel == 'biomas (pooled)', ], 'figs_manuscript_pbl/supp_model_seasonal.png', pt_size = 1.6)
plot_seasonal_cycle(mmm$seasonal, 'figs_manuscript_pbl/supp_seasonal_cycle.png')

############################################################
# Antes x depois (mesma referencia: net_flux do LNLGIS)    #
############################################################
# "atual" = f_a 0.2 e tau 2.5 d (main.R); "tau = 2" = o mesmo com tau 2 (so multiplica por 2.5/2, exato porque o fluxo e linear em 1/tau)
old <- readRDS('output/metrics_temporal.rds')
cmp_metrics <- function(cells, scale = 1, label) {
  cells$est <- cells$est * scale
  w <- cells$area_m2 * cells$ndays / 1e12
  one <- function(g) {
    i <- if (g == 'Brasil') rep(TRUE, nrow(cells)) else cells$biome == g
    x <- cells$ref[i]; y <- cells$est[i]; yr <- cells$year[i]
    by <- tapply(y * w[i], yr, sum); bx <- tapply(x * w[i], yr, sum)
    data.frame(scheme = label, group = g, n = sum(i), r = stats::cor(x, y), slope = unname(stats::coef(stats::lm(y ~ x))[2]),
               bias = mean(y - x), rmse = sqrt(mean((y - x)^2)), est_TgC = sum(by), ref_TgC = sum(bx), ratio = sum(by) / sum(bx),
               r_annual = if (length(by) >= 3) stats::cor(by, bx) else NA_real_)
  }
  do.call(rbind, lapply(c('Brasil', sort(unique(cells$biome))), one))
}
tab <- rbind(cmp_metrics(old$cells, 1, 'current (f_a 0.2, tau 2.5)'),
             cmp_metrics(old$cells, 2.5 / 2, 'tau 2 only'),
             cmp_metrics(res$cells, 1, 'PBL f_a x2, tau 2'))
print(tab, digits = 3, row.names = FALSE)
readr::write_csv(tab, 'output/pbl_vs_current.csv')
