suppressMessages({library(dplyr)}); options(width = 220, digits = 3)
setwd("D:/R/xco2_flux_br")
for (f in list.files("r", full.names = TRUE)) if (!grepl("prodes_aggregate|aq1km_download|setup_pack", f)) try(source(f), silent = TRUE)
biomas <- geobr::read_biomes(year = 2019) |> filter(name_biome != "Sistema Costeiro")
df_flux <- readr::read_rds(list.files("data/Fluxes_MIP", full.names = TRUE))
df_std <- readr::read_rds(list.files("data/Fluxes_MIP_std", full.names = TRUE)) |> dplyr::transmute(lon = round(lon, 4), lat = round(lat, 4), year, month, net_flux_sd = net_flux)
biome_of <- function(xy) { s2 <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(s2)); b <- sf::st_make_valid(sf::st_transform(biomas, 4326))
  hit <- sf::st_intersects(sf::st_as_sf(xy[, c("lon", "lat")], coords = c("lon", "lat"), crs = 4326), b); idx <- vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
  as.character(sf::st_drop_geometry(b)$name_biome)[idx] }
# df: sondagens com delta e sigma ja calculados
run_variant <- function(df, ref_flux = df_flux, flux_df_std = df_std) {
  df$biome <- biome_of(df); df$flux_fraction <- 2 * ifelse(df$biome %in% "Amazônia", 0.05, 0.18)
  d <- fco2_from_delta(df, delta = "delta", sigma_delta = "sigma", method = "massa", tau_days = 2, flux_fraction = df$flux_fraction)
  dc <- agg_xco2_flux(d, ref_flux, vars = "fco2", method = "dois_estagios", sigma = "sigma_f", min_days = 1, join = "left")
  dc <- dplyr::left_join(dplyr::mutate(dc |> filter(fco2_mean > -15), lon = round(lon, 4), lat = round(lat, 4)), flux_df_std, by = c("lon", "lat", "year", "month"))
  res <- compare_flux(dc, flux_all = ref_flux, biomes = biomas, min_days = 1, years = 2015:2023, ref_factor = 1, ref_sigma = "net_flux_sd", sigma_bars = "indep", plot = FALSE)
  res$n_sound <- nrow(df); res
}
summ <- function(res, label) {
  by <- res$balance_year; fm <- res$flux_monthly; m <- res$metrics
  a <- by |> group_by(group) |> summarise(r_ann = cor(est_TgC, ref_TgC), rmse_ann = sqrt(mean((est_TgC - ref_TgC)^2)), diff_ann = mean(est_TgC - ref_TgC), est_ann = mean(est_TgC), .groups = "drop")
  b <- fm |> group_by(group) |> summarise(r_mon = cor(est_flux, ref_flux), .groups = "drop")
  out <- a |> left_join(b, by = "group") |> left_join(m |> select(group, n_cellmonths, rho_cm = rho), by = "group"); out$variant <- label; out$n_sound <- res$n_sound; out
}
