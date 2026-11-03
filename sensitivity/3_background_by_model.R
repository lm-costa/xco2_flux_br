source("C:/Users/luism/AppData/Local/Temp/claude/D--R-xco2-flux-br/aca67c70-3b11-41db-b590-5ec4e323c2af/scratchpad/variant_lib.R")
ens <- readRDS("data/xco2_pre_processed/ensemble/ensemble__IS.rds")
ens$biome <- biome_of(ens); key <- function(x) paste(x$time, round(x$lon, 4), round(x$lat, 4)); ens$k <- key(ens)
bio_cache <- setNames(ens$biome, ens$k)
run_v2 <- function(df, ref_flux) {   # igual run_variant, mas com o bioma ja calculado
  df$flux_fraction <- 2 * ifelse(df$biome %in% "Amazônia", 0.05, 0.18)
  d <- fco2_from_delta(df, delta = "delta", sigma_delta = "sigma", method = "massa", tau_days = 2, flux_fraction = df$flux_fraction)
  dc <- agg_xco2_flux(d, ref_flux, vars = "fco2", method = "dois_estagios", sigma = "sigma_f", min_days = 1, join = "left")
  dc <- dplyr::left_join(dplyr::mutate(dc |> filter(fco2_mean > -15), lon = round(lon, 4), lat = round(lat, 4)), df_std, by = c("lon", "lat", "year", "month"))
  res <- compare_flux(dc, flux_all = ref_flux, biomes = biomas, min_days = 1, years = 2015:2023, ref_factor = 1, ref_sigma = "net_flux_sd", sigma_bars = "indep", plot = FALSE)
  res$n_sound <- nrow(df); res }
models <- sub("__IS[.]rds$", "", list.files("data/xco2_pre_processed", pattern = "__IS[.]rds$"))
fl <- function(m) { f <- file.path("data/Fluxes_MIP_completo", paste0(m, "__LNLGIS.rds")); if (!file.exists(f)) return(NULL); z <- readr::read_rds(f); if (!"net_flux" %in% names(z)) z$net_flux <- z$land_flux + z$ocean_flux; z }
smm <- function(res, m, bg) { by <- res$balance_year; b6 <- by |> filter(group != "Brasil"); fm <- res$flux_monthly
  br <- by |> filter(group == "Brasil"); mo <- fm |> filter(group == "Brasil")
  data.frame(model = m, background = bg, n_sound = res$n_sound,
    r_pooled = cor(b6$est_TgC, b6$ref_TgC), rmse_pooled = sqrt(mean((b6$est_TgC - b6$ref_TgC)^2)), bias_pooled = mean(b6$est_TgC - b6$ref_TgC),
    r_brazil = cor(br$est_TgC, br$ref_TgC), bias_brazil = mean(br$est_TgC - br$ref_TgC), r_month_brazil = cor(mo$est_flux, mo$ref_flux)) }
out <- list()
for (m in models) {
  fx <- fl(m); if (is.null(fx)) { cat(m, ": sem fluxo LNLGIS\n"); next }
  x <- readRDS(file.path("data/xco2_pre_processed", paste0(m, "__IS.rds"))); x$k <- key(x); x <- x[x$k %in% ens$k, ]
  x$biome <- bio_cache[x$k]; e <- ens[match(x$k, ens$k), ]
  own <- x |> mutate(delta = xco2 - xco2_sim, sigma = sqrt(uncertanty^2 + model_error^2))
  ensb <- x |> mutate(delta = xco2 - e$xco2_ens_mean, sigma = sqrt(e$xco2_ens_se^2 + uncertanty^2 + model_error^2))
  r1 <- tryCatch(run_v2(ensb, fx), error = function(err) NULL); r2 <- tryCatch(run_v2(own, fx), error = function(err) NULL)
  if (!is.null(r1) && !is.null(r2)) { out[[m]] <- bind_rows(smm(r1, m, "ensemble mean"), smm(r2, m, "own IS")); cat(m, "ok\n") } else cat(m, "falhou\n")
  saveRDS(bind_rows(out), "C:/Users/luism/AppData/Local/Temp/sens333.rds")
}
S <- bind_rows(out); print(as.data.frame(S |> arrange(model, background)))
w <- S |> tidyr::pivot_wider(names_from = background, values_from = c(r_pooled, rmse_pooled, bias_pooled, r_brazil, bias_brazil, r_month_brazil, n_sound)); 
cat("\nmedias entre modelos:\n"); print(as.data.frame(S |> group_by(background) |> summarise(across(c(r_pooled, rmse_pooled, bias_pooled, r_brazil, bias_brazil, r_month_brazil), mean))))
cat("FIM\n")
