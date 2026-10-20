#' Incremento do OCO-2 nos fluxos do MIP: LNLGIS - IS
#'
#' A estimativa por delta usa como fundo o Xco2 do experimento IS (so in situ), entao representa o que o OCO-2 acrescenta
#' a rede in situ. O net_flux do LNLGIS contem esse incremento e tambem a parte ja restrita pelo in situ. A comparacao
#' "lado a lado" usa, portanto, a diferenca LNLGIS - IS, calculada grupo a grupo (mesmo modelo de transporte e mesmo prior
#' nos dois experimentos) e depois resumida pela media e pelo desvio entre os grupos que tem os dois experimentos.
#'
#' @param model_dir pasta com os rds por grupo e experimento (`<Grupo>__<Experimento>.rds`)
#' @param exp_a,exp_b experimentos (A - B)
#' @param flux coluna de fluxo (padrao "net_flux"; se nao existir, usa land_flux + ocean_flux)
#' @return tibble: lon, lat, year, month, inc_mean, inc_sd (entre grupos), n_models, a_mean, b_mean
mip_increment <- function(model_dir = "data/Fluxes_MIP_completo", exp_a = "LNLGIS", exp_b = "IS", flux = "net_flux") {
  lf <- function(exp) {
    f <- list.files(model_dir, paste0("__", exp, "\\.rds$"), full.names = TRUE)
    f[!grepl("EnsMean|EnsStd", f)]
  }
  grp <- function(f) sub("__.*$", "", basename(f))
  fa <- lf(exp_a); fb <- lf(exp_b); common <- intersect(grp(fa), grp(fb))
  if (!length(common)) stop("nenhum grupo com os dois experimentos")
  rd <- function(f) {
    d <- readRDS(f)
    if (!flux %in% names(d)) d[[flux]] <- d$land_flux + d$ocean_flux
    d <- d[, c("lon", "lat", "year", "month", flux)]; d$lon <- round(d$lon, 4); d$lat <- round(d$lat, 4); names(d)[5] <- "v"; d
  }
  one <- lapply(common, function(g) {
    a <- rd(fa[grp(fa) == g]); b <- rd(fb[grp(fb) == g])
    m <- merge(a, b, by = c("lon", "lat", "year", "month"), suffixes = c("_a", "_b"))
    m$inc <- m$v_a - m$v_b; m$group <- g; m
  })
  all <- do.call(rbind, one)
  all <- all[is.finite(all$inc), ]
  agg <- stats::aggregate(cbind(inc, v_a, v_b) ~ lon + lat + year + month, all, mean)
  sdv <- stats::aggregate(inc ~ lon + lat + year + month, all, stats::sd)
  nn <- stats::aggregate(inc ~ lon + lat + year + month, all, length)
  out <- merge(merge(agg, sdv, by = c("lon", "lat", "year", "month"), suffixes = c("_mean", "_sd")), nn, by = c("lon", "lat", "year", "month"))
  names(out) <- c("lon", "lat", "year", "month", "inc_mean", "a_mean", "b_mean", "inc_sd", "n_models")
  message(sprintf("incremento %s - %s: %d grupos (%s)", exp_a, exp_b, length(common), paste(common, collapse = ", ")))
  tibble::as_tibble(out)
}

#' Compara a estimativa por delta com o incremento do OCO-2 (LNLGIS - IS) e, para contraste, com o net_flux do LNLGIS
#' Usa compare_flux() nas duas referencias, sobre as mesmas celulas x mes.
#' @param df_comp saida de agg_xco2_flux() (com fco2_mean, sigma_f_cell etc.), ja filtrada como no main
#' @param inc saida de mip_increment()
#' @param full_ref net_flux do LNLGIS (ensemble) por celula x mes para o contraste (coluna net_flux)
#' @return lista: inc (resultado de compare_flux com o incremento), total (com o net_flux do LNLGIS) e table (resumo)
compare_increment <- function(df_comp, inc, full_ref, biomes, years = 2015:2023, ...) {
  d <- merge(df_comp, inc[, c("lon", "lat", "year", "month", "inc_mean", "inc_sd")], by = c("lon", "lat", "year", "month"))
  d <- d[is.finite(d$inc_mean), ]
  d_inc <- d; d_inc$net_flux <- d$inc_mean; d_inc$net_flux_sd <- d$inc_sd
  r_inc <- compare_flux(d_inc, biomes = biomes, years = years, ref_sigma = "net_flux_sd", sigma_bars = "indep", plot = FALSE, ...)
  d_tot <- d; d_tot$net_flux <- NULL
  d_tot <- merge(d_tot, full_ref[, c("lon", "lat", "year", "month", "net_flux")], by = c("lon", "lat", "year", "month"))
  r_tot <- compare_flux(d_tot, biomes = biomes, years = years, sigma_bars = "indep", plot = FALSE, ...)
  pick <- function(r, nm) {
    by <- r$balance_year; m <- r$metrics
    data.frame(ref = nm, m)
  }
  list(inc = r_inc, total = r_tot, n_cellmonths = c(inc = nrow(d_inc), total = nrow(d_tot)))
}
