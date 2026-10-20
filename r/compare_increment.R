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

#' Figuras da comparacao com o incremento (LNLGIS - IS) e com o total (LNLGIS)
#' @param cmp saida de compare_increment() (precisa de $inc e $total com balance_year, flux_monthly e cells)
#' @param out_dir pasta de saida
plot_increment <- function(cmp, out_dir = "figs_increment") {
  library(ggplot2); dir.create(out_dir, showWarnings = FALSE)
  en <- c("Brasil" = "Brazil", "Amazônia" = "Amazon", "Caatinga" = "Caatinga", "Cerrado" = "Cerrado", "Mata Atlântica" = "Atlantic Forest", "Pampa" = "Pampa", "Pantanal" = "Pantanal")
  lab_g <- function(x) factor(unname(en[x]), levels = unname(en))
  nm <- c(est = "Estimated (delta)", inc = "LNLGIS − IS (OCO-2 increment)", tot = "LNLGIS (total)")
  cols <- setNames(c("#0B6E6E", "#D55E00", "#4D4D4D"), nm); shp <- setNames(c(16, 17, 15), nm); lty <- setNames(c("solid", "solid", "22"), nm)
  th <- theme_bw(base_size = 10, base_family = "Times New Roman") + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold"))
  save <- function(p, f, w, h) { ragg::agg_png(file.path(out_dir, f), width = w, height = h, units = "in", res = 300); print(p); grDevices::dev.off() }
  fp <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  rcor <- function(x, y) { ct <- stats::cor.test(x, y); list(r = unname(ct$estimate), p = ct$p.value, rmse = sqrt(mean((y - x)^2)), bias = mean(y - x)) }

  # ---- 1. balanco anual: series
  by <- merge(cmp$inc$balance_year[, c("group", "year", "est_TgC", "ref_TgC")], cmp$total$balance_year[, c("group", "year", "ref_TgC")], by = c("group", "year"), suffixes = c("_inc", "_tot"))
  l1 <- rbind(data.frame(g = by$group, year = by$year, v = by$est_TgC, s = nm["est"]), data.frame(g = by$group, year = by$year, v = by$ref_TgC_inc, s = nm["inc"]),
              data.frame(g = by$group, year = by$year, v = by$ref_TgC_tot, s = nm["tot"]))
  l1$g <- lab_g(l1$g); l1$s <- factor(l1$s, levels = nm)
  p1 <- ggplot(l1, aes(year, v, colour = s, shape = s, linetype = s)) + geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.3) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.8) + facet_wrap(~g, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = cols) + scale_shape_manual(values = shp) + scale_linetype_manual(values = lty) + scale_x_continuous(breaks = c(2015, 2019, 2023)) +
    labs(x = NULL, y = expression("Annual carbon balance (Tg C yr"^-1*")"), colour = NULL, shape = NULL, linetype = NULL) + guides(colour = guide_legend(nrow = 1)) + th
  save(p1, "inc_annual_series.png", 11, 6)

  # ---- 2. balanco anual: dispersao por grupo
  s2 <- do.call(rbind, lapply(split(by, by$group), function(d) { m <- rcor(d$ref_TgC_inc, d$est_TgC); m2 <- rcor(d$ref_TgC_tot, d$est_TgC)
    data.frame(g = d$group[1], lab = sprintf("increment: r = %.2f, %s\ntotal: r = %.2f, %s\nRMSE = %.0f / %.0f", m$r, fp(m$p), m2$r, fp(m2$p), m$rmse, m2$rmse)) }))
  d2 <- rbind(data.frame(g = by$group, x = by$ref_TgC_inc, y = by$est_TgC, s = nm["inc"]), data.frame(g = by$group, x = by$ref_TgC_tot, y = by$est_TgC, s = nm["tot"]))
  d2$g <- lab_g(d2$g); s2$g <- lab_g(s2$g); d2$s <- factor(d2$s, levels = nm)
  p2 <- ggplot(d2, aes(x, y, colour = s, shape = s)) + geom_abline(linetype = "22", colour = "grey35") + geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) + geom_point(size = 2, alpha = 0.9) +
    geom_text(data = s2, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -0.05, vjust = 1.2, size = 2.5, family = "Times New Roman") +
    facet_wrap(~g, scales = "free", ncol = 4) + scale_colour_manual(values = cols[c(2, 3)]) + scale_shape_manual(values = shp[c(2, 3)]) +
    labs(x = expression("MIP annual balance (Tg C yr"^-1*")"), y = expression("Estimated annual balance (Tg C yr"^-1*")"), colour = NULL, shape = NULL) + th
  save(p2, "inc_annual_scatter.png", 11, 6.5)

  # ---- 3. media mensal: series
  fm <- merge(cmp$inc$flux_monthly[, c("group", "year", "month", "est_flux", "ref_flux")], cmp$total$flux_monthly[, c("group", "year", "month", "ref_flux")], by = c("group", "year", "month"), suffixes = c("_inc", "_tot"))
  fm$date <- as.Date(sprintf("%d-%02d-15", fm$year, fm$month))
  l3 <- rbind(data.frame(g = fm$group, date = fm$date, v = fm$est_flux, s = nm["est"]), data.frame(g = fm$group, date = fm$date, v = fm$ref_flux_inc, s = nm["inc"]),
              data.frame(g = fm$group, date = fm$date, v = fm$ref_flux_tot, s = nm["tot"]))
  l3$g <- lab_g(l3$g); l3$s <- factor(l3$s, levels = nm)
  m3 <- do.call(rbind, lapply(split(fm, fm$group), function(d) { a <- rcor(d$ref_flux_inc, d$est_flux); b <- rcor(d$ref_flux_tot, d$est_flux)
    data.frame(g = d$group[1], lab = sprintf("r (increment) = %.2f; r (total) = %.2f", a$r, b$r)) })); m3$g <- lab_g(m3$g)
  p3 <- ggplot(l3, aes(date, v, colour = s, linetype = s)) + geom_hline(yintercept = 0, colour = "grey65", linewidth = 0.3) + geom_line(linewidth = 0.45) +
    geom_text(data = m3, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -0.03, vjust = 1.3, size = 2.5, family = "Times New Roman") +
    facet_wrap(~g, scales = "free_y", ncol = 2) + scale_colour_manual(values = cols) + scale_linetype_manual(values = lty) +
    labs(x = NULL, y = expression("Monthly mean flux (g C m"^-2*" d"^-1*")"), colour = NULL, linetype = NULL) + th
  save(p3, "inc_monthly_series.png", 10, 9)

  # ---- 4. media mensal: dispersao por grupo
  d4 <- rbind(data.frame(g = fm$group, x = fm$ref_flux_inc, y = fm$est_flux, s = nm["inc"]), data.frame(g = fm$group, x = fm$ref_flux_tot, y = fm$est_flux, s = nm["tot"]))
  d4$g <- lab_g(d4$g); d4$s <- factor(d4$s, levels = nm)
  s4 <- do.call(rbind, lapply(split(fm, fm$group), function(d) { a <- rcor(d$ref_flux_inc, d$est_flux)
    data.frame(g = d$group[1], lab = sprintf("increment: r = %.2f, %s, RMSE = %.2f", a$r, fp(a$p), a$rmse)) })); s4$g <- lab_g(s4$g)
  p4 <- ggplot(d4, aes(x, y, colour = s, shape = s)) + geom_abline(linetype = "22", colour = "grey35") + geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) + geom_point(size = 1.3, alpha = 0.7) +
    geom_text(data = s4, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE, hjust = -0.03, vjust = 1.3, size = 2.5, family = "Times New Roman") +
    facet_wrap(~g, ncol = 4) + scale_colour_manual(values = cols[c(2, 3)]) + scale_shape_manual(values = shp[c(2, 3)]) +
    labs(x = expression("MIP monthly mean flux (g C m"^-2*" d"^-1*")"), y = expression("Estimated monthly mean flux (g C m"^-2*" d"^-1*")"), colour = NULL, shape = NULL) + th
  save(p4, "inc_monthly_scatter.png", 11, 6.5)

  # ---- 5. pares celula x mes (todos), por grupo
  ci <- cmp$inc$cells; ci$g <- ci$biome; cb <- rbind(ci, transform(ci, g = "Brasil"))
  s5 <- do.call(rbind, lapply(split(cb, cb$g), function(d) { a <- stats::cor(d$ref, d$est, method = "spearman")
    data.frame(g = d$g[1], lab = sprintf("rho = %.2f, RMSE = %.2f, bias = %+.2f", a, sqrt(mean((d$est - d$ref)^2)), mean(d$est - d$ref))) })); s5$g <- lab_g(s5$g)
  cb$g <- lab_g(cb$g); lim <- stats::quantile(c(cb$ref, cb$est), c(0.005, 0.995))
  p5 <- ggplot(cb, aes(ref, est)) + geom_abline(linetype = "22", colour = "grey35") + geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.25) +
    geom_point(colour = "#D55E00", alpha = 0.15, size = 0.6) + geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.5, formula = y ~ x) +
    geom_text(data = s5, aes(x = lim[1], y = lim[2], label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1.2, size = 2.5, family = "Times New Roman") +
    facet_wrap(~g, ncol = 4) + coord_cartesian(xlim = lim, ylim = lim) +
    labs(x = expression("LNLGIS − IS (g C m"^-2*" d"^-1*")"), y = expression("Estimated flux (g C m"^-2*" d"^-1*")")) + th
  save(p5, "inc_cellmonth_scatter.png", 11, 6.5)

  # ---- 6. resumo: r por grupo, incremento x total, anual e mensal
  rr <- function(x, f, y, scale) do.call(rbind, lapply(split(x, x$group), function(d) data.frame(g = d$group[1], scale = scale, ref = f, r = rcor(d[[y[1]]], d[[y[2]]])$r)))
  a1 <- rr(by, "LNLGIS − IS", c("ref_TgC_inc", "est_TgC"), "Annual balance"); a2 <- rr(by, "LNLGIS (total)", c("ref_TgC_tot", "est_TgC"), "Annual balance")
  b1 <- rr(fm, "LNLGIS − IS", c("ref_flux_inc", "est_flux"), "Monthly means"); b2 <- rr(fm, "LNLGIS (total)", c("ref_flux_tot", "est_flux"), "Monthly means")
  d6 <- rbind(a1, a2, b1, b2); d6$g <- lab_g(d6$g); d6$ref <- factor(d6$ref, levels = c("LNLGIS (total)", "LNLGIS − IS"))
  p6 <- ggplot(d6, aes(g, r, fill = ref)) + geom_hline(yintercept = 0, colour = "grey50") + geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_text(aes(label = sprintf("%.2f", r), y = r + ifelse(r >= 0, 0.05, -0.07)), position = position_dodge(0.8), size = 2.4, family = "Times New Roman") +
    facet_wrap(~scale, ncol = 1) + scale_fill_manual(values = c("LNLGIS (total)" = "#4D4D4D", "LNLGIS − IS" = "#D55E00")) +
    labs(x = NULL, y = "Pearson correlation (r) with the estimate", fill = "Reference") + th + theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save(p6, "inc_r_summary.png", 8, 7)
  invisible(list(p1, p2, p3, p4, p5, p6))
}
