#' Comparação modelo a modelo do balanço anual: estimado (delta) x net_flux de cada grupo do OCO-2 MIP
#'
#' Para cada modelo (arquivo `<Modelo>__<Experimento>.rds` em `model_dir`) usa as mesmas células x mês pareadas
#' do estimado (res$cells, saída de compare_flux), exige net_flux finito do modelo na célula x mês e soma
#' fluxo x área x dias do mês em Tg C, por bioma e ano. O Brasil é a soma dos biomas.
#'
#' @param cells      `res$cells` de compare_flux(): lon, lat, year, month, biome, est (g C m-2 d-1), area_m2...
#' @param model_dir  pasta com os rds por modelo
#' @param experiment experimento (padrão "LNLGIS"); arquivos EnsMean/EnsStd são ignorados
#' @param min_cells  mínimo de células pareadas no bioma-ano para o ponto entrar
#' @param km_per_degree km por grau (igual ao compare_flux)
#' @return lista: balance (modelo, grupo, ano, est_TgC, mod_TgC, n_cells), metrics (r, p, rho, RMSE, bias por modelo)
compare_models <- function(cells, model_dir = "data/Fluxes_MIP_completo", experiment = "LNLGIS",
                           min_cells = 1, km_per_degree = 110) {
  fs <- list.files(model_dir, paste0("__", experiment, "\\.rds$"), full.names = TRUE)
  fs <- fs[!grepl("EnsMean|EnsStd", fs)]
  cells$lon <- round(cells$lon, 4); cells$lat <- round(cells$lat, 4)
  ndays <- function(y, m) as.integer(as.Date(sprintf("%d-%02d-01", ifelse(m == 12, y + 1, y), ifelse(m == 12, 1, m + 1))) -
                                       as.Date(sprintf("%d-%02d-01", y, m)))
  w_area <- (km_per_degree * 1e3)^2 / 1e12                           # g -> Tg e m2
  bal <- lapply(fs, function(f) {
    mod <- sub("__.*$", "", basename(f))
    m <- readRDS(f)
    if (!"net_flux" %in% names(m)) m$net_flux <- m$land_flux + m$ocean_flux   # net = land + ocean (sem fóssil)
    m <- m[, c("lon", "lat", "year", "month", "net_flux")]
    m$lon <- round(m$lon, 4); m$lat <- round(m$lat, 4)
    d <- merge(cells[, c("lon", "lat", "year", "month", "biome", "est")], m, by = c("lon", "lat", "year", "month"))
    d <- d[is.finite(d$est) & is.finite(d$net_flux), ]
    d$w <- w_area * ndays(d$year, d$month)
    d$est_TgC <- d$est * d$w; d$mod_TgC <- d$net_flux * d$w
    a <- aggregate(cbind(est_TgC, mod_TgC) ~ biome + year, d, sum)
    a$n_cells <- as.vector(table(paste(d$biome, d$year))[paste(a$biome, a$year)])
    br <- aggregate(cbind(est_TgC, mod_TgC, n_cells) ~ year, a, sum); br$biome <- "Brasil"
    cbind(model = mod, rbind(a, br[names(a)]))
  })
  bal <- do.call(rbind, bal); bal <- bal[bal$n_cells >= min_cells, ]
  names(bal)[names(bal) == "biome"] <- "group"
  mt <- function(d) {
    ct <- suppressWarnings(stats::cor.test(d$mod_TgC, d$est_TgC)); cs <- suppressWarnings(stats::cor.test(d$mod_TgC, d$est_TgC, method = "spearman"))
    data.frame(n = nrow(d), r = unname(ct$estimate), p = ct$p.value, rho = unname(cs$estimate), p_rho = cs$p.value,
               rmse = sqrt(mean((d$est_TgC - d$mod_TgC)^2)), bias = mean(d$est_TgC - d$mod_TgC))
  }
  bio <- bal[bal$group != "Brasil", ]
  list(balance = tibble::as_tibble(bal),
       metrics = do.call(rbind, lapply(split(bio, bio$model), function(d) cbind(model = d$model[1], mt(d)))),
       metrics_brasil = do.call(rbind, lapply(split(bal[bal$group == "Brasil", ], bal$model[bal$group == "Brasil"]),
                                              function(d) cbind(model = d$model[1], mt(d)))))
}

#' Figura em painéis (um por modelo): balanço anual por bioma e ano, modelo (x) x estimado (y), com r e p
#' @param bal saída de compare_models()
plot_models <- function(bal, out = "figs_manuscript/supp_model_by_model.png", ncol = 4, width = 11, height = 11) {
  library(ggplot2)
  b <- bal$balance[bal$balance$group != "Brasil", ]
  en <- c("Amazônia" = "Amazon", "Caatinga" = "Caatinga", "Cerrado" = "Cerrado", "Mata Atlântica" = "Atlantic Forest",
          "Pampa" = "Pampa", "Pantanal" = "Pantanal")
  b$bioma <- factor(unname(en[b$group]), levels = unname(en))
  cores <- c("Amazon" = "#009E73", "Caatinga" = "#E69F00", "Cerrado" = "#D55E00", "Atlantic Forest" = "#0072B2", "Pampa" = "#CC79A7", "Pantanal" = "#56B4E9")
  formas <- c("Amazon" = 16, "Caatinga" = 17, "Cerrado" = 15, "Atlantic Forest" = 18, "Pampa" = 3, "Pantanal" = 4)
  fp <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  mt <- bal$metrics
  mt$lab <- sprintf("r = %.2f, %s\nRMSE = %.0f, bias = %+.0f", mt$r, fp(mt$p), mt$rmse, mt$bias)
  lim <- range(c(b$mod_TgC, b$est_TgC)); lim <- lim + c(-1, 1) * diff(lim) * 0.05
  b$model <- factor(b$model, levels = sort(unique(b$model))); mt$model <- factor(mt$model, levels = levels(b$model))
  p <- ggplot(b, aes(mod_TgC, est_TgC, colour = bioma, shape = bioma)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey35", linewidth = 0.4) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25) +
    geom_point(size = 1.5, stroke = 0.6, alpha = 0.9) +
    geom_text(data = mt, aes(x = lim[1] + 0.02 * diff(lim), y = lim[2] - 0.02 * diff(lim), label = lab), inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.5, family = "Times New Roman") +
    facet_wrap(~model, ncol = ncol) + coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
    scale_colour_manual(values = cores) + scale_shape_manual(values = formas) +
    labs(x = expression("OCO-2 MIP model annual net flux (Tg C yr"^-1*")"), y = expression("Estimated annual flux (Tg C yr"^-1*")"), colour = NULL, shape = NULL) +
    theme_bw(base_size = 9, base_family = "Times New Roman") +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold"))
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  ragg::agg_png(out, width = width, height = height, units = "in", res = 300); print(p); grDevices::dev.off()
  invisible(p)
}

#' Painel (b) em inglês: balanço anual por bioma e ano, ensemble do MIP (x) x estimado (y), com r e p (Pearson, pooled)
#' @param by `balance_year` de compare_flux() (group, year, est_TgC, ref_TgC)
plot_biome_scatter <- function(by, out = "figs_manuscript/supp_biome_scatter.png", width = 90, height = 80) {
  library(ggplot2)
  en <- c("Amazônia" = "Amazon", "Caatinga" = "Caatinga", "Cerrado" = "Cerrado", "Mata Atlântica" = "Atlantic Forest",
          "Pampa" = "Pampa", "Pantanal" = "Pantanal")
  b <- by[by$group != "Brasil", ]; b$bioma <- factor(unname(en[b$group]), levels = unname(en))
  ct <- stats::cor.test(b$ref_TgC, b$est_TgC)
  lab <- sprintf("r = %.2f, %s (n = %d)", ct$estimate, ifelse(ct$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", ct$p.value)), nrow(b))
  cores <- c("Amazon" = "#009E73", "Caatinga" = "#E69F00", "Cerrado" = "#D55E00", "Atlantic Forest" = "#0072B2", "Pampa" = "#CC79A7", "Pantanal" = "#56B4E9")
  formas <- c("Amazon" = 16, "Caatinga" = 17, "Cerrado" = 15, "Atlantic Forest" = 18, "Pampa" = 3, "Pantanal" = 4)
  lim <- c(-62, 140)
  p <- ggplot(b, aes(ref_TgC, est_TgC, colour = bioma, shape = bioma)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey35", linewidth = 0.4) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_point(size = 1.7, stroke = 0.6, alpha = 0.9) +
    annotate("text", x = lim[1] + 4, y = lim[2] - 4, label = lab, hjust = 0, vjust = 1, size = 2.8, family = "Times New Roman") +
    scale_colour_manual(values = cores) + scale_shape_manual(values = formas) + coord_equal(xlim = lim, ylim = lim, expand = FALSE) +
    labs(x = expression("OCO-2 MIP (Tg C yr"^-1*")"), y = expression("Estimated (Tg C yr"^-1*")"), colour = NULL, shape = NULL) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_classic(base_size = 8, base_family = "Times New Roman") +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.25), legend.position = "bottom", legend.text = element_text(size = 7.5),
          legend.key.size = grid::unit(3.2, "mm"), plot.margin = margin(2, 4, 2, 2))
  ragg::agg_png(out, width = width, height = height, units = "mm", res = 600); print(p); grDevices::dev.off()
  print(ct); invisible(p)
}

#' Médias mensais e ciclo sazonal, modelo a modelo (OCO-2 MIP, experimento LNLGIS) x estimado
#'
#' Fluxo médio (g C m-2 d-1) das células pareadas, por modelo, em duas escalas:
#'   monthly  = grupo x ano x mês (série temporal mensal; grupo = cada bioma e o Brasil)
#'   seasonal = grupo x mês do calendário (média dos anos da média mensal; ciclo sazonal)
#' O estimado e cada modelo usam as mesmas células x mês (pares com net_flux finito do modelo).
#' O arquivo EnsMean (se existir na pasta) entra como "EnsMean" e serve de referência/validação.
#'
#' @param cells res$cells de compare_flux()
#' @return lista: monthly, seasonal (model, group, ..., est, mod) e metrics (r, p, rho, RMSE, bias por escala, nível e modelo)
compare_models_monthly <- function(cells, model_dir = "data/Fluxes_MIP_completo", experiment = "LNLGIS") {
  fs <- list.files(model_dir, paste0("__", experiment, "\\.rds$"), full.names = TRUE)
  fs <- fs[!grepl("EnsStd", fs)]
  cells$lon <- round(cells$lon, 4); cells$lat <- round(cells$lat, 4)
  mon <- lapply(fs, function(f) {
    m <- readRDS(f)
    if (!"net_flux" %in% names(m)) m$net_flux <- m$land_flux + m$ocean_flux
    m$lon <- round(m$lon, 4); m$lat <- round(m$lat, 4)
    d <- merge(cells[, c("lon", "lat", "year", "month", "biome", "est")], m[, c("lon", "lat", "year", "month", "net_flux")],
               by = c("lon", "lat", "year", "month"))
    d <- d[is.finite(d$est) & is.finite(d$net_flux), ]
    a <- aggregate(cbind(est, mod = net_flux) ~ biome + year + month, d, mean)
    br <- aggregate(cbind(est, mod = net_flux) ~ year + month, d, mean); br$biome <- "Brasil"
    cbind(model = sub("__.*$", "", basename(f)), rbind(a, br[names(a)]))
  })
  mon <- do.call(rbind, mon); names(mon)[names(mon) == "biome"] <- "group"
  sea <- aggregate(cbind(est, mod) ~ model + group + month, mon, mean)
  mt <- function(d) {
    ct <- suppressWarnings(stats::cor.test(d$mod, d$est)); cs <- suppressWarnings(stats::cor.test(d$mod, d$est, method = "spearman"))
    data.frame(n = nrow(d), r = unname(ct$estimate), p = ct$p.value, rho = unname(cs$estimate), p_rho = cs$p.value,
               rmse = sqrt(mean((d$est - d$mod)^2)), bias = mean(d$est - d$mod))
  }
  met <- function(d, escala) {                                      # biomas juntos (pooled) e Brasil à parte
    bio <- d[d$group != "Brasil", ]
    rbind(do.call(rbind, lapply(split(bio, bio$model), function(x) cbind(escala = escala, nivel = "biomas (pooled)", model = x$model[1], mt(x)))),
          do.call(rbind, lapply(split(d[d$group == "Brasil", ], d$model[d$group == "Brasil"]),
                                function(x) cbind(escala = escala, nivel = "Brasil", model = x$model[1], mt(x)))))
  }
  list(monthly = tibble::as_tibble(mon), seasonal = tibble::as_tibble(sea),
       metrics = tibble::as_tibble(rbind(met(mon, "mensal"), met(sea, "sazonal"))))
}

.pt_en <- c("Amazônia" = "Amazon", "Caatinga" = "Caatinga", "Cerrado" = "Cerrado", "Mata Atlântica" = "Atlantic Forest",
            "Pampa" = "Pampa", "Pantanal" = "Pantanal")
.biome_cols <- c("Amazon" = "#009E73", "Caatinga" = "#E69F00", "Cerrado" = "#D55E00", "Atlantic Forest" = "#0072B2", "Pampa" = "#CC79A7", "Pantanal" = "#56B4E9")
.biome_shp <- c("Amazon" = 16, "Caatinga" = 17, "Cerrado" = 15, "Atlantic Forest" = 18, "Pampa" = 3, "Pantanal" = 4)

#' Dispersão por modelo (como plot_models), em g C m-2 d-1, para a escala mensal ou sazonal
#' @param d `monthly` ou `seasonal` de compare_models_monthly()
#' @param met métricas já filtradas (escala e nível "biomas (pooled)")
plot_models_flux <- function(d, met, out, ncol = 4, width = 11, height = 11, pt_size = 1.1, clip = TRUE) {
  library(ggplot2)
  d <- d[d$group != "Brasil" & d$model != "EnsMean", ]; met <- met[met$model != "EnsMean", ]
  d$bioma <- factor(unname(.pt_en[d$group]), levels = unname(.pt_en))
  fp <- function(p) ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))
  met$lab <- sprintf("r = %.2f, %s\nRMSE = %.2f, bias = %+.2f", met$r, fp(met$p), met$rmse, met$bias)
  lim <- if (clip) stats::quantile(c(d$mod, d$est), c(0.005, 0.995)) else range(c(d$mod, d$est))   # clip: poucos extremos esticam o eixo
  lim <- lim + c(-1, 1) * diff(lim) * 0.05
  n_out <- sum(d$mod < lim[1] | d$mod > lim[2] | d$est < lim[1] | d$est > lim[2])
  if (n_out > 0) message(sprintf("%d ponto(s) fora dos limites do eixo (incluídos nas métricas): %s", n_out, basename(out)))
  p <- ggplot(d, aes(mod, est, colour = bioma, shape = bioma)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey35", linewidth = 0.4) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25) +
    geom_point(size = pt_size, stroke = 0.5, alpha = 0.8) +
    geom_text(data = met, aes(x = lim[1] + 0.02 * diff(lim), y = lim[2] - 0.02 * diff(lim), label = lab), inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 2.5, family = "Times New Roman") +
    facet_wrap(~model, ncol = ncol) + coord_equal(xlim = lim, ylim = lim, expand = FALSE, clip = "on") +
    scale_colour_manual(values = .biome_cols) + scale_shape_manual(values = .biome_shp) +
    labs(x = expression("OCO-2 MIP model net flux (g C m"^-2*" d"^-1*")"), y = expression("Estimated flux (g C m"^-2*" d"^-1*")"), colour = NULL, shape = NULL) +
    theme_bw(base_size = 9, base_family = "Times New Roman") +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold"))
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  ragg::agg_png(out, width = width, height = height, units = "in", res = 300); print(p); grDevices::dev.off()
  invisible(p)
}

#' Ciclo sazonal por bioma: estimado, média do ensemble do MIP e faixa min-max dos modelos
plot_seasonal_cycle <- function(sea, out = "figs_manuscript/supp_seasonal_cycle.png", width = 10, height = 6) {
  library(ggplot2)
  ens <- sea[sea$model == "EnsMean", ]; mods <- sea[sea$model != "EnsMean", ]
  rng <- aggregate(mod ~ group + month, mods, range)
  rng <- data.frame(group = rng$group, month = rng$month, lo = rng$mod[, 1], hi = rng$mod[, 2])
  lv <- c("Brasil", names(.pt_en)); nm <- c("Brazil", unname(.pt_en))
  fac <- function(x) factor(nm[match(x, lv)], levels = nm)
  rng$g <- fac(rng$group); ens$g <- fac(ens$group)
  ln <- rbind(data.frame(g = ens$g, month = ens$month, v = ens$est, serie = "Estimated (delta)"),
              data.frame(g = ens$g, month = ens$month, v = ens$mod, serie = "OCO-2 MIP ensemble mean"))
  p <- ggplot() + geom_ribbon(data = rng, aes(month, ymin = lo, ymax = hi, fill = "OCO-2 MIP models (min-max)"), alpha = 0.25) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(data = ln, aes(month, v, colour = serie, linetype = serie), linewidth = 0.7) +
    geom_point(data = ln, aes(month, v, colour = serie), size = 1.3) +
    facet_wrap(~g, ncol = 4, scales = "free_y") + scale_x_continuous(breaks = c(1, 4, 7, 10), labels = c("Jan", "Apr", "Jul", "Oct")) +
    scale_colour_manual(values = c("Estimated (delta)" = "#0B6E6E", "OCO-2 MIP ensemble mean" = "#4D4D4D")) +
    scale_linetype_manual(values = c("solid", "22")) + scale_fill_manual(values = c("OCO-2 MIP models (min-max)" = "grey50")) +
    labs(x = NULL, y = expression("Net flux (g C m"^-2*" d"^-1*")"), colour = NULL, linetype = NULL, fill = NULL) +
    theme_bw(base_size = 9, base_family = "Times New Roman") +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold"))
  ragg::agg_png(out, width = width, height = height, units = "in", res = 300); print(p); grDevices::dev.off()
  invisible(p)
}

