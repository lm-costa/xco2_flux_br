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

