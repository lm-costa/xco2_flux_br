#' Compara o beta do seu repositório (XCO2_derivative_flux) com o net_flux do MIP e com o fluxo
#' estimado por sounding (delta), na grade de 1 grau do MIP.
#'
#' Três séries por célula x ano, todas em g C m-2 d-1:
#'   beta : tendência intra-anual do XCO2 (beta_line, ppm/mês, grade de 0,5 grau do repositório)
#'          reagrupada em células de 1 grau (média das subcélulas) e convertida em fluxo;
#'   est  : fluxo estimado por sounding (df_comp), média anual dos meses com estimado e net_flux;
#'   ref  : net_flux do MIP, média anual dos mesmos meses.
#' O beta é uma taxa média do ano (positivo = XCO2 subindo = fonte); por isso é comparado com a
#' média anual dos fluxos, não com valores mensais.
#'
#' @param beta_df   data.frame do seu `dfall` (ou output/beta_significant.xlsx lido): lon, lat,
#'                  year, beta_line (ppm por mês, pois linear_reg usa x = mês) e, opcional,
#'                  beta_error e a coluna de classe `xco2`.
#' @param df_comp   Saída de agg_xco2_flux(): lon, lat, year, month, `est` e `ref`.
#' @param biomes    sf dos biomas (geobr::read_biomes). Habilita as métricas e gráficos por bioma.
#' @param est,ref   Colunas do fluxo estimado e do net_flux (g C m-2 d-1; ref_factor ajusta o ref).
#' @param unit      "physical": beta (ppm/mês) -> mg CO2 m-2 mês-1 (1e4 * beta * 44/24,45) -> g CO2
#'                  m-2 mês-1 (/1000) -> g C m-2 d-1 (x 12,011/44,009 / 30).
#'                  "repo": usa o beta_fco2 do repositório como está, que multiplica por 30
#'                  (beta_molm * 30/1000) e depois o lê como g CO2 m-2 mês-1. Os dois diferem por
#'                  um fator 30 (veja a mensagem de diagnóstico). Não afeta Spearman.
#' @param classes   Mantém só estas classes do beta (ex.: c("Source","Sink")). NULL = todas.
#' @param class_col Coluna de classe do beta_df (padrão "xco2").
#' @param min_months  Meses pareados mínimos para a média anual do est/ref (padrão 6).
#' @param min_sub   Mínimo de subcélulas de 0,5 grau por célula de 1 grau (padrão 1).
#' @param res       Resolução da grade final em graus (padrão 1).
#' @param maps      TRUE = gera os mapas anuais (fluxo, bias, r, RMSE, slope) com map_flux_annual() para
#'                  beta x net_flux e beta x estimado (precisa de map_flux_annual.R carregado).
#' @param map_args  Lista de argumentos extras para map_flux_annual(), ex.:
#'                  list(south_america = south_america, map_theme = map_theme, ncol = 3).
#' @param min_years Anos mínimos por célula para a regressão célula a célula dos mapas.
#' @return lista: pairs (lon, lat, year, beta, est, ref, ...), metrics, plots.
compare_beta <- function(beta_df,
                         df_comp,
                         biomes       = NULL,
                         biome_col    = "name_biome",
                         est          = "fco2_mean",
                         ref          = "net_flux",
                         ref_factor   = 1,
                         unit         = c("physical", "repo"),
                         classes      = NULL,
                         class_col    = "xco2",
                         min_months   = 6,
                         min_sub      = 1,
                         res          = 1,
                         min_days     = 1,
                         plot         = TRUE,
                         maps         = FALSE,
                         map_args     = list(),
                         min_years    = 5) {

  unit <- match.arg(unit)
  if (!all(c("lon", "lat", "year", "beta_line") %in% names(beta_df)))
    stop("beta_df precisa de lon, lat, year e beta_line.")
  if (!all(c("lon", "lat", "year", "month", est, ref) %in% names(df_comp)))
    stop("df_comp precisa de lon, lat, year, month, ", est, " e ", ref, ".")

  # ---- beta -> fluxo (g C m-2 d-1) e reagrupamento 0,5 -> 1 grau --------------
  b <- as.data.frame(beta_df)
  if (!is.null(classes)) b <- b[b[[class_col]] %in% classes, ]
  g_co2_month_phys <- b$beta_line * 1e4 * 44 / 24.45 / 1000          # g CO2 m-2 mes-1 (beta por mes)
  g_co2_month_repo <- b$beta_line * 1e4 * 44 / 24.45 * 30 / 1000     # beta_fco2 do repositorio
  to_gC_d <- function(x) x * 12.011 / 44.009 / 30
  message(sprintf("Mediana do beta em g C m-2 d-1: fisico = %.2e | repo (x30) = %.2e",
                  stats::median(to_gC_d(g_co2_month_phys), na.rm = TRUE),
                  stats::median(to_gC_d(g_co2_month_repo), na.rm = TRUE)))
  b$beta <- to_gC_d(if (unit == "physical") g_co2_month_phys else g_co2_month_repo)
  b$lon1 <- round(floor(b$lon / res) * res + res / 2, 4)
  b$lat1 <- round(floor(b$lat / res) * res + res / 2, 4)
  bg <- b |>
    dplyr::group_by(lon = lon1, lat = lat1, year) |>
    dplyr::summarise(beta = mean(beta, na.rm = TRUE), n_sub = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n_sub >= min_sub)

  # ---- est e ref: média anual dos meses pareados ------------------------------
  d <- as.data.frame(df_comp)
  d$est <- d[[est]]; d$ref <- d[[ref]] * ref_factor
  d$lon <- round(d$lon, 4); d$lat <- round(d$lat, 4)
  if ("n_days" %in% names(d)) d <- d[is.na(d$n_days) | d$n_days >= min_days, ]
  d <- d[is.finite(d$est) & is.finite(d$ref), ]
  ann <- d |>
    dplyr::group_by(lon, lat, year) |>
    dplyr::summarise(est = mean(est), ref = mean(ref), n_months = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n_months >= min_months)

  pairs <- dplyr::inner_join(bg, ann, by = c("lon", "lat", "year"))
  if (nrow(pairs) == 0) stop("Nenhuma célula x ano em comum entre beta e df_comp (confira as grades).")
  message(sprintf("%d células x ano em comum (%d células, anos %s).", nrow(pairs),
                  nrow(unique(pairs[c("lon", "lat")])), paste(range(pairs$year), collapse = "-")))

  # ---- bioma ----------------------------------------------------------------------
  if (!is.null(biomes)) {
    s2_old <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(s2_old), add = TRUE)
    bb  <- sf::st_make_valid(sf::st_transform(biomes, 4326))
    xy  <- unique(pairs[, c("lon", "lat")])
    pts <- sf::st_as_sf(xy, coords = c("lon", "lat"), crs = 4326)
    hit <- sf::st_intersects(pts, bb)
    idx <- vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
    xy$biome <- as.character(sf::st_drop_geometry(bb)[[biome_col]])[idx]
    pairs <- merge(pairs, xy, by = c("lon", "lat"), all.x = TRUE, sort = FALSE)
    pairs <- tibble::as_tibble(pairs[!is.na(pairs$biome), ])
  }

  # ---- métricas -----------------------------------------------------------------------
  cmp <- list(c("beta", "ref"), c("beta", "est"), c("est", "ref"))
  nm  <- c(beta = "Beta flux", est = "Delta-based estimate", ref = "MIP net flux")
  one <- function(x, y, group, a, b_) {
    ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
    if (n < 5 || stats::sd(x) == 0 || stats::sd(y) == 0)
      return(tibble::tibble(group = group, comparison = paste(nm[a], "vs", nm[b_]), n = n,
                            rho = NA_real_, rho_p = NA_real_, r = NA_real_, rmse = NA_real_, bias = NA_real_))
    ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
    tibble::tibble(group = group, comparison = paste(nm[a], "vs", nm[b_]), n = n,
                   rho = unname(ct$estimate), rho_p = ct$p.value, r = stats::cor(x, y),
                   rmse = sqrt(mean((x - y)^2)), bias = mean(x - y))   # primeiro - segundo
  }
  groups <- c("Brazil", if (!is.null(biomes)) sort(unique(pairs$biome)))
  metrics <- dplyr::bind_rows(lapply(groups, function(gp) {
    pp <- if (gp == "Brazil") pairs else pairs[pairs$biome == gp, ]
    dplyr::bind_rows(lapply(cmp, function(cc) one(pp[[cc[1]]], pp[[cc[2]]], gp, cc[1], cc[2])))
  }))

  plots <- NULL
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- .make_beta_plots(pairs, metrics, cmp, nm)
  }
  maps_out <- NULL
  if (maps) {
    if (!exists("map_flux_annual", mode = "function"))
      stop("Carregue map_flux_annual.R (source) antes de usar maps = TRUE.")
    pdf <- as.data.frame(pairs[, c("lon", "lat", "year", "beta", "est", "ref")]); pdf$month <- 1L
    run <- function(e, r)
      do.call(map_flux_annual, c(list(df = pdf, biomes = biomes, biome_col = biome_col, est = e, ref = r,
                                      min_months = 1, min_years = min_years), map_args))
    maps_out <- list(beta_vs_ref = run("beta", "ref"), beta_vs_est = run("beta", "est"))
  }
  list(pairs = pairs, metrics = metrics, plots = plots, maps = maps_out)
}

# ---------------------------------------------------------------------------
# Gráficos (inglês, sem título)
# ---------------------------------------------------------------------------
.make_beta_plots <- function(pairs, metrics, cmp, nm) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  unit <- quote(g ~ C ~ m^{-2} ~ d^{-1})
  th <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0))
  fp <- function(p) ifelse(!is.finite(p), "p = NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
  lab <- function(m) paste(paste0("rho = ", ifelse(is.finite(m$rho), sprintf("%.2f", m$rho), "NA")),
                           fp(m$rho_p), sprintf("RMSE = %.3f", m$rmse), sprintf("Bias = %.3f", m$bias), sep = "\n")
  sc <- function(d, xv, yv, facet, labs_df, xl, yl) {
    d$X <- d[[xv]]; d$Y <- d[[yv]]
    gg(d, aes(X, Y)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.4) +
      ggplot2::geom_point(colour = "#2a78d6", alpha = 0.35, size = 1.1) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey15", linewidth = 0.6) +
      ggplot2::geom_text(data = labs_df, aes(x = -Inf, y = Inf, label = lab), hjust = -0.08, vjust = 1.25,
                         size = 3.2, lineheight = 1.05, colour = "grey15", inherit.aes = FALSE) +
      ggplot2::facet_wrap(facet, scales = "free") +
      ggplot2::labs(x = bquote(.(xl) ~ "(" * .(unit) * ")"), y = bquote(.(yl) ~ "(" * .(unit) * ")")) + th
  }
  out <- list()
  # Brasil: os três pares lado a lado
  long <- dplyr::bind_rows(lapply(cmp, function(cc) data.frame(
    pair = paste(nm[cc[1]], "vs", nm[cc[2]]), x = pairs[[cc[2]]], y = pairs[[cc[1]]])))
  lb <- metrics[metrics$group == "Brazil", ]; lb$pair <- lb$comparison; lb$lab <- lab(lb)
  out$scatter_brazil <- sc(long, "x", "y", ggplot2::vars(pair), lb, "Reference series", "First series")
  # por bioma: um gráfico por par
  if ("biome" %in% names(pairs)) {
    for (cc in cmp) {
      mm <- metrics[metrics$comparison == paste(nm[cc[1]], "vs", nm[cc[2]]) & metrics$group != "Brazil", ]
      if (nrow(mm) == 0) next
      mm$biome <- mm$group; mm$lab <- lab(mm)
      out[[paste0("scatter_", cc[1], "_", cc[2])]] <-
        sc(pairs, cc[2], cc[1], ggplot2::vars(biome), mm, nm[[cc[2]]], nm[[cc[1]]])
    }
  }
  out
}

#' Balanço anual e total (Tg C) dos três métodos: beta, estimativa por delta e net_flux
#'
#' Usa as células x anos pareadas de compare_beta() (mesmas células para os três métodos).
#' balanço = fluxo médio anual (g C m-2 d-1) x área da célula x dias do ano.
#'
#' @param x saída de compare_beta() (lista com $pairs) ou o próprio data frame de pares
#'          (lon, lat, year, beta, est, ref e, opcionalmente, biome).
#' @param res resolução da grade em graus (1).
#' @param km_per_degree km por grau (110 -> célula de 1 grau = 110 x 110 km).
#' @param plot se TRUE, devolve também os gráficos (3 barras: beta, estimado, net_flux).
#' @return lista: year (Tg C por grupo x ano), total (Tg C por grupo), plots (balance_year, balance_total).
balance_beta <- function(x, res = 1, km_per_degree = 110, plot = TRUE) {
  pw <- tibble::as_tibble(if (is.data.frame(x)) x else x$pairs)
  if (!all(c("year", "beta", "est", "ref") %in% names(pw))) stop("x precisa de year, beta, est e ref.")
  yr_days <- function(y) ifelse((y %% 4 == 0 & y %% 100 != 0) | y %% 400 == 0, 366, 365)
  pw$w <- (res * km_per_degree * 1e3)^2 * yr_days(pw$year) / 1e12     # g C m-2 d-1 -> Tg C ano-1
  groups <- c("Brazil", if ("biome" %in% names(pw)) sort(unique(pw$biome)))
  one <- function(d, gp, by_year) {
    d$beta_TgC <- d$beta * d$w; d$est_TgC <- d$est * d$w; d$ref_TgC <- d$ref * d$w
    o <- if (by_year) {
      dplyr::summarise(dplyr::group_by(d, year), n_cells = dplyr::n(), beta_TgC = sum(beta_TgC),
                       est_TgC = sum(est_TgC), ref_TgC = sum(ref_TgC), .groups = "drop")
    } else {
      dplyr::summarise(d, n_cell_years = dplyr::n(), beta_TgC = sum(beta_TgC),
                       est_TgC = sum(est_TgC), ref_TgC = sum(ref_TgC))
    }
    dplyr::bind_cols(tibble::tibble(group = gp), o)
  }
  sel <- function(gp) if (gp == "Brazil") pw else pw[pw$biome == gp, ]
  by <- dplyr::bind_rows(lapply(groups, function(gp) one(sel(gp), gp, TRUE)))
  bt <- dplyr::bind_rows(lapply(groups, function(gp) one(sel(gp), gp, FALSE)))
  plots <- NULL
  if (plot && requireNamespace("ggplot2", quietly = TRUE))
    plots <- .make_beta_balance_plots(by, bt,
                                      c(beta = "Beta flux", est = "Delta-based estimate", ref = "MIP net flux"))
  list(year = by, total = bt, plots = plots)
}

# Balanço anual e total (Tg C): três barras por grupo/ano (beta, estimado, net_flux)
.make_beta_balance_plots <- function(by, bt, nm) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  cols <- c("Beta flux" = "#2a9d8f", "Delta-based estimate" = "#2a78d6", "MIP net flux" = "#e07b39")
  th <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "top")
  lv <- unique(by$group)
  long <- function(d, id) {
    o <- do.call(rbind, lapply(c("beta", "est", "ref"), function(k)
      data.frame(d[id], series = unname(nm[k]), TgC = d[[paste0(k, "_TgC")]])))
    o$series <- factor(o$series, levels = unname(nm[c("beta", "est", "ref")]))
    o$group <- factor(o$group, levels = lv); o
  }
  dodge <- ggplot2::position_dodge(width = 0.8)
  ly <- long(as.data.frame(by), c("group", "year"))
  lt <- long(as.data.frame(bt), "group")
  list(
    balance_year = gg(ly, aes(factor(year), TgC, fill = series)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
      ggplot2::geom_col(position = dodge, width = 0.75) +
      ggplot2::facet_wrap(~ group, scales = "free_y", ncol = 2) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(x = NULL, y = "Annual balance (Tg C)", fill = NULL) + th +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
    balance_total = gg(lt, aes(group, TgC, fill = series)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
      ggplot2::geom_col(position = dodge, width = 0.75) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(x = NULL, y = "Total balance (Tg C, sum over cells and years)", fill = NULL) + th +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  )
}

#' Salva os gráficos de balance_beta() em PNG
save_beta_balance <- function(bal, dir = ".", dpi = 300) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(dir, "beta_balance_year.png"), bal$plots$balance_year, width = 11, height = 9, dpi = dpi, bg = "white")
  ggplot2::ggsave(file.path(dir, "beta_balance_total.png"), bal$plots$balance_total, width = 8, height = 4.5, dpi = dpi, bg = "white")
}

#' Salva os gráficos de compare_beta() em PNG
save_beta_plots <- function(res, dir = ".", dpi = 300) {
  stopifnot(!is.null(res$plots)); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(res$plots)) {
    sz <- if (nm == "scatter_brazil") c(11, 4) else c(11, 8)
    ggplot2::ggsave(file.path(dir, paste0("beta_", nm, ".png")), res$plots[[nm]],
                    width = sz[1], height = sz[2], dpi = dpi, bg = "white")
  }
}

#' Salva os mapas de compare_beta(maps = TRUE): uma subpasta por comparação
save_beta_maps <- function(res, dir = ".", dpi = 300, panel = 6) {
  stopifnot(!is.null(res$maps))
  for (nm in names(res$maps)) save_flux_maps(res$maps[[nm]], file.path(dir, nm), dpi = dpi, panel = panel)
}

# ---------------------------------------------------------------------------
# Exemplo de uso
# ---------------------------------------------------------------------------
# dfall   <- readxl::read_xlsx("output/beta_significant.xlsx")     # ou o dfall do seu loop
# df_comp <- readRDS("teste_comp.rds")                             # agg_xco2_flux() com net_flux
# biomas  <- geobr::read_biomes(year = 2019)
# res <- compare_beta(dfall, df_comp, biomes = biomas, unit = "physical")
# res$metrics        # rho, p, r, RMSE, bias: Brasil e biomas, para os 3 pares
# res$plots$scatter_brazil; res$plots$scatter_beta_ref
# save_beta_plots(res, "figs_beta/")
#
# # balanço anual e total (Tg C), 3 barras: beta, estimado, net_flux
# bal <- balance_beta(res)        # bal$year, bal$total, bal$plots$balance_year / balance_total
# save_beta_balance(bal, "figs_beta/")
#
# # mapas anuais (mesmo estilo dos outros), beta x net_flux e beta x estimado:
# source("map_flux_annual.R")
# res <- compare_beta(dfall, df_comp, biomes = biomas, maps = TRUE,
#                     map_args = list(south_america = south_america, map_theme = map_theme))
# save_beta_maps(res, "figs_beta_maps/")   # subpastas beta_vs_ref/ e beta_vs_est/
# res$maps$beta_vs_ref$summary_biome       # métricas da regressão célula a célula por bioma
