#' Compara o fluxo estimado (por sondagem -> célula x mês) com o net_flux do MIP
#'
#' Compara, nas células x mês em que existem as duas informações (pareadas):
#'  - Brasil todo e por bioma: métricas de concordância (viés, RMSE, r, inclinação...);
#'  - balanço total (Tg C) mensal, anual e acumulado: estimado x net_flux, com a
#'    incerteza propagada do estimado.
#'
#' Só células x mês com os DOIS valores entram nas métricas e nos balanços, para que
#' as duas séries cubram exatamente a mesma área/período. Para ver quanto do net_flux
#' total do Brasil essa amostra representa, informe `flux_all`.
#'
#' @param df         Saída de agg_xco2_flux() (uma linha por célula x mês), com
#'                   lon, lat, year, month e as colunas `est`, `ref`.
#' @param flux_all   (opcional) Tabela completa de fluxo (saída do flux_extractor(),
#'                   mesmo experimento/submission do `df`) com lon, lat, year, month e
#'                   `ref`. Serve para calcular o net_flux de TODAS as células e a
#'                   fração da área x tempo coberta pelas sondagens.
#' @param biomes     sf com os polígonos dos biomas. NULL = tenta geobr::read_biomes(2019).
#'                   Se não houver geobr/biomas, só o Brasil é analisado.
#' @param biome_col  Coluna com o nome do bioma em `biomes` (padrão "name_biome").
#' @param est        Coluna do fluxo estimado (padrão "fco2_mean", g C m-2 d-1).
#' @param est_sigma  Coluna da incerteza do estimado na célula x mês (padrão
#'                   "sigma_f_cell"). NULL ou inexistente = sem incerteza.
#' @param ref_sigma  Coluna com a incerteza (1 desvio padrão, mesma unidade de `ref`) do
#'                   fluxo de referência por célula x mês, ex.: o desvio do ensemble do MIP
#'                   (arquivo EnsStd). NULL = a incerteza da referência não entra.
#' @param ref        Coluna do fluxo de referência (padrão "net_flux").
#' @param ref_factor Fator que converte `ref` para g C m-2 d-1 (padrão 1 = já está
#'                   nessa unidade). CONFIRA a unidade do net_flux no README do MIP.
#' @param res        Resolução da grade em graus (1 ou 2 valores). NULL = inferida.
#' @param drop_outside Com `biomes`, exclui as células cujo centro cai fora de todos os biomas
#'                   (fora do Brasil): "Brasil" passa a ser só o que está dentro dos biomas e
#'                   o grupo "Outside biomes" não existe. FALSE = mantém essas células no Brasil.
#' @param km_per_degree Km por grau usado na área da célula (padrão 110: célula de 1 grau =
#'                   110 km x 110 km = 12.100 km2, igual para todas as células, sem ponderação
#'                   por latitude).
#' @param min_days   Mínimo de dias com sondagem na célula x mês (usa `n_days`).
#' @param years      Anos a considerar. NULL = todos.
#' @param sigma_bars Barras de incerteza no gráfico do balanço: "both" (linha grossa = células
#'                   independentes + haste = totalmente correlacionadas), "indep" ou "corr".
#'                   Com erros correlacionados a haste costuma ser muito maior que o balanço;
#'                   "indep" deixa o gráfico legível (mostra o limite inferior).
#' @param plot       Gera os gráficos (requer ggplot2).
#'
#' @return lista com: cells (dados pareados + bioma + área), metrics (por grupo),
#'         balance_monthly, balance_year, balance_total (Tg C), flux_monthly
#'         (fluxo médio mensal das células pareadas, g C m-2 d-1) e plots (em inglês).
compare_flux <- function(df,
                         flux_all   = NULL,
                         biomes     = NULL,
                         biome_col  = "name_biome",
                         est        = "fco2_mean",
                         est_sigma  = "sigma_f_cell",
                         ref_sigma  = NULL,
                         ref        = "net_flux",
                         ref_factor = 1,
                         res        = NULL,
                         km_per_degree = 110,
                         drop_outside = TRUE,
                         min_days   = 1,
                         years      = NULL,
                         sigma_bars = c("both", "indep", "corr"),
                         plot       = TRUE) {

  sigma_bars <- match.arg(sigma_bars)

  # ---- checagens e preparo --------------------------------------------------
  need <- c("lon", "lat", "year", "month", est, ref)
  if (!all(need %in% names(df)))
    stop("df precisa das colunas: ", paste(setdiff(need, names(df)), collapse = ", "))
  if (!is.null(flux_all) && !all(c("lon", "lat", "year", "month", ref) %in% names(flux_all)))
    stop("flux_all precisa de lon, lat, year, month e '", ref, "'.")

  if (ref_factor == 1)
    message("Assumindo que '", ref, "' já está em g C m-2 d-1 (ref_factor = 1). ",
            "Confira a unidade no README do MIP.")

  df <- as.data.frame(df)
  df$est <- df[[est]]
  df$ref <- df[[ref]] * ref_factor
  df$sig <- if (!is.null(est_sigma) && est_sigma %in% names(df)) df[[est_sigma]] else NA_real_
  if (!is.null(ref_sigma) && !ref_sigma %in% names(df))
    stop("ref_sigma '", ref_sigma, "' não existe em df.")
  df$rsig <- if (!is.null(ref_sigma)) df[[ref_sigma]] * ref_factor else NA_real_
  df$lon <- round(df$lon, 4); df$lat <- round(df$lat, 4)

  if ("n_days" %in% names(df)) df <- df[is.na(df$n_days) | df$n_days >= min_days, ]
  if (!is.null(years)) df <- df[df$year %in% years, ]
  paired <- df[is.finite(df$est) & is.finite(df$ref), ]
  if (nrow(paired) == 0) stop("Nenhuma célula x mês com estimado e net_flux ao mesmo tempo.")

  if (!is.null(flux_all)) {
    flux_all <- as.data.frame(flux_all)
    flux_all$ref <- flux_all[[ref]] * ref_factor
    flux_all$lon <- round(flux_all$lon, 4); flux_all$lat <- round(flux_all$lat, 4)
    if (!is.null(years)) flux_all <- flux_all[flux_all$year %in% years, ]
    flux_all <- flux_all[is.finite(flux_all$ref), ]
  }

  # ---- grade, área da célula e dias do mês ---------------------------------
  if (is.null(res)) {
    step <- function(v) { v <- sort(unique(round(v, 4))); if (length(v) < 2) NA_real_ else min(diff(v)) }
    xy <- rbind(paired[, c("lon", "lat")], if (!is.null(flux_all)) flux_all[, c("lon", "lat")])
    res <- c(step(xy$lon), step(xy$lat))
    if (anyNA(res)) stop("Não consegui inferir a resolução da grade; informe 'res'.")
  }
  if (length(res) == 1) res <- rep(res, 2)
  message(sprintf("Grade: %g x %g graus; area da celula = %.0f km2 (%g km/grau). Confira se esta correto.",
                  res[1], res[2], (res[1] * km_per_degree) * (res[2] * km_per_degree), km_per_degree))

  rad <- pi / 180
  # área da célula constante: km_per_degree x km_per_degree por grau (padrao 110 km x 110 km
  # para 1 grau); sem ponderacao por latitude
  area_m2 <- function(lat) rep((res[1] * km_per_degree * 1e3) * (res[2] * km_per_degree * 1e3), length(lat))
  days_in_month <- function(y, m) {
    d1  <- as.Date(sprintf("%04d-%02d-01", y, m))
    nxt <- as.Date(ifelse(m == 12, sprintf("%04d-01-01", y + 1L),
                          sprintf("%04d-%02d-01", y, m + 1L)))
    as.integer(nxt - d1)
  }
  paired$area_m2 <- area_m2(paired$lat)
  paired$ndays   <- days_in_month(paired$year, paired$month)
  if (!is.null(flux_all)) {
    flux_all$area_m2 <- area_m2(flux_all$lat)
    flux_all$ndays   <- days_in_month(flux_all$year, flux_all$month)
  }

  # ---- biomas ---------------------------------------------------------------
  if (is.null(biomes) && requireNamespace("geobr", quietly = TRUE)) {
    biomes <- tryCatch(geobr::read_biomes(year = 2019, simplified = TRUE, showProgress = FALSE),
                       error = function(e) { message("Não consegui baixar os biomas: ", conditionMessage(e)); NULL })
  }
  has_biome <- !is.null(biomes)
  if (!has_biome) message("Sem polígonos de biomas: só a análise do Brasil será feita.")

  assign_biome <- function(xy) {
    s2_old <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(s2_old), add = TRUE)
    b   <- sf::st_make_valid(sf::st_transform(biomes, 4326))
    pts <- sf::st_as_sf(xy, coords = c("lon", "lat"), crs = 4326)
    hit <- sf::st_intersects(pts, b)
    idx <- vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
    out <- as.character(sf::st_drop_geometry(b)[[biome_col]])[idx]
    out[is.na(out)] <- "Fora dos biomas"
    out
  }
  if (has_biome) {
    cells <- unique(rbind(paired[, c("lon", "lat")], if (!is.null(flux_all)) flux_all[, c("lon", "lat")]))
    cells$biome <- assign_biome(cells)
    paired <- merge(paired, cells, by = c("lon", "lat"), all.x = TRUE, sort = FALSE)
    if (!is.null(flux_all)) flux_all <- merge(flux_all, cells, by = c("lon", "lat"), all.x = TRUE, sort = FALSE)
    if (drop_outside) {
      # Brasil = apenas celulas dentro de algum bioma; nada fora dos biomas entra em tabelas/graficos
      n_out <- length(unique(paste(paired$lon, paired$lat)[paired$biome == "Fora dos biomas"]))
      if (n_out > 0) message(sprintf("%d celula(s) fora dos biomas (centro fora dos poligonos) excluidas.", n_out))
      paired <- paired[paired$biome != "Fora dos biomas", ]
      if (!is.null(flux_all)) flux_all <- flux_all[flux_all$biome != "Fora dos biomas", ]
      if (nrow(paired) == 0) stop("Nenhuma celula pareada dentro dos biomas.")
    }
  }

  # ---- grupos: Brasil (todas as células) + cada bioma -----------------------
  add_groups <- function(d) {
    out <- cbind(d, group = "Brasil", stringsAsFactors = FALSE)
    if (has_biome) out <- rbind(out, cbind(d, group = d$biome, stringsAsFactors = FALSE))
    out
  }
  P <- add_groups(paired)
  A <- if (!is.null(flux_all)) add_groups(flux_all) else NULL

  # ---- métricas de concordância --------------------------------------------
  metr <- function(d) {
    x <- d$ref; y <- d$est; e <- y - x
    ok <- nrow(d) >= 3 && stats::sd(x) > 0 && stats::sd(y) > 0
    tibble::tibble(
      n_cellmonths = nrow(d),
      n_cells      = nrow(unique(d[, c("lon", "lat")])),
      mean_ref     = mean(x),
      mean_est     = mean(y),
      bias         = mean(e),                     # estimado - net_flux
      mae          = mean(abs(e)),
      rmse         = sqrt(mean(e^2)),
      r            = if (ok) stats::cor(x, y) else NA_real_,
      rho          = if (ok) stats::cor(x, y, method = "spearman") else NA_real_,   # Spearman
      rho_p        = if (ok) suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE)$p.value) else NA_real_,
      slope        = if (ok) unname(stats::coef(stats::lm(y ~ x))[2]) else NA_real_,
      frac_within_2sigma = if (all(is.na(d$sig))) NA_real_ else mean(abs(e) <= 2 * d$sig, na.rm = TRUE)
    )
  }
  metrics <- P |>
    dplyr::group_by(group) |>
    dplyr::group_modify(~ metr(.x)) |>
    dplyr::ungroup()

  # ---- balanços (Tg C) ------------------------------------------------------
  TG <- 1e12   # g por Tg
  sum_or_na <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
  balance_monthly <- P |>
    dplyr::mutate(w = area_m2 * ndays / TG) |>
    dplyr::group_by(group, year, month) |>
    dplyr::summarise(
      n_cells  = dplyr::n(),
      est_TgC  = sum(est * w),
      ref_TgC  = sum(ref * w),
      est_sigma_indep_TgC = if (all(is.na(sig))) NA_real_ else sqrt(sum((sig * w)^2, na.rm = TRUE)),
      est_sigma_corr_TgC  = sum_or_na(sig * w),
      ref_sigma_indep_TgC = if (all(is.na(rsig))) NA_real_ else sqrt(sum((rsig * w)^2, na.rm = TRUE)),
      ref_sigma_corr_TgC  = sum_or_na(rsig * w),
      area_ndays_paired   = sum(w),
      .groups = "drop"
    ) |>
    dplyr::mutate(diff_TgC = est_TgC - ref_TgC)

  if (!is.null(A)) {
    ref_all <- A |>
      dplyr::mutate(w = area_m2 * ndays / TG) |>
      dplyr::group_by(group, year, month) |>
      dplyr::summarise(ref_all_TgC = sum(ref * w), area_ndays_all = sum(w), .groups = "drop")
    balance_monthly <- balance_monthly |>
      dplyr::left_join(ref_all, by = c("group", "year", "month")) |>
      dplyr::mutate(coverage = area_ndays_paired / area_ndays_all)   # fração área x tempo amostrada
  }

  # incerteza de (estimado - referencia): estimado e referencia independentes entre si (quadratura);
  # se a incerteza da referencia nao foi informada (NA), so a do estimado entra
  add_diff_sigma <- function(b) {
    r0 <- function(x) ifelse(is.na(x), 0, x)
    b$diff_sigma_indep_TgC <- sqrt(b$est_sigma_indep_TgC^2 + r0(b$ref_sigma_indep_TgC)^2)
    b$diff_sigma_corr_TgC  <- sqrt(b$est_sigma_corr_TgC^2  + r0(b$ref_sigma_corr_TgC)^2)
    b
  }
  balance_monthly <- add_diff_sigma(balance_monthly)

  agg_bal <- function(b, ...) {
    b |>
      dplyr::group_by(group, ...) |>
      dplyr::summarise(
        n_months = dplyr::n(),
        est_TgC  = sum(est_TgC),
        ref_TgC  = sum(ref_TgC),
        diff_TgC = sum(diff_TgC),
        # meses independentes: soma em quadratura; totalmente correlacionados: soma direta
        est_sigma_indep_TgC = if (all(is.na(est_sigma_indep_TgC))) NA_real_ else sqrt(sum(est_sigma_indep_TgC^2, na.rm = TRUE)),
        est_sigma_corr_TgC  = sum_or_na(est_sigma_corr_TgC),
        ref_sigma_indep_TgC = if (all(is.na(ref_sigma_indep_TgC))) NA_real_ else sqrt(sum(ref_sigma_indep_TgC^2, na.rm = TRUE)),
        ref_sigma_corr_TgC  = sum_or_na(ref_sigma_corr_TgC),
        dplyr::across(dplyr::any_of(c("ref_all_TgC")), sum),
        .groups = "drop"
      ) |>
      add_diff_sigma()
  }
  balance_year  <- agg_bal(balance_monthly, year)
  balance_total <- agg_bal(balance_monthly)
  if ("area_ndays_all" %in% names(balance_monthly)) {
    cov_tot <- balance_monthly |>
      dplyr::group_by(group) |>
      dplyr::summarise(coverage = sum(area_ndays_paired) / sum(area_ndays_all), .groups = "drop")
    balance_total <- dplyr::left_join(balance_total, cov_tot, by = "group")
  }

  # ---- fluxo mensal médio (g C m-2 d-1) das células pareadas (área igual) ----
  # média espacial das células pareadas; incerteza do estimado na média:
  #  indep = células independentes (limite inferior); corr = totalmente correlacionadas (superior)
  flux_monthly <- P |>
    dplyr::group_by(group, year, month) |>
    dplyr::summarise(
      n_cells   = dplyr::n(),
      est_flux  = stats::weighted.mean(est, area_m2),
      ref_flux  = stats::weighted.mean(ref, area_m2),
      est_sigma_indep = if (all(is.na(sig))) NA_real_ else sqrt(sum((sig * area_m2)^2, na.rm = TRUE)) / sum(area_m2),
      est_sigma_corr  = if (all(is.na(sig))) NA_real_ else sum(sig * area_m2, na.rm = TRUE) / sum(area_m2),
      ref_sigma_indep = if (all(is.na(rsig))) NA_real_ else sqrt(sum((rsig * area_m2)^2, na.rm = TRUE)) / sum(area_m2),
      ref_sigma_corr  = if (all(is.na(rsig))) NA_real_ else sum(rsig * area_m2, na.rm = TRUE) / sum(area_m2),
      .groups = "drop"
    )

  # ---- gráficos (rótulos em inglês) ----------------------------------------
  plots <- NULL
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- make_flux_plots(P, balance_total, flux_monthly, metrics, ref = ref, sigma_bars = sigma_bars,
                             balance_year = balance_year)
  }

  list(cells = tibble::as_tibble(paired), metrics = metrics,
       balance_monthly = balance_monthly, balance_year = balance_year,
       balance_total = balance_total, flux_monthly = flux_monthly, plots = plots)
}

# ---------------------------------------------------------------------------
# Gráficos (figuras em inglês)
# ---------------------------------------------------------------------------
# Cores categóricas validadas (azul/laranja: ΔE CVD 24.7, visão normal 33.6, contraste >= 3:1)
.flux_cols <- c(est = "#2a78d6", ref = "#eb6834")

# Nomes de exibição em inglês (grupos); o que não estiver na lista fica como está
.group_en <- function(x) {
  map <- c("Brasil" = "Brazil",
           "Amazônia" = "Amazon", "Amazonia" = "Amazon",
           "Cerrado" = "Cerrado", "Caatinga" = "Caatinga",
           "Mata Atlântica" = "Atlantic Forest", "Mata Atlantica" = "Atlantic Forest",
           "Pampa" = "Pampa", "Pantanal" = "Pantanal",
           "Sistema Costeiro" = "Coastal System",
           "Fora dos biomas" = "Outside biomes")
  out <- unname(map[x]); out[is.na(out)] <- x[is.na(out)]; out
}
.group_levels <- function(x) {                       # Brazil primeiro, depois ordem alfabética
  u <- unique(.group_en(x)); c(intersect("Brazil", u), sort(setdiff(u, "Brazil")))
}

make_flux_plots <- function(P, balance_total, flux_monthly, metrics, ref = "net_flux",
                            sigma_bars = "both", balance_year = NULL) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  ref_lab <- "MIP Net Flux"
  lev     <- .group_levels(c(P$group, balance_total$group))
  fac     <- function(x) factor(.group_en(x), levels = lev)
  unit    <- quote(g ~ C ~ m^{-2} ~ d^{-1})
  base_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "top")
  col_scale <- ggplot2::scale_colour_manual(values = c(Estimated = unname(.flux_cols["est"]),
                                                       setNames(unname(.flux_cols["ref"]), ref_lab)))

  # 1) dispersão estimado x referência, por grupo (Brasil + biomas)
  Pp <- P; Pp$grp <- fac(Pp$group)
  # correlacao de Pearson (r), RMSE e vies (estimado - referencia)
  stats_lab <- data.frame(grp = fac(metrics$group))
  fmt_p <- function(p) ifelse(!is.finite(p), "p = NA", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
  stats_lab$lab <- paste(
    paste0("rho = ", ifelse(is.finite(metrics$rho), sprintf("%.2f", metrics$rho), "NA")),   # Spearman
    fmt_p(metrics$rho_p),
    sprintf("RMSE = %.2f", metrics$rmse),
    sprintf("Bias = %.2f", metrics$bias), sep = "\n")
  p_scatter <- gg(Pp, aes(ref, est)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_point(colour = .flux_cols["est"], alpha = 0.35, size = 1.1) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey15", linewidth = 0.6) +
    ggplot2::geom_text(data = stats_lab, aes(x = -Inf, y = Inf, label = lab),
                       hjust = -0.08, vjust = 1.25, size = 3.2, lineheight = 1.05,
                       colour = "grey15", inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ grp, scales = "free") +
    ggplot2::labs(x = bquote(.(ref_lab) ~ "(" * .(unit) * ")"),
                  y = bquote("Estimated flux (" * .(unit) * ")")) +
    base_theme

  # 2) balanço total (Tg C) por grupo
  bt <- balance_total; bt$grp <- fac(bt$group)
  long_tot <- rbind(
    data.frame(grp = bt$grp, series = "Estimated", TgC = bt$est_TgC,
               s_in = bt$est_sigma_indep_TgC, s_out = bt$est_sigma_corr_TgC),
    data.frame(grp = bt$grp, series = ref_lab, TgC = bt$ref_TgC,
               s_in = bt$ref_sigma_indep_TgC, s_out = bt$ref_sigma_corr_TgC))
  long_tot$series <- factor(long_tot$series, levels = c("Estimated", ref_lab))
  dodge <- ggplot2::position_dodge(width = 0.75)
  ref_has_sig <- any(is.finite(long_tot$s_out[long_tot$series == ref_lab]))
  p_balance <- gg(long_tot, aes(grp, TgC, fill = series)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_col(position = dodge, width = 0.7) +
    { if (sigma_bars %in% c("both", "corr"))
      ggplot2::geom_errorbar(aes(ymin = TgC - s_out, ymax = TgC + s_out), position = dodge,
                             width = 0.25, colour = "grey15", linewidth = 0.4, na.rm = TRUE) } +
    { if (sigma_bars %in% c("both", "indep"))
      ggplot2::geom_errorbar(aes(ymin = TgC - s_in, ymax = TgC + s_in), position = dodge,
                             width = if (sigma_bars == "indep") 0.25 else 0,
                             colour = "grey15", linewidth = if (sigma_bars == "indep") 0.4 else 1.1,
                             na.rm = TRUE) } +
    ggplot2::scale_fill_manual(values = c(Estimated = unname(.flux_cols["est"]),
                                          setNames(unname(.flux_cols["ref"]), ref_lab))) +
    ggplot2::labs(x = NULL, y = "Tg C (sum over paired cells and months)",
                  fill = NULL) +
    base_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  # 2b) balanço anual (Tg C) por grupo: barras estimado x referência por ano
  p_balance_year <- NULL
  if (!is.null(balance_year)) {
    by <- balance_year; by$grp <- fac(by$group)
    long_y <- rbind(
      data.frame(grp = by$grp, year = by$year, series = "Estimated", TgC = by$est_TgC,
                 s_in = by$est_sigma_indep_TgC, s_out = by$est_sigma_corr_TgC),
      data.frame(grp = by$grp, year = by$year, series = ref_lab, TgC = by$ref_TgC,
                 s_in = by$ref_sigma_indep_TgC, s_out = by$ref_sigma_corr_TgC))
    long_y$series <- factor(long_y$series, levels = c("Estimated", ref_lab))
    yrs <- sort(unique(long_y$year))
    p_balance_year <- gg(long_y, aes(factor(year), TgC, fill = series)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
      ggplot2::geom_col(position = dodge, width = 0.7) +
      { if (sigma_bars %in% c("both", "corr"))
        ggplot2::geom_errorbar(aes(ymin = TgC - s_out, ymax = TgC + s_out), position = dodge,
                               width = 0.25, colour = "grey15", linewidth = 0.35, na.rm = TRUE) } +
      { if (sigma_bars %in% c("both", "indep"))
        ggplot2::geom_errorbar(aes(ymin = TgC - s_in, ymax = TgC + s_in), position = dodge,
                               width = if (sigma_bars == "indep") 0.25 else 0,
                               colour = "grey15", linewidth = if (sigma_bars == "indep") 0.35 else 0.9,
                               na.rm = TRUE) } +
      ggplot2::facet_wrap(~ grp, scales = "free_y", ncol = 2) +
      ggplot2::scale_fill_manual(values = c(Estimated = unname(.flux_cols["est"]),
                                            setNames(unname(.flux_cols["ref"]), ref_lab))) +
      ggplot2::labs(x = NULL, y = "Annual balance (Tg C)", fill = NULL) +
      base_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  # 3) série temporal do FLUXO médio mensal (g C m-2 d-1), por grupo; meses sem dado = lacuna
  fm <- flux_monthly
  fm$date <- as.Date(sprintf("%04d-%02d-15", fm$year, fm$month))
  fm$grp  <- fac(fm$group)
  grid <- expand.grid(grp = factor(lev, levels = lev),
                      date = seq(min(fm$date), max(fm$date), by = "month"))
  grid$date <- as.Date(format(grid$date, "%Y-%m-15"))
  fm <- merge(grid, fm, by = c("grp", "date"), all.x = TRUE)
  fm <- fm[order(fm$grp, fm$date), ]
  long_f <- rbind(
    data.frame(grp = fm$grp, date = fm$date, series = "Estimated", flux = fm$est_flux,
               lo = fm$est_flux - fm$est_sigma_corr, hi = fm$est_flux + fm$est_sigma_corr),
    data.frame(grp = fm$grp, date = fm$date, series = ref_lab, flux = fm$ref_flux,
               lo = fm$ref_flux - fm$ref_sigma_corr, hi = fm$ref_flux + fm$ref_sigma_corr))
  long_f$series <- factor(long_f$series, levels = c("Estimated", ref_lab))
  has_sig <- any(is.finite(long_f$lo))
  p_monthly <- gg(long_f, aes(date, flux, colour = series, group = series)) +
    { if (has_sig) ggplot2::geom_ribbon(aes(ymin = lo, ymax = hi, fill = series),
                                        data = function(d) d[is.finite(d$lo), ],
                                        colour = NA, alpha = 0.15, show.legend = FALSE) } +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.6, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.3, na.rm = TRUE) +
    ggplot2::facet_wrap(~ grp, scales = "free_y", ncol = 2) +
    col_scale +
    ggplot2::scale_fill_manual(values = c(Estimated = unname(.flux_cols["est"]),
                                          setNames(unname(.flux_cols["ref"]), ref_lab))) +
    ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ggplot2::labs(x = NULL, y = bquote("Mean flux (" * .(unit) * ")"), colour = NULL) +
    base_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  out <- list(scatter = p_scatter, balance = p_balance, balance_year = p_balance_year, monthly = p_monthly)
  out[!vapply(out, is.null, logical(1))]
}

#' Salva os gráficos de compare_flux() em PNG
#' @param res  Saída de compare_flux().
#' @param dir  Pasta de saída.
#' @param dpi  Resolução (padrão 300).
save_flux_plots <- function(res, dir = ".", dpi = 300) {
  stopifnot(!is.null(res$plots))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  n <- length(unique(res$flux_monthly$group))
  rows <- ceiling(n / 2)
  spec <- list(scatter = c(11, 3 + 2.6 * ceiling(n / 3)),
               balance = c(8, 5),
               balance_year = c(11, 1.6 + 2.4 * rows),
               monthly = c(11, 1.6 + 2.4 * rows))
  for (nm in names(res$plots)) {
    ggplot2::ggsave(file.path(dir, paste0("flux_", nm, ".png")), res$plots[[nm]],
                    width = spec[[nm]][1], height = spec[[nm]][2], dpi = dpi, bg = "white")
  }
  invisible(file.path(dir, paste0("flux_", names(res$plots), ".png")))
}

# ---------------------------------------------------------------------------
# Exemplo de uso
# ---------------------------------------------------------------------------
# df_comp  <- readRDS("teste_comp.rds")
# df_flux  <- readRDS("teste_flux/EnsMean__LNLGIS.rds")   # fluxo completo (opcional, p/ cobertura)
# biomas   <- geobr::read_biomes(year = 2019)
#
# res <- compare_flux(
#   df_comp,
#   flux_all   = df_flux,
#   biomes     = biomas,
#   min_days   = 3,            # exige >= 3 dias com sondagem na célula x mês
#   years      = 2015:2023,    # período em que existe net_flux
#   ref_factor = 1             # ajuste se o net_flux não estiver em g C m-2 d-1
# )
#
# # Incerteza do MIP (desvio entre modelos, arquivo EnsStd), quando for usar:
# df_std <- flux_extractor(pasta, "teste_flux_std/", geometry = geom,
#                         submissions = "EnsStd", experiments = "LNLGIS") |> invisible()
# df_std <- readRDS("teste_flux_std/EnsStd__LNLGIS.rds") |>
#   dplyr::transmute(lon = round(lon, 4), lat = round(lat, 4), year, month, net_flux_sd = net_flux)
# df_comp2 <- dplyr::left_join(dplyr::mutate(df_comp, lon = round(lon, 4), lat = round(lat, 4)),
#                              df_std, by = c("lon", "lat", "year", "month"))
# res <- compare_flux(df_comp2, ref_sigma = "net_flux_sd", sigma_bars = "indep")
#
# res$metrics          # Brasil + biomas: bias, rmse, r, slope...
# res$balance_total    # balanço acumulado (Tg C): estimado x net_flux
# res$balance_year     # por ano
# res$flux_monthly     # fluxo médio mensal (g C m-2 d-1) por grupo
# res$plots$scatter; res$plots$balance (total); res$plots$balance_year (anual); res$plots$monthly   # por bioma, em inglês
# save_flux_plots(res, dir = "figs/")                       # salva PNGs
