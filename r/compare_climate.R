#' Anomalias climáticas (NASA POWER) x fluxo do MIP e fluxo estimado: mapas e scatter por bioma
#'
#'
#'  1. Coordenadas: usa só as células (lon, lat) presentes em `df_comp` com estimado e net_flux
#'     (não cria mais a grade do Brasil nem usa def_pol).
#'  2. Download (nasapower::get_power, comunidade "ag", diário) só para essas células, de 1 jan do
#'     primeiro ano a 31 dez do último ano de `df_comp`. Cada célula vira um CSV em `data_dir`;
#'     CSVs já existentes são reaproveitados (retoma de onde parou). Tenta até `max_tries` vezes.
#'  3. Pré-processamento mensal por célula: prec (soma), qg, rh, temp, ws (médias).
#'  4. Anomalia só do CLIMA: valor - mediana do período, por célula, com a SÉRIE TODA do NASA
#'     POWER (não só as datas em comum). O fluxo (estimado e net_flux) NÃO vira anomalia: é
#'     comparado direto. Só as datas em que as três séries existem entram na comparação.
#'       - "annual": uma linha por célula x ano (média dos meses do ano, >= `min_months`);
#'         anomalia climática = valor do ano - mediana dos anos da célula; fluxo = média diária
#'         do ano (g C m-2 d-1).
#'       - "monthly": uma linha por célula x mês; anomalia climática = valor - mediana do MESMO mês
#'         do calendário nos anos da célula (senão a sazonalidade domina); fluxo = valor do mês.
#'  5. Comparação: correlação célula a célula (fluxo x anomalia da variável) para o estimado e
#'     para o MIP; mapas dessas correlações e scatter de todas as células x tempo, por bioma
#'     (e Brasil), com r, rho e inclinação.
#'
#' @param df_comp   Saída de agg_xco2_flux(): lon, lat, year, month, `est` e `ref`.
#' @param biomes    sf com os biomas (geobr::read_biomes(year = 2019)). Necessário para o
#'                  scatter por bioma e para os contornos; NULL = só Brasil, sem contornos.
#' @param biome_col Coluna do nome do bioma.
#' @param est,ref   Colunas do fluxo estimado e do MIP (g C m-2 d-1).
#' @param ref_factor Fator que converte `ref` para g C m-2 d-1 (padrão 1).
#' @param vars      Variáveis: subconjunto de "prec", "qg", "rh", "temp", "ws".
#' @param time_scale "annual" (padrão) ou "monthly".
#' @param power_data data.frame já baixado/processado (bruto do POWER com LON, LAT, YEAR, MM, DD,
#'                  ou mensal com cell_lon, cell_lat, year, month e as variáveis). NULL = baixa.
#' @param data_dir  Pasta dos CSVs brutos do POWER (cache do download).
#' @param download  FALSE = não baixa, usa só o que já está em `data_dir`.
#' @param max_tries Tentativas por célula no download.
#' @param pause     Segundos entre tentativas.
#' @param progress  Mostra o progresso do download numa linha: barra, k/total, lon/lat da célula atual,
#'                  quantas faltam, tempo decorrido e ETA (e a tentativa, se houver retentativa).
#' @param min_months Meses mínimos no ano, em cada série, para formar o valor anual ("annual").
#' @param min_years Anos mínimos por célula para a anomalia climática e para a correlação.
#' @param min_days  Mínimo de dias com sondagem na célula x mês (usa `n_days`, se existir).
#' @param years     Anos usados só na COMPARAÇÃO (a anomalia climática usa a série toda do POWER). NULL = todos.
#' @param drop_outside Exclui células fora dos biomas.
#' @param plot,maps Gera scatter e mapas (ggplot2). Os mapas aceitam as mesmas opções de estilo
#'                  de map_flux_annual() em `map_args` (south_america, map_theme, basemap,
#'                  exclude_biomes, xlim, ylim, decorations, biome_linewidth, biome_colour).
#' @param map_args  Lista com opções de estilo dos mapas.
#' @param res       Resolução da grade em graus (1 ou 2 valores). NULL = inferida.
#' @param out_dir   Pasta onde salvar o NASA POWER mensal por célula (prec, qg, rh, temp, ws), antes do
#'                  pareamento com o df_comp e antes da anomalia. NULL = não salva.
#' @param out_file  Nome do arquivo salvo em `out_dir`.
#'
#' @return lista: cells (coordenadas usadas), power_monthly, data (célula x tempo: fluxo e anomalias climáticas), cell_cor (correlação por célula), summary (por bioma, variável, fonte),
#'         pooled (métricas do scatter), plots (scatter por variável), maps (r por variável).
compare_climate <- function(df_comp,
                            biomes       = NULL,
                            biome_col    = "name_biome",
                            est          = "fco2_mean",
                            ref          = "net_flux",
                            ref_factor   = 1,
                            vars         = c("prec", "qg", "rh", "temp", "ws"),
                            time_scale   = c("annual", "monthly"),
                            power_data   = NULL,
                            data_dir     = "nasa_power_data/data-raw",
                            download     = TRUE,
                            max_tries    = 5,
                            pause        = 5,
                            progress     = TRUE,
                            min_months   = 6,
                            min_years    = 5,
                            min_days     = 1,
                            years        = NULL,
                            drop_outside = TRUE,
                            plot         = TRUE,
                            maps         = TRUE,
                            map_args     = list(),
                            res          = NULL,
                            out_dir      = NULL,
                            out_file     = "power_cell_month.rds") {

  time_scale <- match.arg(time_scale)
  vars <- match.arg(vars, c("prec", "qg", "rh", "temp", "ws"), several.ok = TRUE)
  need <- c("lon", "lat", "year", "month", est, ref)
  if (!all(need %in% names(df_comp)))
    stop("df_comp precisa de: ", paste(setdiff(need, names(df_comp)), collapse = ", "))
  src_lab <- c(est = "Delta-based estimate", ref = "MIP net flux")

  # ---- 1) séries completas de cada fonte (as anomalias usam a série toda) -----
  d0 <- as.data.frame(df_comp)
  d0$lon <- round(d0$lon, 4); d0$lat <- round(d0$lat, 4)
  if ("n_days" %in% names(d0)) d0 <- d0[is.na(d0$n_days) | d0$n_days >= min_days, ]
  est_tab <- data.frame(d0[, c("lon", "lat", "year", "month")], est = d0[[est]])
  est_tab <- est_tab[is.finite(est_tab$est), ]
  ref_tab <- data.frame(d0[, c("lon", "lat", "year", "month")], ref = d0[[ref]] * ref_factor)
  ref_tab <- ref_tab[is.finite(ref_tab$ref), ]
  ckey <- function(x) paste(x$lon, x$lat)
  common <- intersect(ckey(est_tab), ckey(ref_tab))
  if (length(common) == 0) stop("Nenhuma célula com estimado e net_flux.")
  est_tab <- est_tab[ckey(est_tab) %in% common, ]; ref_tab <- ref_tab[ckey(ref_tab) %in% common, ]

  has_biome <- !is.null(biomes)
  b <- NULL; cell_biome <- NULL
  if (has_biome) {
    s2_old <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(s2_old), add = TRUE)
    b   <- sf::st_make_valid(sf::st_transform(biomes, 4326))
    xy  <- unique(est_tab[, c("lon", "lat")])
    pts <- sf::st_as_sf(xy, coords = c("lon", "lat"), crs = 4326)
    hit <- sf::st_intersects(pts, b)
    idx <- vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
    xy$biome <- as.character(sf::st_drop_geometry(b)[[biome_col]])[idx]
    if (drop_outside) {
      n_out <- sum(is.na(xy$biome))
      if (n_out > 0) message(sprintf("%d célula(s) fora dos biomas excluída(s).", n_out))
      xy <- xy[!is.na(xy$biome), ]
      if (nrow(xy) == 0) stop("Nenhuma célula dentro dos biomas.")
    } else xy$biome[is.na(xy$biome)] <- "Outside biomes"
    cell_biome <- xy
    est_tab <- est_tab[ckey(est_tab) %in% ckey(xy), ]; ref_tab <- ref_tab[ckey(ref_tab) %in% ckey(xy), ]
  }
  cells <- unique(est_tab[, c("lon", "lat")])
  yr <- range(c(est_tab$year, ref_tab$year))
  message(sprintf("%d células; série do estimado %d-%d, do net_flux %d-%d.", nrow(cells),
                  min(est_tab$year), max(est_tab$year), min(ref_tab$year), max(ref_tab$year)))

  # ---- 2) NASA POWER: download só das células (ou dados fornecidos) ----------
  pw <- if (!is.null(power_data)) .power_prepare(power_data, cells) else {
    if (download) .power_download_cells(cells, yr, data_dir, max_tries, pause, progress)
    .power_read_cells(cells, data_dir)
  }
  if (nrow(pw) == 0) stop("Sem dados do NASA POWER para as células (confira data_dir / download).")

  # ---- 3) mensal por célula (série completa do POWER) --------------------------
  pm <- .power_monthly(pw)
  if (!is.null(out_dir)) {                       # dataset mensal processado, antes do pareamento
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(pm, file.path(out_dir, out_file))
    message("NASA POWER mensal salvo em '", file.path(out_dir, out_file), "'.")
  }
  pm <- pm[ckey(pm) %in% common, c("lon", "lat", "year", "month", vars), drop = FALSE]

  # ---- 4) anomalia só do CLIMA (série toda do POWER); fluxo usado direto; datas comuns --
  keyb <- c("lon", "lat", if (has_biome) "biome")
  tkey <- c("lon", "lat", "year", if (time_scale == "monthly") "month")
  a_est <- .anom_series(est_tab, "est", time_scale, min_months, min_years, anomaly = FALSE)
  a_ref <- .anom_series(ref_tab, "ref", time_scale, min_months, min_years, anomaly = FALSE)
  a_clm <- .anom_series(pm, vars, time_scale, min_months, min_years)      # anomalia
  if (nrow(a_est) == 0 || nrow(a_ref) == 0 || nrow(a_clm) == 0)
    stop("Alguma série ficou vazia: confira min_years / min_months.")
  names(a_clm)[names(a_clm) %in% vars] <- paste0(vars, "_anom")
  anom <- merge(merge(a_est, a_ref, by = tkey), a_clm, by = tkey)   # só datas comuns às 3 séries
  if (!is.null(years)) anom <- anom[anom$year %in% years, ]
  if (has_biome) anom <- merge(anom, cell_biome, by = c("lon", "lat"))
  anom <- anom |> dplyr::group_by(lon, lat) |> dplyr::filter(dplyr::n() >= min_years) |> dplyr::ungroup()
  if (nrow(anom) == 0) stop("Nenhuma data em comum entre estimado, net_flux e NASA POWER.")
  anom <- tibble::as_tibble(anom)
  message(sprintf("%d registros célula x %s em comum.", nrow(anom),
                  if (time_scale == "annual") "ano" else "mês"))

  # ---- 5) correlação célula a célula e métricas por bioma ------------------------
  long <- do.call(rbind, lapply(vars, function(v) do.call(rbind, lapply(c("est", "ref"), function(s)
    data.frame(anom[, c(keyb, "year", if (time_scale == "monthly") "month")],
               variable = v, source = unname(src_lab[s]),
               x = anom[[paste0(v, "_anom")]], y = anom[[s]],
               stringsAsFactors = FALSE)))))
  long <- long[is.finite(long$x) & is.finite(long$y), ]

  cor_one <- function(x, y) {
    n <- length(x)
    if (n < 3 || stats::sd(x) == 0 || stats::sd(y) == 0)
      return(c(n = n, r = NA_real_, p = NA_real_, rho = NA_real_, slope = NA_real_))
    ct <- stats::cor.test(x, y)
    c(n = n, r = unname(ct$estimate), p = ct$p.value,
      rho = stats::cor(x, y, method = "spearman"), slope = stats::cov(x, y) / stats::var(x))
  }
  cell_cor <- long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(keyb, "variable", "source")))) |>
    dplyr::group_modify(~ tibble::as_tibble(as.list(cor_one(.x$x, .x$y)))) |>
    dplyr::ungroup()

  summ <- function(g) tibble::tibble(
    n_cells = nrow(g), median_r = stats::median(g$r, na.rm = TRUE),
    median_rho = stats::median(g$rho, na.rm = TRUE),
    frac_r_positive = mean(g$r > 0, na.rm = TRUE), frac_p_lt_0.05 = mean(g$p < 0.05, na.rm = TRUE))
  gb <- function(x, by) x |> dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::group_modify(~ summ(.x)) |> dplyr::ungroup()
  summary <- dplyr::bind_rows(
    cbind(group = "Brazil", gb(cell_cor, c("variable", "source"))),
    if (has_biome) {
      s <- gb(cell_cor, c("biome", "variable", "source")); names(s)[names(s) == "biome"] <- "group"; s
    })

  # métricas do scatter (todas as células x tempo juntas), Brasil + biomas
  grp_list <- c("Brazil", if (has_biome) sort(unique(long$biome)))
  pooled <- dplyr::bind_rows(lapply(grp_list, function(gp) {
    lg <- if (gp == "Brazil") long else long[long$biome == gp, ]
    lg |> dplyr::group_by(variable, source) |>
      dplyr::group_modify(~ tibble::as_tibble(as.list(cor_one(.x$x, .x$y)))) |>
      dplyr::ungroup() |> dplyr::mutate(group = gp, .before = 1)
  }))

  plots <- NULL; maps_out <- NULL
  if (plot && requireNamespace("ggplot2", quietly = TRUE))
    plots <- .climate_scatter(long, pooled, vars, has_biome, src_lab, time_scale)
  if (maps && requireNamespace("ggplot2", quietly = TRUE)) {
    if (is.null(res)) {
      step <- function(v) { v <- sort(unique(round(v, 4))); if (length(v) < 2) NA_real_ else min(diff(v)) }
      res <- c(step(cells$lon), step(cells$lat))
    }
    maps_out <- do.call(.climate_maps, c(list(cell_cor = cell_cor, vars = vars, b = b,
                                              biome_col = biome_col, res = res,
                                              src_lab = src_lab), map_args))
  }

  list(cells = cells, power_monthly = tibble::as_tibble(pm), data = tibble::as_tibble(anom),
       cell_cor = cell_cor, summary = summary, pooled = pooled, plots = plots, maps = maps_out)
}

# anomalia = valor - mediana da série da célula (série toda, não só datas comuns).
#  annual : média dos meses do ano (>= min_months) e mediana dos anos da célula
#  monthly: mediana do mesmo mês do calendário nos anos da célula
# anomaly = FALSE: só agrega (anual = média dos meses do ano; mensal = valor mensal), sem tirar a mediana.
.anom_series <- function(tab, cols, scale, min_months, min_years, anomaly = TRUE) {
  tab <- tibble::as_tibble(tab)
  if (scale == "annual") {
    tab <- tab |>
      dplyr::group_by(lon, lat, year) |>
      dplyr::filter(dplyr::n() >= min_months) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(cols), mean), .groups = "drop")
    grp <- c("lon", "lat")
  } else grp <- c("lon", "lat", "month")
  if (!anomaly) return(tab)
  tab |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp))) |>
    dplyr::filter(dplyr::n() >= min_years) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(cols), ~ .x - stats::median(.x, na.rm = TRUE))) |>
    dplyr::ungroup()
}

# ---------------------------------------------------------------------------
# NASA POWER: download, leitura e processamento mensal
# ---------------------------------------------------------------------------
.power_pars <- c("ALLSKY_SFC_SW_DWN", "RH2M", "T2M", "PRECTOTCORR", "WS2M", "WD2M")

# nome do arquivo de cada célula: <lon>_<lat>.csv (como no seu power_download)
.power_file <- function(lon, lat, dir) file.path(dir, paste0(lon, "_", lat, ".csv"))

# baixa só as células cujo CSV ainda não existe; tenta `max_tries` vezes (sem loop infinito)
.power_download_cells <- function(cells, yr, dir, max_tries, pause, progress = TRUE) {
  if (!requireNamespace("nasapower", quietly = TRUE))
    stop("Instale o pacote 'nasapower' (ou use download = FALSE / power_data).")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  start <- as.Date(sprintf("%d-01-01", yr[1]))
  end   <- min(as.Date(sprintf("%d-12-31", yr[2])), Sys.Date() - 10)   # POWER tem defasagem de alguns dias
  todo  <- which(!file.exists(.power_file(cells$lon, cells$lat, dir)))
  n <- length(todo)
  message(sprintf("NASA POWER: %d de %d células a baixar (%d já estavam em '%s'); %s a %s.",
                  n, nrow(cells), nrow(cells) - n, dir, start, end))
  if (n == 0) return(invisible(NULL))

  fmt_time <- function(sec) {
    if (!is.finite(sec)) return("--")
    sec <- round(sec); sprintf("%02d:%02d:%02d", sec %/% 3600, (sec %% 3600) %/% 60, sec %% 60)
  }
  # linha única atualizada no lugar: [barra] k/n (%) | lon, lat atuais | faltam | tempo decorrido | ETA
  show <- function(k, i, t0, tag = "") {
    if (!progress) return(invisible(NULL))
    done <- k - 1; el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- if (done > 0) el / done * (n - done) else NA_real_
    w <- 20; fill <- floor(w * done / n)
    cat(sprintf("\r[%s%s] %d/%d (%3.0f%%) | lon %s, lat %s | faltam %d | %s | ETA %s%s   ",
                strrep("#", fill), strrep("-", w - fill), k, n, 100 * done / n,
                cells$lon[i], cells$lat[i], n - k, fmt_time(el), fmt_time(eta), tag))
    utils::flush.console()
  }
  t0 <- Sys.time(); failed <- character(0)
  for (k in seq_len(n)) {
    i <- todo[k]; ok <- FALSE
    for (tr in seq_len(max_tries)) {
      show(k, i, t0, if (tr > 1) sprintf(" | tentativa %d/%d", tr, max_tries) else "")
      dw <- try(nasapower::get_power(community = "ag", lonlat = c(cells$lon[i], cells$lat[i]),
                                     pars = .power_pars, dates = c(start, end),
                                     temporal_api = "daily"), silent = TRUE)
      if (!inherits(dw, "try-error")) { ok <- TRUE; break }
      if (tr < max_tries) Sys.sleep(pause)
    }
    if (ok) {
      utils::write.csv(as.data.frame(dw), .power_file(cells$lon[i], cells$lat[i], dir), row.names = FALSE)
    } else {
      failed <- c(failed, sprintf("(%s, %s)", cells$lon[i], cells$lat[i]))
      if (progress) cat("\n"); message("Falha no download da célula lon ", cells$lon[i], ", lat ", cells$lat[i])
    }
  }
  if (progress) {
    cat(sprintf("\r[%s] %d/%d (100%%) | concluído em %s%s\n", strrep("#", 20), n, n,
                fmt_time(as.numeric(difftime(Sys.time(), t0, units = "secs"))), strrep(" ", 60)))
  }
  message(sprintf("Download concluído: %d ok, %d com falha em %s.", n - length(failed), length(failed),
                  fmt_time(as.numeric(difftime(Sys.time(), t0, units = "secs")))))
  if (length(failed))
    message("Células sem dados: ", paste(failed, collapse = " "), ". Rode de novo para tentar outra vez.")
  invisible(NULL)
}

# lê os CSVs das células; marca a célula pedida (lon/lat do nome do arquivo), não a do ponto
# de grade do POWER (LON/LAT retornados pela API são o ponto mais próximo)
.power_read_cells <- function(cells, dir) {
  out <- lapply(seq_len(nrow(cells)), function(i) {
    f <- .power_file(cells$lon[i], cells$lat[i], dir)
    if (!file.exists(f)) return(NULL)
    x <- utils::read.csv(f)
    x$cell_lon <- cells$lon[i]; x$cell_lat <- cells$lat[i]
    x
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(data.frame())
  message(sprintf("NASA POWER: %d de %d células com dados.", length(out), nrow(cells)))
  dplyr::bind_rows(out)
}

# aceita dado já baixado (bruto diário ou mensal pronto)
.power_prepare <- function(x, cells) {
  x <- as.data.frame(x)
  if (all(c("cell_lon", "cell_lat") %in% names(x))) return(x)
  stop("power_data precisa ter cell_lon e cell_lat (a célula pedida) além das colunas do POWER.")
}

# mensal por célula: prec = soma; qg, rh, temp, ws = médias (igual ao power_pre_processing)
.power_monthly <- function(pw) {
  if (all(c("year", "month", "prec") %in% names(pw)))                     # já mensal
    return(tibble::as_tibble(dplyr::rename(pw, lon = cell_lon, lat = cell_lat)))
  for (v in intersect(.power_pars, names(pw))) pw[[v]][pw[[v]] <= -990] <- NA   # -999 = sem dado
  pw |>
    dplyr::group_by(lon = round(cell_lon, 4), lat = round(cell_lat, 4), year = YEAR, month = MM) |>
    dplyr::summarise(prec = sum(PRECTOTCORR, na.rm = TRUE),
                     qg   = mean(ALLSKY_SFC_SW_DWN, na.rm = TRUE),
                     rh   = mean(RH2M, na.rm = TRUE),
                     temp = mean(T2M, na.rm = TRUE),
                     ws   = mean(WS2M, na.rm = TRUE), .groups = "drop")
}

# ---------------------------------------------------------------------------
# Gráficos (inglês, sem título)
# ---------------------------------------------------------------------------
.climate_labels <- list(
  prec = expression("Precipitation anomaly (mm"~month^{-1}*")"),
  qg   = expression("Shortwave radiation anomaly (MJ"~m^{-2}~d^{-1}*")"),
  rh   = expression("Relative humidity anomaly (%)"),
  temp = expression("Temperature anomaly ("*degree*"C)"),
  ws   = expression("Wind speed anomaly (m"~s^{-1}*")"))

.climate_scatter <- function(long, pooled, vars, has_biome, src_lab, time_scale) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  cols <- c("Delta-based estimate" = "#2a78d6", "MIP net flux" = "#e07b39")
  th <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "top")
  yl <- expression("Flux (g C"~m^{-2}~d^{-1}*")")
  setNames(lapply(vars, function(v) {
    lv <- long[long$variable == v, ]
    lv$group <- "Brazil"
    lg <- lv
    if (has_biome) lg <- rbind(lv, transform(lv, group = biome))
    lg$source <- factor(lg$source, levels = unname(src_lab))
    lg$group <- factor(lg$group, levels = c("Brazil", if (has_biome) sort(unique(lv$biome))))
    pm <- pooled[pooled$variable == v, ]
    pm$group <- factor(pm$group, levels = levels(lg$group))
    pm$source <- factor(pm$source, levels = unname(src_lab))
    pm$lab <- sprintf("r = %.2f", pm$r); pm$lab[!is.finite(pm$r)] <- "r = NA"
    pm$vj  <- ifelse(pm$source == src_lab[["est"]], 1.3, 2.6)
    gg(lg, aes(x, y, colour = source)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
      ggplot2::geom_point(alpha = 0.25, size = 0.9) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
      ggplot2::geom_text(data = pm, aes(x = -Inf, y = Inf, label = lab, vjust = vj, colour = source),
                         hjust = -0.1, size = 3.2, inherit.aes = FALSE, show.legend = FALSE) +
      ggplot2::facet_wrap(~ group, scales = "free") +
      ggplot2::scale_colour_manual(values = cols) +
      ggplot2::labs(x = .climate_labels[[v]], y = yl, colour = NULL) + th
  }), vars)
}

# mapas das correlações célula a célula (estilo de map_flux_annual: tile, contorno, biomas)
.climate_maps <- function(cell_cor, vars, b, biome_col, res, src_lab,
                          south_america = NULL, map_theme = NULL, basemap = TRUE,
                          exclude_biomes = "Sistema Costeiro", xlim = c(-75, -35),
                          ylim = c(-35, 5.5), decorations = TRUE, biome_linewidth = 0.5,
                          biome_colour = "black") {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  has_gsp <- requireNamespace("ggspatial", quietly = TRUE)
  if (anyNA(res)) res <- c(1, 1)
  if (length(res) == 1) res <- rep(res, 2)
  if (is.null(south_america) && requireNamespace("rnaturalearth", quietly = TRUE)) {
    south_america <- try(rnaturalearth::ne_countries(continent = "South America",
                                                     returnclass = "sf"), silent = TRUE)
    if (inherits(south_america, "try-error")) south_america <- NULL
  }
  b_draw <- NULL
  if (!is.null(b)) {
    b_draw <- b
    if (!is.null(exclude_biomes) && biome_col %in% names(b))
      b_draw <- b[!(as.character(b[[biome_col]]) %in% exclude_biomes), ]
  }
  no_year <- function(x) if (!is.null(x) && "year" %in% names(x)) x[, setdiff(names(x), "year")] else x
  b_draw <- no_year(b_draw); south_america <- no_year(south_america)

  base_layers <- list()
  if (basemap && has_gsp)
    base_layers <- c(base_layers, list(ggspatial::annotation_map_tile(type = "cartolight", progress = "none")))
  if (!is.null(south_america))
    base_layers <- c(base_layers, list(ggplot2::geom_sf(data = south_america, inherit.aes = FALSE,
                                                        col = "grey", fill = "white")))
  if (!is.null(b_draw))
    base_layers <- c(base_layers, list(ggplot2::geom_sf(data = b_draw, inherit.aes = FALSE,
                                                        col = "black", fill = "grey50")))
  th_user <- if (is.function(map_theme)) map_theme() else map_theme
  th_user <- if (is.null(th_user)) list() else if (inherits(th_user, "list")) th_user else list(th_user)
  has_geom <- function(cls) any(vapply(th_user, function(x) inherits(x, "LayerInstance") &&
                                         inherits(x$geom, cls), logical(1)))
  top_border <- if (!is.null(b_draw))
    list(ggplot2::geom_sf(data = b_draw, inherit.aes = FALSE, fill = NA,
                          colour = biome_colour, linewidth = biome_linewidth))
  deco <- list()
  if (decorations && has_gsp) {
    if (!has_geom("GeomNorthArrow"))
      deco <- c(deco, list(ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                                             style = ggspatial::north_arrow_nautical(),
                                                             height = grid::unit(1.2, "cm"), width = grid::unit(1.2, "cm"))))
    if (!has_geom("GeomScaleBar"))
      deco <- c(deco, list(ggspatial::annotation_scale(location = "bl", width_hint = 0.3,
                                                       unit_category = "metric")))
  }
  legend_in <- if (utils::packageVersion("ggplot2") >= "3.5.0")
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(1, 0),
                   legend.justification = c(1, 0))
  else ggplot2::theme(legend.position = c(1, 0), legend.justification = c(1, 0))
  th <- c(th_user, list(
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 20),
                   text = ggplot2::element_text(size = 20),
                   legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.8), colour = NA)),
    legend_in))
  bar <- ggplot2::guide_colourbar(barheight = grid::unit(4, "cm"))

  setNames(lapply(vars, function(v) {
    dm <- cell_cor[cell_cor$variable == v & is.finite(cell_cor$r), ]
    dm$source <- factor(dm$source, levels = unname(src_lab))
    p <- gg(dm, aes(lon, lat))
    for (l in base_layers) p <- p + l
    p <- p + ggplot2::geom_tile(aes(color = r, fill = r), width = res[1], height = res[2]) +
      ggplot2::scale_color_gradient2(low = "#2a78d6", mid = "grey93", high = "#eb6834", midpoint = 0,
                                     limits = c(-1, 1), name = "r", guide = bar) +
      ggplot2::scale_fill_gradient2(low = "#2a78d6", mid = "grey93", high = "#eb6834", midpoint = 0,
                                    limits = c(-1, 1), name = "r", guide = bar)
    for (l in top_border) p <- p + l
    p <- p + ggplot2::facet_wrap(~ source, ncol = 2)
    for (l in deco) p <- p + l
    p <- p + ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                               crs = sf::st_crs(4326), default_crs = sf::st_crs(4326)) +
      ggplot2::scale_x_continuous(breaks = seq(-70, -40, by = 10)) +
      ggplot2::scale_y_continuous(breaks = seq(-30, 0, by = 10)) +
      ggplot2::labs(x = "Longitude", y = "Latitude") + ggplot2::theme_bw()
    for (l in th) if (!is.null(l)) p <- p + l
    p
  }), vars)
}

#' Salva os gráficos de compare_climate() em PNG (scatter e mapas, um arquivo por variável)
save_climate_plots <- function(res, dir = ".", dpi = 300, panel = 6) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (v in names(res$plots))
    ggplot2::ggsave(file.path(dir, paste0("climate_scatter_", v, ".png")), res$plots[[v]],
                    width = 11, height = 8, dpi = dpi, bg = "white")
  for (v in names(res$maps))
    ggplot2::ggsave(file.path(dir, paste0("climate_map_r_", v, ".png")), res$maps[[v]],
                    width = panel * 2, height = panel, dpi = dpi, bg = "white")
  invisible(dir)
}

# ---------------------------------------------------------------------------
# Exemplo de uso
# ---------------------------------------------------------------------------
# source("setup_packages.R"); setup_packages(extra = "nasapower")
# df_comp <- readRDS("teste_comp.rds")
# biomas  <- geobr::read_biomes(year = 2019)
#
# res <- compare_climate(df_comp, biomes = biomas,
#                        time_scale = "annual",                  # ou "monthly"
#                        data_dir = "nasa_power_data/data-raw",  # cache do download (retoma se parar)
#                        map_args = list(south_america = south_america, map_theme = map_theme_2))
# res$summary                # mediana de r por bioma, variável e fonte (estimado x MIP)
# res$pooled                 # r, rho, inclinação do scatter geral por bioma
# res$plots$prec             # scatter: fluxo x anomalia da precipitação, por bioma
# res$maps$temp              # mapa de r (estimado e MIP lado a lado) com a temperatura
# save_climate_plots(res, "figs_climate/")
#
# # sem baixar de novo (só o que já está na pasta):
# res <- compare_climate(df_comp, biomes = biomas, download = FALSE)

#' Mapas das medias anuais das variaveis climaticas (NASA POWER) nas celulas de 1 grau do fluxo
#' Usa o dataset salvo por compare_climate(out_dir = ...) (prec, qg, rh, temp, ws por celula x mes).
#' Media dos valores anuais em `years`: precipitacao = soma mensal (mm ano-1); demais variaveis = media mensal.
#' @param power data frame lon, lat, year, month, prec, qg, rh, temp, ws
#' @param biomes sf dos biomas; @param south_america sf de fundo (opcional)
#' @return lista: plot (patchwork), cells (media por celula) e by_biome (media por bioma, 1 linha por variavel)
plot_climate_maps <- function(power, biomes, south_america = NULL, years = 2015:2023, out = "figs_climate/climate_maps.png", width = 9, height = 10) {
  library(ggplot2); library(patchwork)
  p <- power[power$year %in% years, ]
  full <- stats::aggregate(month ~ lon + lat + year, p, length); full <- full[full$month == 12, c("lon", "lat", "year")]
  p <- merge(p, full, by = c("lon", "lat", "year"))                           # so anos com os 12 meses
  ann <- stats::aggregate(cbind(prec = prec, qg = qg, rh = rh, temp = temp, ws = ws) ~ lon + lat + year, p,
                          function(x) c(m = mean(x), s = sum(x)))
  an <- data.frame(lon = ann$lon, lat = ann$lat, year = ann$year, prec = ann$prec[, "s"], qg = ann$qg[, "m"], rh = ann$rh[, "m"], temp = ann$temp[, "m"], ws = ann$ws[, "m"])
  cells <- stats::aggregate(cbind(prec, qg, rh, temp, ws) ~ lon + lat, an, mean)
  b <- sf::st_make_valid(biomes[biomes$name_biome != "Sistema Costeiro", ])
  pts <- sf::st_as_sf(cells, coords = c("lon", "lat"), crs = sf::st_crs(b), remove = FALSE)
  suppressMessages(sf::sf_use_s2(FALSE)); idx <- sf::st_intersects(pts, b)
  cells$biome <- vapply(idx, function(i) if (length(i)) as.character(sf::st_drop_geometry(b)$name_biome[i[1]]) else NA_character_, "")
  vars <- list(rh = list(lab = "Relative humidity (%)", col = "mako", dir = -1),
               prec = list(lab = expression("Precipitation (mm yr"^-1*")"), col = "viridis", dir = -1),
               temp = list(lab = "Air temperature (°C)", col = "inferno", dir = 1),
               qg = list(lab = expression("Shortwave radiation (MJ m"^-2*" d"^-1*")"), col = "plasma", dir = 1),
               ws = list(lab = expression("Wind speed at 2 m (m s"^-1*")"), col = "cividis", dir = 1))
  tags <- c("(a)", "(b)", "(c)", "(d)", "(e)")
  mk <- function(v, i) {
    g <- ggplot() + (if (!is.null(south_america)) geom_sf(data = south_america, fill = "grey92", colour = "grey70", linewidth = 0.2)) +
      geom_tile(data = cells, aes(lon, lat, fill = .data[[v]]), width = 1, height = 1) + geom_sf(data = b, fill = NA, colour = "black", linewidth = 0.25) +
      scale_fill_viridis_c(option = vars[[v]]$col, direction = vars[[v]]$dir, name = vars[[v]]$lab, guide = guide_colourbar(barheight = grid::unit(2.6, "cm"), barwidth = grid::unit(0.35, "cm"))) +
      coord_sf(xlim = c(-74.5, -34), ylim = c(-34, 5.5), expand = FALSE) + labs(title = tags[i], x = NULL, y = NULL) +
      theme_bw(base_size = 9, base_family = "Times New Roman") + theme(plot.title = element_text(face = "bold", size = 10, hjust = 0), legend.title = element_text(size = 8),
            legend.text = element_text(size = 7.5), panel.grid.major = element_line(colour = "grey85", linewidth = 0.2), plot.margin = margin(2, 2, 2, 2))
    if (i == 1) g <- g + ggspatial::annotation_north_arrow(location = "tr", style = ggspatial::north_arrow_minimal, height = grid::unit(0.8, "cm"), width = grid::unit(0.8, "cm")) +
      ggspatial::annotation_scale(location = "bl", height = grid::unit(0.15, "cm"), text_cex = 0.6)
    g
  }
  pl <- Map(mk, names(vars), seq_along(vars))
  fig <- (pl[[1]] + pl[[2]]) / (pl[[3]] + pl[[4]]) / (pl[[5]] + patchwork::plot_spacer())
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  ragg::agg_png(out, width = width, height = height, units = "in", res = 300); print(fig); grDevices::dev.off()
  bb <- do.call(rbind, lapply(c("rh", "prec", "temp", "qg", "ws"), function(v) { a <- stats::aggregate(cells[[v]], list(biome = cells$biome), mean); data.frame(var = v, a) }))
  bw <- stats::reshape(bb, idvar = "biome", timevar = "var", direction = "wide"); names(bw) <- sub("^x\\.", "", names(bw))
  list(plot = fig, cells = tibble::as_tibble(cells), by_biome = tibble::as_tibble(bw))
}
