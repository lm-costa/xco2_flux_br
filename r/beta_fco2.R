#' Beta (tendência de XCO2 por célula) na grade do net_flux do MIP, convertido em fluxo de CO2
#'
#' Reproduz a sua metodologia do beta, ano a ano:
#'   1. soundings -> célula da grade do MIP (1 grau, centros em x.5) x mês
#'      (média, sd, incerteza média, nobs, erro padrão, cv);
#'   2. modelo geral do ano  lm(xco2_mean ~ x)  e regionalização
#'      xco2r = (b0 - delta) - (mean(xco2_mean) - b0),  delta = xco2_est - xco2_mean;
#'   3. beta regional: lm(média mensal de xco2r ~ date)  (beta_r, erro, ilbr, slbr);
#'   4. beta por célula: lm(xco2r ~ date) com >= `min_obs` meses  (ppm por dia);
#'   5. classificação por ano: > Q3 "Source", < Q1 "Sink", resto "Non Significant";
#'   6. conversão:  beta_molm = (1e4 * beta) * 44 / 24.45
#'                  beta_fco2 = beta_molm * 30 / 1000      (g CO2 m-2 month-1)
#'                  beta_fco2_gC_d = beta_fco2 * (12.011/44.009) / 30   (g C m-2 d-1; é o que a função devolve)
#'
#' Unidades (beta em ppm/dia, como lm() sobre datas):
#'   1e4 (altura da coluna, m) * 1e-6 (ppm) * 1e3 (g -> mg) = 1e4 * beta  =>  beta_molm em
#'   mg CO2 m-2 d-1;  x 30 d / 1000  =>  beta_fco2 em g CO2 m-2 month-1.
#'
#' @param xco2   data.frame com uma linha por sounding (use o ensemble/obs, sem repetir
#'               soundings por submission): lon, lat, `xco2_col`, `uncertainty_col` (opcional)
#'               e year/month OU datetime (POSIXct) para derivá-los.
#' @param flux_grid  (opcional) data.frame com lon, lat do net_flux do MIP (ex.: df_flux). Mantém
#'               só as células que existem nessa grade.
#' @param res    Resolução da grade em graus (padrão 1 = grade do MIP).
#' @param years  Anos a processar. NULL = todos presentes.
#' @param xco2_col,uncertainty_col  Nomes das colunas de XCO2 e incerteza.
#' @param dist_col  (opcional) coluna do seu filtro `dist_xco2`; mantém dist < `max_dist`.
#' @param max_dist  Limite do filtro (padrão 0.25).
#' @param min_obs   Mínimo de meses por célula x ano (padrão 5, equivale a n_obs > 4).
#' @param detrend   "index" = exatamente o seu código (x = 1:n() das linhas célula x mês);
#'                  "time" = remove a tendência regional temporal (lm de xco2_mean ~ date, pooled);
#'                  "none" = sem remoção (xco2r = xco2_mean).
#' @param column_height,mw,molar_volume,days_per_month  Constantes da conversão
#'                  (1e4 m, 44 g/mol, 24.45 L/mol, 30 dias).
#' @param linear_reg  (opcional) sua função linear_reg(data, output=); se NULL usa a interna:
#'                  beta1 = inclinação (ppm/dia), p_value, n, betaerror = erro padrão da
#'                  inclinação, modelerror = erro padrão residual.
#'
#' @param full     FALSE (padrão) = devolve só lon, lat, year e beta_fco2_gC_d. TRUE = lista com
#'                  todas as colunas intermediárias (beta_line, erros, p, n_obs, class, beta_molm,
#'                  beta_fco2, betaerror_fco2), `regional` (por ano) e `cell_month`.
#'
#' @return tibble com lon, lat, year e beta_fco2_gC_d (g C m-2 d-1, mesma unidade do net_flux),
#'         uma linha por célula x ano. Com full = TRUE, uma lista (ver `full`).
beta_xco2_flux <- function(xco2,
                           flux_grid       = NULL,
                           res             = 1,
                           years           = NULL,
                           xco2_col        = "xco2",
                           uncertainty_col = "uncertanty",
                           dist_col        = NULL,
                           max_dist        = 0.25,
                           min_obs         = 5,
                           detrend         = c("index", "time", "none"),
                           column_height   = 1e4,
                           mw              = 44,
                           molar_volume    = 24.45,
                           days_per_month  = 30,
                           linear_reg      = NULL,
                           full            = FALSE) {

  detrend <- match.arg(detrend)
  d <- as.data.frame(xco2)
  if (!all(c("lon", "lat", xco2_col) %in% names(d)))
    stop("xco2 precisa de lon, lat e '", xco2_col, "'.")
  if (!all(c("year", "month") %in% names(d))) {
    if (!"datetime" %in% names(d)) stop("Forneça year/month ou datetime.")
    d$year  <- as.integer(format(d$datetime, "%Y", tz = "UTC"))
    d$month <- as.integer(format(d$datetime, "%m", tz = "UTC"))
  }
  d$xco2_v <- d[[xco2_col]]
  d$unc_v  <- if (uncertainty_col %in% names(d)) d[[uncertainty_col]] else NA_real_
  if (!is.null(dist_col)) d <- d[is.finite(d[[dist_col]]) & d[[dist_col]] < max_dist, ]
  d <- d[is.finite(d$xco2_v) & is.finite(d$lon) & is.finite(d$lat), ]

  # ---- célula da grade do MIP (centros em x.5 para res = 1) -----------------
  d$lon <- round(floor(d$lon / res) * res + res / 2, 4)
  d$lat <- round(floor(d$lat / res) * res + res / 2, 4)
  if (!is.null(flux_grid)) {
    g  <- unique(paste(round(flux_grid$lon, 4), round(flux_grid$lat, 4)))
    ok <- paste(d$lon, d$lat) %in% g
    if (!all(ok)) message(sprintf("%.1f%% dos soundings caem fora da grade do net_flux e foram descartados.",
                                  100 * mean(!ok)))
    d <- d[ok, ]
  }
  if (nrow(d) == 0) stop("Nenhum sounding restante.")
  if (is.null(years)) years <- sort(unique(d$year))

  # ---- regressão linear simples (forma fechada) -----------------------------
  lr <- function(date, y) {
    n <- length(y); x <- as.numeric(date)
    sxx <- sum((x - mean(x))^2)
    if (n < 3 || sxx == 0) return(c(beta1 = NA, betaerror = NA, modelerror = NA, p_value = NA, n = n))
    b   <- sum((x - mean(x)) * (y - mean(y))) / sxx
    res <- y - (mean(y) + b * (x - mean(x)))
    s   <- sqrt(sum(res^2) / (n - 2))
    se  <- s / sqrt(sxx)
    p   <- if (se > 0) 2 * stats::pt(-abs(b / se), n - 2) else NA_real_
    c(beta1 = b, betaerror = se, modelerror = s, p_value = p, n = n)
  }
  cell_fit <- function(g) {
    if (is.null(linear_reg)) return(lr(g$date, g$xco2r))
    dat <- data.frame(date = g$date, xco2 = g$xco2r, id_time = g$date)
    c(beta1 = linear_reg(dat, output = "beta1"), betaerror = linear_reg(dat, output = "betaerror"),
      modelerror = linear_reg(dat, output = "modelerror"), p_value = linear_reg(dat, output = "p_value"),
      n = linear_reg(dat, output = "n"))
  }

  out <- list(); reg <- list(); cm <- list()
  for (i in years) {
    di <- d[d$year == i, ]
    if (nrow(di) == 0) next

    # 1) célula x mês
    cmi <- di |>
      dplyr::group_by(lon, lat, year, month) |>
      dplyr::summarise(xco2_mean = mean(xco2_v, na.rm = TRUE),
                       xco2_sd   = stats::sd(xco2_v, na.rm = TRUE),
                       xco2_uncertanty = mean(unc_v, na.rm = TRUE),
                       nobs      = dplyr::n(), .groups = "drop") |>
      dplyr::mutate(xco2_ep = xco2_sd / sqrt(nobs),
                    cv      = 100 * xco2_sd / xco2_mean,
                    date    = as.Date(sprintf("%04d-%02d-15", year, month))) |>
      dplyr::arrange(lon, lat, month)
    cm[[as.character(i)]] <- cmi

    # 2) modelo geral e regionalização
    n <- nrow(cmi)
    if (detrend == "index") {
      x   <- seq_len(n)
      co  <- stats::coef(stats::lm(xco2_mean ~ x, data = cmi))
      est <- co[1] + co[2] * x
      delta <- est - cmi$xco2_mean
      cmi$xco2r <- (co[1] - delta) - (mean(cmi$xco2_mean) - co[1])        # seu código
    } else if (detrend == "time") {
      xt  <- as.numeric(cmi$date)
      b1  <- if (length(unique(xt)) > 1) stats::coef(stats::lm(cmi$xco2_mean ~ xt))[2] else 0
      cmi$xco2r <- cmi$xco2_mean - b1 * (xt - mean(xt))
    } else cmi$xco2r <- cmi$xco2_mean

    # 3) beta regional
    ra <- cmi |> dplyr::group_by(date) |> dplyr::summarise(xco2 = mean(xco2r), .groups = "drop")
    rf <- lr(ra$date, ra$xco2)
    reg[[as.character(i)]] <- data.frame(year = i, beta_r = rf[["beta1"]], ep = rf[["betaerror"]],
                                         ilbr = rf[["beta1"]] - rf[["betaerror"]],
                                         slbr = rf[["beta1"]] + rf[["betaerror"]])

    # 4) beta por célula
    key <- paste(cmi$lon, cmi$lat)
    fits <- lapply(split(cmi, key), cell_fit)
    ft   <- as.data.frame(do.call(rbind, fits))
    cells <- cmi[!duplicated(key), c("lon", "lat")]
    ft$key <- names(fits)
    ft <- ft[match(paste(cells$lon, cells$lat), ft$key), ]
    bi <- data.frame(cells, beta_line = ft$beta1, beta_error = ft$betaerror,
                     model_error = ft$modelerror, p_value = ft$p_value, n_obs = ft$n, year = i)
    bi <- bi[bi$n_obs >= min_obs & is.finite(bi$beta_line), ]
    if (nrow(bi) == 0) next

    # 5) classificação por ano (quartis)
    q3 <- stats::quantile(bi$beta_line, .75); q1 <- stats::quantile(bi$beta_line, .25)
    bi$class <- ifelse(bi$beta_line > q3, "Source", ifelse(bi$beta_line < q1, "Sink", "Non Significant"))
    out[[as.character(i)]] <- bi
  }
  if (length(out) == 0) stop("Nenhuma célula com >= ", min_obs, " meses em nenhum ano.")

  beta <- dplyr::bind_rows(out)
  # 6) conversão para fluxo
  k_molm <- column_height * mw / molar_volume                  # (1e4 * beta) * 44 / 24.45
  beta$beta_molm      <- beta$beta_line * k_molm               # mg CO2 m-2 d-1
  beta$beta_fco2      <- beta$beta_molm * days_per_month / 1000   # g CO2 m-2 month-1
  beta$betaerror_fco2 <- beta$beta_error * k_molm * days_per_month / 1000
  beta$beta_fco2_gC_d <- beta$beta_fco2 * 12.011 / 44.009 / days_per_month   # g C m-2 d-1 (comparável ao net_flux)

  if (!full) return(tibble::as_tibble(beta[, c("lon", "lat", "year", "beta_fco2_gC_d")]))
  list(beta = tibble::as_tibble(beta), regional = dplyr::bind_rows(reg),
       cell_month = tibble::as_tibble(dplyr::bind_rows(cm)))
}

# ---------------------------------------------------------------------------
# Exemplo de uso
# ---------------------------------------------------------------------------
# obs   <- readRDS("teste/ensemble/ensemble__IS.rds")        # 1 linha por sounding (obs)
# flux  <- readRDS("teste_flux/EnsMean__LNLGIS.rds")         # grade do net_flux
# out <- beta_xco2_flux(obs, flux_grid = flux, years = 2015:2022, detrend = "index")
# out   # lon, lat, year, beta_fco2_gC_d  (g C m-2 d-1)
# full <- beta_xco2_flux(obs, flux_grid = flux, full = TRUE)   # com beta_fco2 (g CO2 m-2 month-1), erros etc.

#' Pareia o beta (g C m-2 d-1) com o net_flux do MIP, por célula x ano
#'
#' O beta é a tendência intra-anual do XCO2 em cada célula, convertida em fluxo: equivale a um
#' fluxo MÉDIO do ano (positivo = XCO2 subindo = fonte). Por isso é comparado com a média anual do
#' net_flux, usando só os meses em que a célula tem XCO2 (os mesmos meses do ajuste do beta).
#'
#' @param beta_full  Saída de beta_xco2_flux(..., full = TRUE).
#' @param flux       Fluxo do MIP (df_flux): lon, lat, year, month e `ref`.
#' @param ref,ref_factor  Coluna do net_flux e fator para g C m-2 d-1 (padrão 1).
#' @param anomaly    TRUE = subtrai, em cada ano, a média de todas as células de cada série
#'                   (use com detrend = "time", que remove a tendência regional).
#' @return data.frame (lon, lat, year, month = 1, fco2_mean = beta, net_flux = média anual pareada,
#'         n_months), pronto para map_flux_annual(min_months = 1) e para os scatter/métricas.
pair_beta_flux <- function(beta_full, flux, ref = "net_flux", ref_factor = 1, anomaly = FALSE) {
  cm <- as.data.frame(beta_full$cell_month[, c("lon", "lat", "year", "month")])
  f  <- as.data.frame(flux)[, c("lon", "lat", "year", "month", ref)]
  cm$lon <- round(cm$lon, 4); cm$lat <- round(cm$lat, 4)
  f$lon  <- round(f$lon, 4);  f$lat  <- round(f$lat, 4)
  if (anyDuplicated(f[c("lon", "lat", "year", "month")]))
    stop("`flux` tem linhas repetidas por célula x mês; use um único experimento/submission (ex.: EnsMean).")
  m <- merge(cm, f, by = c("lon", "lat", "year", "month"))
  m <- m[is.finite(m[[ref]]), ]
  r <- stats::aggregate(m[[ref]] * ref_factor, m[c("lon", "lat", "year")], mean)
  names(r)[4] <- "net_flux"
  r$n_months <- stats::aggregate(m[[ref]], m[c("lon", "lat", "year")], length)$x
  b <- as.data.frame(beta_full$beta); b$lon <- round(b$lon, 4); b$lat <- round(b$lat, 4)
  out <- merge(b[, c("lon", "lat", "year", "beta_fco2_gC_d")], r, by = c("lon", "lat", "year"))
  names(out)[names(out) == "beta_fco2_gC_d"] <- "fco2_mean"
  if (anomaly) {
    out$fco2_mean <- ave(out$fco2_mean, out$year, FUN = function(x) x - mean(x))
    out$net_flux  <- ave(out$net_flux,  out$year, FUN = function(x) x - mean(x))
  }
  out$month <- 1L
  out[order(out$year, out$lat, out$lon), c("lon", "lat", "year", "month", "fco2_mean", "net_flux", "n_months")]
}

# Uso:
# full   <- beta_xco2_flux(obs, flux_grid = df_flux, full = TRUE)
# pareado <- pair_beta_flux(full, df_flux)                    # beta x média anual do net_flux
# cor(pareado$fco2_mean, pareado$net_flux, method = "spearman")
# res <- map_flux_annual(pareado, biomes = biomas, south_america = south_america,
#                        map_theme = map_theme, min_months = 1, min_years = 5)
