#' Agrega o XCO2 / fluxo estimado (sondagens) na grade/mês dos fluxos e junta com o fluxo
#'
#' Cada sondagem é atribuída à célula da grade de fluxo mais próxima (centro
#' da célula) e ao mês do seu datetime; depois é agregada por célula x mês e
#' juntada com a tabela de fluxos.
#'
#' Métodos de agregação (`method`):
#'  - "dois_estagios" (padrão): média por célula x DIA e depois média dos dias do
#'    mês. Cada dia pesa igual (sondagens da mesma passagem são quase uma amostra
#'    só) e o nº de amostras independentes é n_days.
#'  - "direto": média simples de todas as sondagens da célula x mês.
#'
#' @param xco2     data.frame/tibble de XCO2 (ex.: saída do mip_ensemble() ou do
#'                 mip_extractor()) ou vetor de caminhos .rds. Precisa de
#'                 lon, lat e datetime (ou time).
#' @param flux     data.frame/tibble de fluxo (saída do flux_extractor()) ou
#'                 vetor de caminhos .rds. Precisa de lon, lat, year, month.
#' @param vars     Colunas a agregar. NULL = detecta entre xco2, xco2_ens_mean,
#'                 xco2_ens_sd, xco2_sim, uncertanty, model_error, fco2, spt.
#'                 A coluna `delta` é sempre agregada se existir; se não existir e
#'                 houver xco2 e xco2_ens_mean, é calculada como
#'                 delta = xco2 - xco2_ens_mean em cada sondagem, antes de agregar.
#' @param method   "dois_estagios" (padrão) ou "direto". Ver acima.
#' @param res      Resolução da grade de fluxo em graus (1 ou 2 valores: lon, lat).
#'                 NULL = inferida a partir das células de `flux`.
#' @param min_n    Nº mínimo de sondagens por célula x mês.
#' @param min_days Nº mínimo de dias distintos com sondagem na célula x mês.
#' @param sigma    Nome(s) da(s) coluna(s) com a incerteza (1 desvio padrão) por
#'                 sondagem, ex.: "sigma_f". Gera a incerteza da média da célula x mês:
#'                 * "dois_estagios": <sigma>_cell. Dentro do dia os erros são
#'                   totalmente correlacionados (σ_dia = média dos σ); entre dias
#'                   são independentes: σ_mês = sqrt(Σ σ_dia²) / n_days.
#'                 * "direto": <sigma>_cell_indep (rho = 0), <sigma>_cell_corr
#'                   (rho = 1) e, se `rho` for informado, <sigma>_cell_rho.
#' @param rho      Só para method = "direto": correlação média dos erros entre
#'                 sondagens da célula x mês (0 a 1).
#' @param by_experiment Se FALSE (padrão), o experimento NÃO entra na junção: o
#'                 delta pode vir de um experimento (ex.: IS) e o fluxo de outro
#'                 (ex.: LNLGIS). O experimento do fluxo fica em `experiment_flux`.
#'                 Se TRUE, junta também por `experiment` (mesmo nome nas duas bases).
#'                 Filtre antes o `flux` para uma submission/experimento só, senão as
#'                 linhas se repetem por submission/experimento.
#' @param join     "inner" (só células/meses com XCO2 e fluxo) ou "left"
#'                 (mantém todas as células/meses com XCO2).
#' @param out_file Se informado, salva o resultado em .rds neste caminho.
#'
#' @return tibble com uma linha por célula x mês (x experimento, se existir nas
#'         duas tabelas): lon, lat (centro da célula), year, month, n_sound,
#'         n_days, <var>_mean, <var>_sd (entre dias, no "dois_estagios"),
#'         <var>_sem ("dois_estagios": sd/sqrt(n_days), checagem empírica da
#'         incerteza da média) e as colunas de fluxo.
agg_xco2_flux <- function(xco2,
                          flux,
                          vars     = NULL,
                          method   = c("dois_estagios", "direto"),
                          res      = NULL,
                          min_n    = 1,
                          min_days = 1,
                          sigma    = NULL,
                          rho      = NULL,
                          by_experiment = FALSE,
                          join     = c("inner", "left"),
                          out_file = NULL) {

  method <- match.arg(method)
  join   <- match.arg(join)
  rd <- function(x) {
    if (is.character(x)) dplyr::bind_rows(lapply(x, readRDS)) else as.data.frame(x)
  }
  xco2 <- rd(xco2)
  flux <- rd(flux)

  # ---- checagens -----------------------------------------------------------
  if (!"datetime" %in% names(xco2)) {
    if (!"time" %in% names(xco2)) stop("xco2 precisa de 'datetime' ou 'time'.")
    xco2$datetime <- as.POSIXct(xco2$time, origin = "1970-01-01", tz = "UTC")
  }
  need_x <- c("lon", "lat", "datetime")
  need_f <- c("lon", "lat", "year", "month")
  if (!all(need_x %in% names(xco2))) stop("xco2 precisa das colunas: ", paste(need_x, collapse = ", "))
  if (!all(need_f %in% names(flux))) stop("flux precisa das colunas: ", paste(need_f, collapse = ", "))

  # ---- grade de fluxo ------------------------------------------------------
  if (is.null(res)) {
    step <- function(v) {
      v <- sort(unique(round(v, 6)))
      if (length(v) < 2) NA_real_ else min(diff(v))
    }
    res <- c(step(flux$lon), step(flux$lat))
    if (anyNA(res)) stop("Não consegui inferir a resolução (flux tem uma só linha/célula). Informe 'res'.")
  }
  if (length(res) == 1) res <- rep(res, 2)

  lon0 <- min(flux$lon)   # centro de célula conhecido (origem da malha)
  lat0 <- min(flux$lat)

  # índices inteiros da célula (evita problemas de ponto flutuante no join)
  xco2$ix <- round((xco2$lon - lon0) / res[1])
  xco2$iy <- round((xco2$lat - lat0) / res[2])
  flux$ix <- round((flux$lon - lon0) / res[1])
  flux$iy <- round((flux$lat - lat0) / res[2])

  xco2$year  <- as.integer(format(xco2$datetime, "%Y"))
  xco2$month <- as.integer(format(xco2$datetime, "%m"))
  xco2$day   <- as.Date(xco2$datetime)

  # ---- variáveis a agregar -------------------------------------------------
  if (is.null(vars)) {
    vars <- intersect(c("xco2", "xco2_ens_mean", "xco2_ens_sd", "xco2_sim",
                        "uncertanty", "model_error", "fco2", "spt"), names(xco2))
  }
  # delta: só é calculado se a coluna não existir (não sobrescreve o delta do usuário)
  if (!"delta" %in% names(xco2) && all(c("xco2", "xco2_ens_mean") %in% names(xco2))) {
    xco2$delta <- xco2$xco2 - xco2$xco2_ens_mean
  }
  if ("delta" %in% names(xco2)) vars <- union(vars, "delta")
  vars <- intersect(vars, names(xco2))
  if (length(vars) == 0) stop("Nenhuma variável para agregar.")
  if (!is.null(sigma) && !all(sigma %in% names(xco2)))
    stop("Coluna(s) de incerteza não encontrada(s) em xco2: ",
         paste(setdiff(sigma, names(xco2)), collapse = ", "))
  sigma_cols <- if (is.null(sigma)) character(0) else sigma

  # ---- chaves de agrupamento/junção ---------------------------------------
  # keys = chaves da junção; gkeys = chaves de agrupamento (mantém o experimento do XCO2)
  keys <- c("ix", "iy", "year", "month")
  if (by_experiment) {
    if (!all(c("experiment") %in% names(xco2)) || !("experiment" %in% names(flux)))
      stop("by_experiment = TRUE exige a coluna 'experiment' em xco2 e em flux.")
    keys <- c(keys, "experiment")
  } else if ("experiment" %in% names(flux)) {
    # o experimento do fluxo pode ser diferente do usado no delta (ex.: delta vs IS,
    # fluxo LNLGIS): não entra na chave e é renomeado para não confundir
    names(flux)[names(flux) == "experiment"] <- "experiment_flux"
  }
  gkeys <- c(keys, if ("experiment" %in% names(xco2) && !"experiment" %in% keys) "experiment")

  if (method == "dois_estagios") {

    vars_d <- setdiff(vars, sigma_cols)

    # estágio 1: média por célula x DIA
    # (dentro do dia os erros são tratados como totalmente correlacionados:
    #  σ_dia = média dos σ das sondagens)
    daily <- xco2 |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(gkeys, "day")))) |>
      dplyr::summarise(
        n_sound_day = dplyr::n(),
        dplyr::across(dplyr::all_of(vars_d), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
        dplyr::across(dplyr::all_of(sigma_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
        .groups = "drop"
      )

    # estágio 2: média dos dias do mês (dias independentes)
    #   σ_mês = sqrt(Σ σ_dia²) / n_days
    agg <- daily |>
      dplyr::group_by(dplyr::across(dplyr::all_of(gkeys))) |>
      dplyr::summarise(
        n_sound = sum(n_sound_day),
        n_days  = dplyr::n(),
        dplyr::across(
          dplyr::all_of(vars_d),
          list(
            mean = ~ mean(.x, na.rm = TRUE),
            sd   = ~ stats::sd(.x, na.rm = TRUE),
            sem  = ~ stats::sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
          ),
          .names = "{.col}_{.fn}"
        ),
        dplyr::across(
          dplyr::all_of(sigma_cols),
          ~ sqrt(sum(.x^2, na.rm = TRUE)) / dplyr::n(),
          .names = "{.col}_cell"
        ),
        .groups = "drop"
      ) |>
      dplyr::filter(n_sound >= min_n, n_days >= min_days)

  } else {

    # agregação direta: todas as sondagens da célula x mês
    agg <- xco2 |>
      dplyr::group_by(dplyr::across(dplyr::all_of(gkeys))) |>
      dplyr::summarise(
        n_sound = dplyr::n(),
        n_days  = dplyr::n_distinct(day),
        dplyr::across(
          dplyr::all_of(vars),
          list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ stats::sd(.x, na.rm = TRUE)),
          .names = "{.col}_{.fn}"
        ),
        # raiz da média dos σ² = incerteza da média com erros totalmente correlacionados
        dplyr::across(
          dplyr::all_of(sigma_cols),
          ~ sqrt(mean(.x^2, na.rm = TRUE)),
          .names = "{.col}_cell_corr"
        ),
        .groups = "drop"
      ) |>
      dplyr::filter(n_sound >= min_n, n_days >= min_days)

    # incerteza da média da célula x mês: σ_média² = σ² [1/n + (1 - 1/n) ρ]
    for (s in sigma_cols) {
      corr <- agg[[paste0(s, "_cell_corr")]]                        # rho = 1
      agg[[paste0(s, "_cell_indep")]] <- corr / sqrt(agg$n_sound)   # rho = 0
      if (!is.null(rho)) {
        agg[[paste0(s, "_cell_rho")]] <-
          corr * sqrt(1 / agg$n_sound + (1 - 1 / agg$n_sound) * rho)
      }
    }
  }

  # ---- junção com o fluxo --------------------------------------------------
  if (anyDuplicated(flux[, keys]) > 0)
    warning("flux tem mais de uma linha por célula x mês (várias submissions/",
            "experimentos?): as linhas do XCO2 serão repetidas. Filtre o flux antes.")
  jf  <- if (join == "inner") dplyr::inner_join else dplyr::left_join
  out <- jf(agg, flux, by = keys)

  if (nrow(out) == 0) {
    warning("Junção vazia: confira se xco2 e flux cobrem a mesma região/período.")
  }

  # células sem fluxo (join = "left"): preenche lon/lat com o centro da célula
  # calculado a partir dos índices, para não perder a localização
  out$lon <- ifelse(is.na(out$lon), lon0 + out$ix * res[1], out$lon)
  out$lat <- ifelse(is.na(out$lat), lat0 + out$iy * res[2], out$lat)

  out <- out |>
    dplyr::select(-ix, -iy) |>
    dplyr::relocate(dplyr::any_of(c("lon", "lat", "year", "month", "experiment",
                                    "experiment_flux", "submission", "n_sound", "n_days"))) |>
    tibble::as_tibble()

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, out_file)
    message("Salvo: ", out_file, " (", format(nrow(out), big.mark = "."), " linhas)")
  }

  out
}

# ---------------------------------------------------------------------------
# Exemplo de uso
# ---------------------------------------------------------------------------
# # df_n já com fco2 e sigma_f por sondagem (fco2_from_delta())
# flux <- readRDS("teste_flux/EnsMean__LNLGIS.rds")        # saída do flux_extractor()
#
# # Padrão: média por dia e depois média dos dias do mês
# base <- agg_xco2_flux(df_n, flux, vars = "fco2", sigma = "sigma_f",
#                       min_n = 5, min_days = 3)
#
# # Para comparar: média direta de todas as sondagens da célula x mês
# base_dir <- agg_xco2_flux(df_n, flux, vars = "fco2", sigma = "sigma_f",
#                           method = "direto", rho = 0.3, min_n = 5, min_days = 3)
#
# # Saída: fco2_mean, fco2_sd, fco2_sem, sigma_f_cell, n_sound, n_days + net_flux etc.
