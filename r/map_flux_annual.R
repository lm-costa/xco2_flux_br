#' Mapas de fluxo anual (paleta inferno) e regressão célula a célula estimado x net_flux
#'
#' 1. Fluxo anual por célula: média do fluxo diário dos meses PAREADOS (estimado e net_flux
#'    existentes) do ano x dias do ano, em g C m-2 ano-1. Os dois fluxos usam exatamente os
#'    mesmos meses em cada célula, então são comparáveis. Células x ano com poucos meses
#'    pareados são descartadas (`min_months`).
#' 2. Mapas anuais do estimado e do net_flux (inferno, mesma escala de cores).
#' 3. Regressão célula a célula: em cada célula, estimado anual ~ net_flux anual ao longo dos
#'    anos (precisa de `min_years` anos). Métricas por célula: inclinação, intercepto, r,
#'    R2, p da regressão, rho de Spearman, RMSE, viés. Mapas das métricas e resumos
#'    (todo o Brasil e por bioma).
#'
#' @param df         Saída de agg_xco2_flux() (uma linha por célula x mês), com lon, lat,
#'                   year, month e as colunas `est` e `ref`.
#' @param biomes     sf com os biomas (ex.: geobr::read_biomes(year = 2019)). Desenha as
#'                   fronteiras nos mapas e permite o resumo por bioma. NULL = sem biomas.
#' @param biome_col  Coluna com o nome do bioma (padrão "name_biome").
#' @param drop_outside Com `biomes`, exclui células cujo centro cai fora de todos os biomas.
#' @param est,ref    Colunas do fluxo estimado e do net_flux (g C m-2 d-1).
#' @param ref_factor Fator que converte `ref` para g C m-2 d-1 (padrão 1; confira a unidade).
#' @param years      Anos a usar. NULL = todos.
#' @param min_days   Mínimo de dias com sondagem na célula x mês (usa `n_days`, se existir).
#' @param min_months Mínimo de meses pareados na célula x ano para calcular o fluxo anual.
#' @param min_years  Mínimo de anos por célula para fazer a regressão.
#' @param days_in_year Dias do ano para anualizar o fluxo diário (padrão 365).
#' @param limits     Limites das escalas de cores. NULL = percentis 1% e 99% de CADA mapa (escalas
#'                   independentes). Ou lista com `estimated`, `reference` e/ou `bias`,
#'                   ex.: list(estimated = c(-300, 300), bias = c(-500, 500)).
#' @param south_america sf com o contorno da América do Sul (preenchido de branco). NULL = tenta
#'                   rnaturalearth; se indisponível, omite.
#' @param map_theme  Seu tema de mapa (ex.: map_theme_2) como função ou objeto theme(). NULL = padrão.
#' @param basemap    Desenha o tile `cartolight` (ggspatial::annotation_map_tile; precisa de internet).
#' @param exclude_biomes Biomas que não são desenhados (padrão "Sistema Costeiro").
#' @param xlim,ylim  Janela dos mapas.
#' @param biome_linewidth,biome_colour Contorno dos biomas desenhado POR CIMA das células.
#' @param decorations Adiciona seta norte e barra de escala (ggspatial).
#' @param bias_palette "diverging" (azul-cinza-laranja centrada em 0, padrão) ou "inferno".
#' @param ncol       Colunas de painéis (anos) nos mapas de fluxo.
#' @param res        Resolução da grade em graus (1 ou 2 valores). NULL = inferida.
#' @param plot       Gera os mapas (requer ggplot2).
#'
#' @return lista: annual (célula x ano), regression (por célula), summary (Brasil),
#'         summary_biome (por bioma, se `biomes`), maps (lista de ggplots).
map_flux_annual <- function(df,
                            biomes       = NULL,
                            biome_col    = "name_biome",
                            drop_outside = TRUE,
                            est          = "fco2_mean",
                            ref          = "net_flux",
                            ref_factor   = 1,
                            years        = NULL,
                            min_days     = 1,
                            min_months   = 6,
                            min_years    = 5,
                            days_in_year = 365,
                            limits       = NULL,
                            south_america = NULL,
                            map_theme    = NULL,
                            basemap      = TRUE,
                            exclude_biomes = "Sistema Costeiro",
                            xlim         = c(-75, -35),
                            ylim         = c(-35, 5.5),
                            decorations  = TRUE,
                            biome_linewidth = 0.5,
                            biome_colour = "black",
                            bias_palette = c("diverging", "inferno"),
                            ncol         = 3,
                            res          = NULL,
                            plot         = TRUE) {

  bias_palette <- match.arg(bias_palette)
  need <- c("lon", "lat", "year", "month", est, ref)
  if (!all(need %in% names(df)))
    stop("df precisa das colunas: ", paste(setdiff(need, names(df)), collapse = ", "))
  if (ref_factor == 1)
    message("Assumindo que '", ref, "' já está em g C m-2 d-1 (ref_factor = 1). Confira a unidade no README do MIP.")

  d <- as.data.frame(df)
  d$est <- d[[est]]
  d$ref <- d[[ref]] * ref_factor
  d$lon <- round(d$lon, 4); d$lat <- round(d$lat, 4)
  if ("n_days" %in% names(d)) d <- d[is.na(d$n_days) | d$n_days >= min_days, ]
  if (!is.null(years)) d <- d[d$year %in% years, ]
  d <- d[is.finite(d$est) & is.finite(d$ref), ]
  if (nrow(d) == 0) stop("Nenhuma célula x mês com estimado e net_flux ao mesmo tempo.")

  # ---- biomas ---------------------------------------------------------------
  has_biome <- !is.null(biomes)
  if (has_biome) {
    s2_old <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(s2_old), add = TRUE)
    b   <- sf::st_make_valid(sf::st_transform(biomes, 4326))
    xy  <- unique(d[, c("lon", "lat")])
    pts <- sf::st_as_sf(xy, coords = c("lon", "lat"), crs = 4326)
    hit <- sf::st_intersects(pts, b)
    idx <- vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))
    xy$biome <- as.character(sf::st_drop_geometry(b)[[biome_col]])[idx]
    d <- merge(d, xy, by = c("lon", "lat"), all.x = TRUE, sort = FALSE)
    if (drop_outside) {
      n_out <- length(unique(paste(d$lon, d$lat)[is.na(d$biome)]))
      if (n_out > 0) message(sprintf("%d célula(s) fora dos biomas excluída(s).", n_out))
      d <- d[!is.na(d$biome), ]
      if (nrow(d) == 0) stop("Nenhuma célula dentro dos biomas.")
    } else d$biome[is.na(d$biome)] <- "Outside biomes"
  }

  # ---- fluxo anual por célula (mesmos meses nos dois fluxos) -----------------
  key <- c("lon", "lat", "year", if (has_biome) "biome")
  annual <- d |>
    dplyr::group_by(dplyr::across(dplyr::all_of(key))) |>
    dplyr::summarise(n_months   = dplyr::n(),
                     est_annual = mean(est) * days_in_year,
                     ref_annual = mean(ref) * days_in_year,
                     .groups = "drop") |>
    dplyr::filter(n_months >= min_months) |>
    dplyr::mutate(bias_annual = est_annual - ref_annual)      # estimado - net_flux
  if (nrow(annual) == 0) stop("Nenhuma célula x ano com >= ", min_months, " meses pareados.")

  # ---- regressão célula a célula --------------------------------------------
  reg_one <- function(x) {
    y <- x$est_annual; x <- x$ref_annual; n <- length(x)
    sxx <- sum((x - mean(x))^2); syy <- sum((y - mean(y))^2); sxy <- sum((x - mean(x)) * (y - mean(y)))
    ok  <- n >= 3 && sxx > 0 && syy > 0
    r   <- if (ok) sxy / sqrt(sxx * syy) else NA_real_
    tt  <- if (ok) r * sqrt((n - 2) / (1 - r^2)) else NA_real_
    tibble::tibble(
      n_years   = n,
      mean_ref  = mean(x), mean_est = mean(y),
      slope     = if (sxx > 0) sxy / sxx else NA_real_,
      intercept = if (sxx > 0) mean(y) - sxy / sxx * mean(x) else NA_real_,
      r         = r, r2 = r^2,
      p         = if (ok) 2 * stats::pt(-abs(tt), n - 2) else NA_real_,
      rho       = if (ok) stats::cor(x, y, method = "spearman") else NA_real_,
      rmse      = sqrt(mean((y - x)^2)),
      bias      = mean(y - x))                  # estimado - net_flux
  }
  regression <- annual |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c("lon", "lat", if (has_biome) "biome")))) |>
    dplyr::filter(dplyr::n() >= min_years) |>
    dplyr::group_modify(~ reg_one(.x)) |>
    dplyr::ungroup()
  if (nrow(regression) == 0) stop("Nenhuma célula com >= ", min_years, " anos de dados pareados.")

  # ---- resumos --------------------------------------------------------------
  summ <- function(g) tibble::tibble(
    n_cells        = nrow(g),
    median_r       = stats::median(g$r, na.rm = TRUE),
    median_r2      = stats::median(g$r2, na.rm = TRUE),
    median_rho     = stats::median(g$rho, na.rm = TRUE),
    median_slope   = stats::median(g$slope, na.rm = TRUE),
    median_rmse    = stats::median(g$rmse, na.rm = TRUE),
    median_bias    = stats::median(g$bias, na.rm = TRUE),
    frac_r_positive = mean(g$r > 0, na.rm = TRUE),
    frac_p_lt_0.05  = mean(g$p < 0.05, na.rm = TRUE))
  summary_br <- summ(regression)
  # métricas "empilhadas" (todos os pares célula x ano), como complemento
  ak <- annual[paste(annual$lon, annual$lat) %in% paste(regression$lon, regression$lat), ]
  summary_br$pooled_r    <- stats::cor(ak$ref_annual, ak$est_annual)
  summary_br$pooled_rho  <- stats::cor(ak$ref_annual, ak$est_annual, method = "spearman")
  summary_br$pooled_rmse <- sqrt(mean((ak$est_annual - ak$ref_annual)^2))
  summary_br$pooled_bias <- mean(ak$est_annual - ak$ref_annual)
  summary_biome <- if (has_biome) {
    regression |> dplyr::group_by(biome) |> dplyr::group_modify(~ summ(.x)) |> dplyr::ungroup()
  } else NULL

  # ---- mapas ----------------------------------------------------------------
  maps <- NULL
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    maps <- .make_annual_maps(annual, regression, b = if (has_biome) b else NULL,
                              biome_col = biome_col, limits = limits, ncol = ncol, res = res,
                              south_america = south_america, map_theme = map_theme,
                              basemap = basemap, exclude_biomes = exclude_biomes,
                              xlim = xlim, ylim = ylim, decorations = decorations,
                              biome_linewidth = biome_linewidth, biome_colour = biome_colour,
                              bias_palette = bias_palette)
  }

  list(annual = tibble::as_tibble(annual), regression = regression,
       summary = summary_br, summary_biome = summary_biome, maps = maps)
}

# ---------------------------------------------------------------------------
# Mapas (figuras em inglês, sem título): estilo tile cartolight + contorno + biomas
# ---------------------------------------------------------------------------
.make_annual_maps <- function(annual, regression, b, biome_col, limits, ncol, res,
                              south_america, map_theme, basemap, exclude_biomes,
                              xlim, ylim, decorations, biome_linewidth, biome_colour, bias_palette) {
  gg <- ggplot2::ggplot; aes <- ggplot2::aes
  has_gsp <- requireNamespace("ggspatial", quietly = TRUE)

  step <- function(v) { v <- sort(unique(round(v, 4))); if (length(v) < 2) NA_real_ else min(diff(v)) }
  if (is.null(res)) res <- c(step(regression$lon), step(regression$lat))
  if (anyNA(res)) res <- c(1, 1)
  if (length(res) == 1) res <- rep(res, 2)

  # contorno da América do Sul
  if (is.null(south_america) && requireNamespace("rnaturalearth", quietly = TRUE)) {
    south_america <- try(rnaturalearth::ne_countries(continent = "South America",
                                                     returnclass = "sf"), silent = TRUE)
    if (inherits(south_america, "try-error")) south_america <- NULL
  }
  # biomas desenhados (sem os excluídos)
  b_draw <- NULL
  if (!is.null(b)) {
    b_draw <- b
    if (!is.null(exclude_biomes) && biome_col %in% names(b))
      b_draw <- b[!(as.character(b[[biome_col]]) %in% exclude_biomes), ]
  }

  # remove a coluna `year` das camadas sf: geobr::read_biomes() traz year = 2019, e o facet_wrap(~year)
  # então desenharia os biomas só no painel de 2019 (sem `year`, a camada se repete em todos os painéis)
  no_year <- function(x) if (!is.null(x) && "year" %in% names(x)) x[, setdiff(names(x), "year")] else x
  b_draw <- no_year(b_draw); south_america <- no_year(south_america)

  # camadas de fundo (na ordem do seu código)
  base_layers <- list()
  if (basemap && has_gsp)
    base_layers <- c(base_layers, list(ggspatial::annotation_map_tile(type = "cartolight", progress = "none")))
  else if (basemap)
    message("Pacote 'ggspatial' ausente: mapa sem tile base.")
  if (!is.null(south_america))
    base_layers <- c(base_layers, list(ggplot2::geom_sf(data = south_america, inherit.aes = FALSE,
                                                        col = "grey", fill = "white")))
  if (!is.null(b_draw))
    base_layers <- c(base_layers, list(ggplot2::geom_sf(data = b_draw, inherit.aes = FALSE,
                                                        col = "black", fill = "grey50")))
  # camadas/tema do usuário (map_theme pode devolver theme + annotation_scale/north_arrow)
  th_user <- if (is.function(map_theme)) map_theme() else map_theme
  th_user <- if (is.null(th_user)) list() else if (inherits(th_user, "list")) th_user else list(th_user)
  has_geom <- function(cls) any(vapply(th_user, function(x) inherits(x, "LayerInstance") &&
                                         inherits(x$geom, cls), logical(1)))
  top_border <- if (!is.null(b_draw))
    list(ggplot2::geom_sf(data = b_draw, inherit.aes = FALSE, fill = NA,
                          colour = biome_colour, linewidth = biome_linewidth))
  deco <- list()
  if (decorations && has_gsp) {   # só adiciona o que o seu tema ainda não traz
    if (!has_geom("GeomNorthArrow"))
      deco <- c(deco, list(ggspatial::annotation_north_arrow(location = "tr", which_north = "true",
                                                             style = ggspatial::north_arrow_nautical(),
                                                             height = grid::unit(1.2, "cm"), width = grid::unit(1.2, "cm"))))
    if (!has_geom("GeomScaleBar"))
      deco <- c(deco, list(ggspatial::annotation_scale(location = "bl", width_hint = 0.3,
                                                       unit_category = "metric")))
  }

  # legenda dentro, canto inferior direito (compatível com ggplot2 < e >= 3.5)
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

  build <- function(data, col, scale_fun, name, facet = FALSE) {
    p <- gg(data, aes(lon, lat))
    for (l in base_layers) p <- p + l
    p <- p + ggplot2::geom_tile(aes(color = .data[[col]], fill = .data[[col]]),
                                width = res[1], height = res[2]) +
      scale_fun(name)
    for (l in top_border) p <- p + l      # contorno dos biomas sobre as células
    if (facet) p <- p + ggplot2::facet_wrap(~ year, ncol = ncol)
    for (l in deco) p <- p + l
    p <- p + ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                               crs = sf::st_crs(4326), default_crs = sf::st_crs(4326)) +
      ggplot2::scale_x_continuous(breaks = seq(-70, -40, by = 10)) +
      ggplot2::scale_y_continuous(breaks = seq(-30, 0, by = 10)) +
      ggplot2::labs(x = "Longitude", y = "Latitude") +
      ggplot2::theme_bw()
    for (l in th) if (!is.null(l)) p <- p + l
    p
  }

  # --- escalas ------------------------------------------------------------
  q <- function(v, p = c(0.01, 0.99)) as.numeric(stats::quantile(v, p, na.rm = TRUE))
  lim_of <- function(key, v) if (!is.null(limits[[key]])) limits[[key]] else q(v)
  bar <- ggplot2::guide_colourbar(barheight = grid::unit(4, "cm"))
  inferno <- function(lim) function(name) list(
    ggplot2::scale_color_viridis_c(option = "inferno", limits = lim, oob = scales::squish, name = name, guide = bar),
    ggplot2::scale_fill_viridis_c(option = "inferno", limits = lim, oob = scales::squish, name = name, guide = bar))
  diverging <- function(lim, mid) function(name) list(
    ggplot2::scale_color_gradient2(low = "#2a78d6", mid = "grey93", high = "#eb6834", midpoint = mid,
                                   limits = lim, oob = scales::squish, name = name, guide = bar),
    ggplot2::scale_fill_gradient2(low = "#2a78d6", mid = "grey93", high = "#eb6834", midpoint = mid,
                                  limits = lim, oob = scales::squish, name = name, guide = bar))

  fco2   <- expression('FCO'[2])
  bias_l <- if (!is.null(limits[["bias"]])) limits[["bias"]] else { m <- max(abs(q(annual$bias_annual))); c(-m, m) }

  bias_scale <- if (bias_palette == "diverging") diverging(bias_l, 0) else inferno(bias_l)
  sl <- q(regression$slope, c(0.02, 0.98)); bl <- max(abs(q(regression$bias, c(0.02, 0.98))))
  div_sl <- diverging(c(min(sl[1], 0), max(sl[2], 2)), 1)
  unit <- quote(g ~ C ~ m^{-2} ~ yr^{-1})

  list(
    flux_estimated = build(annual, "est_annual", inferno(lim_of("estimated", annual$est_annual)), fco2, TRUE),
    flux_reference = build(annual, "ref_annual", inferno(lim_of("reference", annual$ref_annual)), fco2, TRUE),
    bias           = build(annual, "bias_annual", bias_scale, "Bias", TRUE),
    r         = build(regression, "r", inferno(c(-1, 1)), "r"),
    rmse      = build(regression, "rmse", inferno(c(0, q(regression$rmse, c(0.02, 0.98))[2])), "RMSE"),
    slope     = build(regression, "slope", div_sl, "Slope"),
    bias_mean = build(regression, "bias", diverging(c(-bl, bl), 0), "Bias")
  )
}

#' Salva os mapas de map_flux_annual() em PNG
#' @param res  Saída de map_flux_annual().
#' @param dir  Pasta de saída.
#' @param dpi  Resolução (padrão 300).
#' @param panel Tamanho (polegadas) de cada painel.
save_flux_maps <- function(res, dir = ".", dpi = 300, panel = 6) {
  stopifnot(!is.null(res$maps))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ny   <- length(unique(res$annual$year))
  ncol <- min(3, ny); nrow <- ceiling(ny / ncol)
  size <- function(nm) if (nm %in% c("flux_estimated", "flux_reference", "bias"))
    c(panel * ncol, panel * nrow) else c(panel, panel)
  for (nm in names(res$maps)) {
    sz <- size(nm)
    ggplot2::ggsave(file.path(dir, paste0("map_", nm, ".png")), res$maps[[nm]],
                    width = sz[1], height = sz[2], dpi = dpi, bg = "white")
  }
  invisible(file.path(dir, paste0("map_", names(res$maps), ".png")))
}

# ---------------------------------------------------------------------------
# Exemplo de uso
# ---------------------------------------------------------------------------
# df_comp <- readRDS("teste_comp.rds")
# biomas  <- geobr::read_biomes(year = 2019)
#
# out <- map_flux_annual(df_comp, biomes = biomas,
#                        south_america = south_america,   # seu objeto sf
#                        map_theme = map_theme_2,         # seu tema
#                        years = 2015:2023,
#                        min_months = 6,   # meses pareados mínimos na célula x ano
#                        min_years  = 5)   # anos mínimos para a regressão da célula
#
# out$summary         # Brasil: medianas de r, slope, rmse, bias; fração de células com p < 0.05
# out$summary_biome   # o mesmo por bioma
# out$regression      # métricas célula a célula (lon, lat, slope, r, p, rmse, bias...)
# out$maps$flux_estimated; out$maps$flux_reference; out$maps$bias; out$maps$r; out$maps$slope
# save_flux_maps(out, "figs/")
