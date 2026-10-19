#' Distribuição latitudinal do fluxo anual: uma linha por ano
#'
#' Usa `annual` (saída de map_flux_annual(): célula x ano com est_annual e ref_annual, g C m-2 ano-1).
#' Agrupa as células em faixas de latitude (`bin` graus) e tira a média de cada ano; a linha preta
#' grossa é a média de todos os anos. Painéis: estimado, net_flux do MIP e diferença (estimado - MIP).
#' Rótulos em inglês, sem título.
#'
#' @param annual     data.frame `annual` de map_flux_annual() (lon, lat, year, est_annual, ref_annual).
#' @param bin        Largura da faixa de latitude em graus.
#' @param min_cells  Mínimo de células por faixa x ano (faixas com menos são omitidas).
#' @param min_months Mínimo de meses pareados na célula x ano (usa `n_months`, se existir).
#' @param years      Anos a usar. NULL = todos.
#' @param flip       TRUE = latitude no eixo y (perfil vertical, como no mapa); FALSE = latitude no x.
#' @param palette    Escala viridis ("viridis", "plasma", "cividis"...).
#'
#' @return lista: plot (ggplot), data (média por faixa x ano x série).
plot_lat_profile <- function(annual,
                             bin        = 2,
                             min_cells  = 5,
                             min_months = 3,
                             years      = NULL,
                             flip       = TRUE,
                             palette    = "viridis") {
  need <- c("lat", "year", "est_annual", "ref_annual")
  if (!all(need %in% names(annual))) stop("annual precisa de: ", paste(need, collapse = ", "))
  d <- as.data.frame(annual)
  if ("n_months" %in% names(d)) d <- d[d$n_months >= min_months, ]
  if (!is.null(years)) d <- d[d$year %in% years, ]
  d$diff_annual <- d$est_annual - d$ref_annual
  d$latb <- floor(d$lat / bin) * bin + bin / 2

  lab <- c(est_annual = "Delta-based estimate", ref_annual = "MIP net flux", diff_annual = "Difference (estimate - MIP)")
  long <- do.call(rbind, lapply(names(lab), function(v)
    data.frame(latb = d$latb, year = d$year, series = unname(lab[v]), value = d[[v]],
               stringsAsFactors = FALSE)))
  prof <- long |>
    dplyr::group_by(series, year, latb) |>
    dplyr::summarise(flux = mean(value), n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n >= min_cells)
  allyr <- long |>
    dplyr::group_by(series, latb) |>
    dplyr::summarise(flux = mean(value), n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n >= min_cells * 2)
  prof$series  <- factor(prof$series,  levels = unname(lab))
  allyr$series <- factor(allyr$series, levels = unname(lab))
  prof$year <- factor(prof$year)

  p <- ggplot2::ggplot(prof, ggplot2::aes(latb, flux, colour = year, group = year)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.9) +
    ggplot2::geom_line(data = allyr, ggplot2::aes(latb, flux, group = 1), inherit.aes = FALSE,
                       colour = "black", linewidth = 1.3) +
    ggplot2::facet_wrap(~ series, nrow = 1) +
    ggplot2::scale_colour_viridis_d(option = palette, end = 0.92, name = NULL) +
    ggplot2::labs(x = "Latitude", y = expression("Annual flux (g C"~m^{-2}~yr^{-1}*")")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(face = "bold", hjust = 0),
                   legend.position = "right")
  if (flip) {
    p <- p + ggplot2::coord_flip() +
      ggplot2::labs(x = "Latitude", y = expression("Annual flux (g C"~m^{-2}~yr^{-1}*")"))
  }
  list(plot = p, data = prof)
}
