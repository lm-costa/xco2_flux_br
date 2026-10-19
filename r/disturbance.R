#' Baixa os GeoTIFFs mensais de área queimada do INPE (AQ1KM, Coleção 2)
#'
#' Fonte: https://dataserver-coids.inpe.br/queimadas/queimadas/area_queimada/colecao2/tif/
#' (arquivos `YYYY_MM_AQM_modis.tif`, pixel de 0,01°, 0 = não queimado, 1 = queimado; ~0,8-1 MB cada).
#' Arquivos que já existem em `dir` não são baixados de novo.
#'
#' @param dir      Pasta de destino (criada se não existir).
#' @param years    Anos a baixar.
#' @param months   Meses a baixar (1 a 12).
#' @param base_url Pasta do servidor com os GeoTIFFs.
#' @param max_tries Tentativas por arquivo.
#' @param pause    Segundos entre arquivos.
#'
#' @return Caminhos dos arquivos baixados ou já existentes (invisível); avisa os que falharam.
aq1km_download <- function(dir      = "data-raw/aq1km",
                           years    = 2015:2023,
                           months   = 1:12,
                           base_url = "https://dataserver-coids.inpe.br/queimadas/queimadas/area_queimada/colecao2/tif",
                           max_tries = 3,
                           pause     = 0.5) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  nm   <- as.vector(outer(years, months, function(y, m) sprintf("%d_%02d_AQM_modis.tif", y, m)))
  nm   <- sort(nm)
  dest <- file.path(dir, nm)
  todo <- which(!(file.exists(dest) & file.size(dest) > 0))
  message(sprintf("aq1km_download: %d arquivos pedidos, %d já existem, %d a baixar.",
                  length(nm), length(nm) - length(todo), length(todo)))
  failed <- character(0)
  for (i in todo) {
    ok <- FALSE
    for (tr in seq_len(max_tries)) {
      ok <- tryCatch({ suppressWarnings(utils::download.file(file.path(base_url, nm[i]), dest[i],
                                                             quiet = TRUE, mode = "wb")); TRUE },
                     error = function(e) FALSE)
      if (ok && file.size(dest[i]) > 0) break
      ok <- FALSE; if (file.exists(dest[i])) unlink(dest[i])
      Sys.sleep(pause * tr)
    }
    if (!ok) failed <- c(failed, nm[i])
    Sys.sleep(pause)
  }
  if (length(failed)) message("Falharam: ", paste(failed, collapse = ", "))
  invisible(dest[file.exists(dest)])
}

#' Área queimada mensal do INPE (AQ1KM, GeoTIFF 0/1) agregada na grade do net_flux: célula x mês
#'
#' Cada GeoTIFF mensal (`YYYY_MM_AQM_modis.tif`, pixel de 0,01°, 1 = queimado) é recortado pelos
#' biomas (só pixels dentro do Brasil) e a área queimada (km2, corrigida pela latitude de cada
#' pixel) é somada por célula de `res` graus (centros em x.5 para res = 1). Cada pixel é atribuído
#' à célula pelo seu centro. Os pixels de centro exatamente na borda da célula (os centros do
#' AQ1KM caem em graus inteiros) vão para a célula ao norte/leste.
#'
#' `monthly` é a área queimada de cada mês. `annual` e `total_biome` usam a UNIÃO dos meses do ano
#' (um pixel queimado em mais de um mês do mesmo ano conta uma vez); a soma simples dos meses
#' fica 8-11% acima nesse produto.
#'
#' Atenção: o AQ1KM (MODIS, 1 km) dá bem mais área queimada que mapas de 30 m (por exemplo o
#' MapBiomas Fogo), em particular nos biomas com fogo pequeno e fragmentado. Os valores em km2
#' não são comparáveis entre os produtos.
#'
#' @param dir        Pasta com os GeoTIFFs mensais.
#' @param years      Anos a usar. NULL = todos os arquivos da pasta.
#' @param biomes     sf com os biomas (geobr::read_biomes(year = 2019)); só pixels dentro deles contam.
#' @param biome_col  Coluna com o nome do bioma.
#' @param res        Tamanho da célula em graus.
#' @param out_dir    Pasta onde salvar os resultados. NULL = não salva.
#' @param out_file   Nome do arquivo da tabela célula x mês em `out_dir`.
#'
#' @return lista: monthly (lon, lat, year, month, burned_km2), annual (lon, lat, year, burned_km2),
#'         land (lon, lat, area_land_km2: área da célula dentro dos biomas),
#'         total_biome (km2 queimados por bioma x ano).
aq1km_aggregate <- function(dir,
                            years     = NULL,
                            biomes,
                            biome_col = "name_biome",
                            res       = 1,
                            out_dir   = NULL,
                            out_file  = "aq1km_cell_month.rds") {
  fs <- list.files(dir, pattern = "^\\d{4}_\\d{2}_AQM_modis\\.tif$", full.names = TRUE)
  yr <- as.integer(substr(basename(fs), 1, 4)); mo <- as.integer(substr(basename(fs), 6, 7))
  if (!is.null(years)) { keep <- yr %in% years; fs <- fs[keep]; yr <- yr[keep]; mo <- mo[keep] }
  if (length(fs) == 0) stop("Nenhum GeoTIFF 'YYYY_MM_AQM_modis.tif' em ", dir)
  message(sprintf("aq1km_aggregate: %d arquivos (%d-%d).", length(fs), min(yr), max(yr)))

  tmpl <- terra::rast(fs[1])
  b <- sf::st_make_valid(sf::st_transform(biomes, 4326))
  b$.id <- seq_len(nrow(b))
  bid  <- terra::rasterize(terra::vect(b), tmpl, field = ".id")          # índice do bioma por pixel
  biome_of <- as.character(sf::st_drop_geometry(b)[[biome_col]])
  area_km2 <- terra::cellSize(tmpl, unit = "km")                         # área de cada pixel
  ids  <- terra::values(bid, mat = FALSE)
  ark  <- terra::values(area_km2, mat = FALSE)
  in_br <- which(!is.na(ids))

  cell_of <- function(cells) {
    xy <- terra::xyFromCell(tmpl, cells)
    data.frame(lon = round(floor(xy[, 1] / res) * res + res / 2, 4),
               lat = round(floor(xy[, 2] / res) * res + res / 2, 4))
  }
  # área de terra (dentro dos biomas) por célula, calculada uma vez
  land <- cbind(cell_of(in_br), a = ark[in_br]) |>
    dplyr::group_by(lon, lat) |>
    dplyr::summarise(area_land_km2 = sum(a), .groups = "drop")

  monthly <- vector("list", length(fs)); annual <- list(); tot <- list()
  for (y in sort(unique(yr))) {
    any_y <- rep(FALSE, length(ids))                       # união dos meses do ano (pixel contado 1 vez)
    for (i in which(yr == y)) {
      v <- terra::values(terra::rast(fs[i]), mat = FALSE) == 1
      v[is.na(v)] <- FALSE
      v <- v & !is.na(ids)
      any_y <- any_y | v
      cells <- which(v)
      if (length(cells) == 0) next
      monthly[[i]] <- cbind(cell_of(cells), a = ark[cells]) |> dplyr::group_by(lon, lat) |>
        dplyr::summarise(burned_km2 = sum(a), .groups = "drop") |>
        dplyr::mutate(year = yr[i], month = mo[i], .before = 1)
    }
    cells <- which(any_y)
    d <- cbind(cell_of(cells), a = ark[cells], bio = biome_of[ids[cells]])
    annual[[as.character(y)]] <- d |> dplyr::group_by(lon, lat) |>
      dplyr::summarise(burned_km2 = sum(a), .groups = "drop") |> dplyr::mutate(year = y, .before = 1)
    tot[[as.character(y)]] <- d |> dplyr::group_by(biome = bio) |>
      dplyr::summarise(burned_km2 = sum(a), .groups = "drop") |> dplyr::mutate(year = y, .before = 1)
  }
  monthly <- dplyr::bind_rows(monthly) |> dplyr::select(lon, lat, year, month, burned_km2)
  annual  <- dplyr::bind_rows(annual)  |> dplyr::select(lon, lat, year, burned_km2)
  total_biome <- dplyr::bind_rows(tot) |> dplyr::select(biome, year, burned_km2)

  out <- list(monthly = tibble::as_tibble(monthly), annual = tibble::as_tibble(annual),
              land = tibble::as_tibble(land), total_biome = tibble::as_tibble(total_biome))
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file.path(out_dir, out_file))
    message("aq1km_aggregate: salvo em '", file.path(out_dir, out_file), "'.")
  }
  out
}

#' Desmatamento anual do PRODES (INPE/TerraBrasilis, WFS) agregado na grade do net_flux
#'
#' Lê os polígonos de incremento anual do PRODES por WFS, em páginas, e soma `area_km` por célula de
#' `res` graus (centros em x.5 para res = 1), atribuindo cada polígono à célula do seu centroide.
#' Os polígonos NÃO são guardados: cada página é baixada para um arquivo temporário, agregada e
#' apagada. Cada camada x ano agregada é salva em `out_dir` (retoma de onde parou).
#'
#' `year` é o ano PRODES (1 de agosto do ano anterior a 31 de julho do ano indicado).
#'
#' @param layers   Vetor nomeado: nome = rótulo do bioma, valor = camada WFS ("workspace:layer").
#' @param years    Anos PRODES.
#' @param res      Tamanho da célula em graus.
#' @param page     Feições por requisição.
#' @param pause    Segundos entre requisições (para não sobrecarregar o servidor).
#' @param out_dir  Pasta onde guardar as tabelas por camada x ano e o resultado final.
#' @param url      Endpoint WFS do TerraBrasilis.
#' @param max_tries Tentativas por página.
#'
#' @return lista: cell (biome, lon, lat, year, def_km2, n_poly), total (biome, year, def_km2, n_poly).
prodes_aggregate <- function(layers = c("Amazônia"       = "prodes-amazon-nb:yearly_deforestation_biome",
                                         "Cerrado"        = "prodes-cerrado-nb:yearly_deforestation",
                                         "Caatinga"       = "prodes-caatinga-nb:yearly_deforestation",
                                         "Mata Atlântica" = "prodes-mata-atlantica-nb:yearly_deforestation",
                                         "Pampa"          = "prodes-pampa-nb:yearly_deforestation",
                                         "Pantanal"       = "prodes-pantanal-nb:yearly_deforestation"),
                             years     = 2015:2023,
                             res       = 1,
                             page      = 10000,
                             pause     = 1,
                             out_dir   = "data/prodes",
                             url       = "https://terrabrasilis.dpi.inpe.br/geoserver/ows",
                             max_tries = 4) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  s2_old <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(s2_old), add = TRUE)
  q <- function(...) paste0(url, "?", paste(names(list(...)), vapply(list(...), function(x) utils::URLencode(as.character(x), reserved = TRUE), ""),
                                             sep = "=", collapse = "&"))

  # Algumas camadas do GeoServer dão erro 400 quando se pede só algumas colunas (propertyName);
  # nesses casos pede todas as colunas (payload maior).
  use_prop <- new.env()
  one_page <- function(layer, y, start) {
    tmp <- tempfile(fileext = ".json"); on.exit(unlink(tmp), add = TRUE)
    for (tr in seq_len(max_tries)) {
      prop <- if (is.null(use_prop[[layer]])) TRUE else use_prop[[layer]]
      args <- list(service = "WFS", version = "2.0.0", request = "GetFeature", typeNames = layer,
                   outputFormat = "application/json", srsName = "EPSG:4326",
                   CQL_FILTER = paste0("year=", y), sortBy = "fid", count = page, startIndex = start)
      if (prop) args$propertyName <- "fid,year,area_km,geom"
      u <- do.call(q, args)
      ok <- tryCatch({ suppressWarnings(utils::download.file(u, tmp, quiet = TRUE, mode = "wb")); TRUE },
                     error = function(e) FALSE)
      if (ok && file.exists(tmp) && file.size(tmp) > 0) {
        g <- tryCatch(sf::read_sf(tmp, quiet = TRUE), error = function(e) NULL)
        if (!is.null(g)) { if (nrow(g) > 0) use_prop[[layer]] <- prop; return(g) }   # página vazia não prova nada
      }
      if (prop) { use_prop[[layer]] <- FALSE; next }   # erro com propertyName: tenta com todas as colunas
      Sys.sleep(pause * tr * 2)
    }
    stop("Falha ao baixar a página (", layer, ", ano ", y, ", início ", start, ").")
  }

  cells <- list(); tots <- list()
  for (bio in names(layers)) for (y in years) {
    f <- file.path(out_dir, sprintf("%s_%d.rds", gsub("[^A-Za-z]", "", iconv(bio, to = "ASCII//TRANSLIT")), y))
    if (file.exists(f)) { r <- readRDS(f) } else {
      acc <- list(); start <- 0; npoly <- 0
      repeat {
        g <- one_page(layers[[bio]], y, start)
        if (nrow(g) == 0) break
        # centro do polígono = média dos vértices (robusto a geometrias inválidas, sem make_valid)
        co  <- sf::st_coordinates(sf::st_geometry(g))
        idx <- co[, ncol(co)]
        sx  <- rowsum(co[, 1], idx); sy <- rowsum(co[, 2], idx)
        nv  <- as.vector(table(idx))
        ids <- as.integer(rownames(sx))
        cx  <- as.vector(sx) / nv; cy <- as.vector(sy) / nv
        d <- data.frame(lon = round(floor(cx / res) * res + res / 2, 4),
                        lat = round(floor(cy / res) * res + res / 2, 4),
                        a = g$area_km[ids])
        acc[[length(acc) + 1]] <- d
        npoly <- npoly + nrow(g)
        if (nrow(g) < page) break
        start <- start + page; Sys.sleep(pause)
      }
      d <- if (length(acc)) dplyr::bind_rows(acc) else data.frame(lon = numeric(), lat = numeric(), a = numeric())
      r <- list(cell = d |> dplyr::group_by(lon, lat) |>
                  dplyr::summarise(def_km2 = sum(a, na.rm = TRUE), n_poly = dplyr::n(), .groups = "drop") |>
                  dplyr::mutate(biome = bio, year = y, .before = 1),
                total = data.frame(biome = bio, year = y, def_km2 = sum(d$a, na.rm = TRUE), n_poly = npoly))
      saveRDS(r, f)
      message(sprintf("[%s] %s %d: %d polígonos, %.0f km2", format(Sys.time(), "%H:%M:%S"), bio, y, npoly,
                      r$total$def_km2))
      Sys.sleep(pause)
    }
    cells[[length(cells) + 1]] <- r$cell; tots[[length(tots) + 1]] <- r$total
  }
  out <- list(cell = tibble::as_tibble(dplyr::bind_rows(cells)), total = tibble::as_tibble(dplyr::bind_rows(tots)))
  saveRDS(out, file.path(out_dir, "prodes_cell_year.rds"))
  out
}

#' Mapas anuais de área queimada ou desmatamento agregados nas células do fluxo (escala em Mha)
#'
#' Um painel por ano (e, opcionalmente, o acumulado do período) com a área por célula de 1°, em Mha
#' (1 Mha = 10 mil km2 = 1 milhão de hectares). Mesmo estilo de map_flux_annual(): tile `cartolight`,
#' contorno da América do Sul, biomas, legenda dentro, seta norte e escala. Figuras em inglês, sem título.
#'
#' @param df         data.frame com lon, lat, year e a coluna `value` em km2 (por exemplo,
#'                   `fire$annual` com `burned_km2` ou `prodes$cell` com `def_km2`). Linhas repetidas
#'                   da mesma célula x ano (camadas de biomas diferentes) são somadas.
#' @param value      Coluna em km2.
#' @param label      Nome da grandeza na legenda (ex.: "Burned area", "Deforestation").
#' @param years      Anos a mostrar. NULL = todos.
#' @param cells      data.frame com lon e lat das células a desenhar (por exemplo `fire$land`). As que não
#'                   aparecem em `df` num ano são desenhadas com 0 (em vez de ficarem sem cor). NULL = só
#'                   as células presentes em `df`.
#' @param coverage   data.frame com biome e year dos biomas x anos que têm dado no produto (por exemplo
#'                   `prodes$total |> dplyr::filter(n_poly > 0)`). Os demais ficam sem cor (não são 0).
#'                   Usa os nomes de bioma de `biome_col` em `biomes`; precisa de `cells` e `biomes`.
#' @param biomes     sf com os biomas (geobr::read_biomes(year = 2019)). NULL = sem biomas.
#' @param biome_col  Coluna com o nome do bioma.
#' @param limits    Limites da escala em Mha para os painéis anuais, c(min, max). NULL = 0 até o
#'                   percentil `q_max` dos valores (a mesma escala para todos os anos).
#' @param q_max      Percentil usado como máximo da escala quando `limits` é NULL.
#' @param palette    Opção viridis ("inferno", "magma", "viridis"...).
#' @param total      Também gera o mapa do acumulado do período (escala própria).
#' @param ncol       Colunas de painéis.
#' @param res        Resolução da grade em graus (1 ou 2 valores). NULL = inferida.
#' @param south_america,map_theme,basemap,exclude_biomes,xlim,ylim,decorations,biome_linewidth,biome_colour
#'                   Como em map_flux_annual().
#'
#' @return lista: maps (annual, total), data (célula x ano com a coluna `Mha`), limits (em Mha).
map_disturbance_annual <- function(df,
                                   value          = "burned_km2",
                                   label          = "Burned area",
                                   years          = NULL,
                                   cells          = NULL,
                                   coverage       = NULL,
                                   biomes         = NULL,
                                   biome_col      = "name_biome",
                                   limits         = NULL,
                                   q_max          = 0.99,
                                   palette        = "inferno",
                                   total          = TRUE,
                                   ncol           = 3,
                                   res            = NULL,
                                   south_america  = NULL,
                                   map_theme      = NULL,
                                   basemap        = TRUE,
                                   exclude_biomes = "Sistema Costeiro",
                                   xlim           = c(-75, -35),
                                   ylim           = c(-35, 5.5),
                                   decorations    = TRUE,
                                   biome_linewidth = 0.5,
                                   biome_colour   = "black") {
  need <- c("lon", "lat", "year", value)
  if (!all(need %in% names(df))) stop("df precisa de: ", paste(setdiff(need, names(df)), collapse = ", "))
  d <- as.data.frame(df)
  d$v <- d[[value]]
  d <- d |>
    dplyr::group_by(lon = round(lon, 4), lat = round(lat, 4), year) |>
    dplyr::summarise(Mha = sum(v, na.rm = TRUE) / 1e4, .groups = "drop")
  if (!is.null(years)) d <- d[d$year %in% years, ]
  if (nrow(d) == 0) stop("Nenhuma linha depois do filtro de anos.")

  # células sem registro viram zero (se `cells` for dado); com `coverage`, os biomas x anos sem
  # dado no produto ficam NA (não desenhados) em vez de zero
  if (!is.null(cells)) {
    yrs <- sort(unique(d$year))
    g <- merge(data.frame(lon = round(cells$lon, 4), lat = round(cells$lat, 4)), data.frame(year = yrs))
    d <- merge(g, d, by = c("lon", "lat", "year"), all.x = TRUE)
    d$Mha[is.na(d$Mha)] <- 0
    if (!is.null(coverage)) {
      if (is.null(biomes)) stop("`coverage` precisa de `biomes` para localizar o bioma de cada célula.")
      s2_old <- sf::sf_use_s2(FALSE); on.exit(sf::sf_use_s2(s2_old), add = TRUE)
      bb  <- sf::st_make_valid(sf::st_transform(biomes, 4326))
      xy  <- unique(d[, c("lon", "lat")])
      hit <- sf::st_intersects(sf::st_as_sf(xy, coords = c("lon", "lat"), crs = 4326), bb)
      xy$biome <- as.character(sf::st_drop_geometry(bb)[[biome_col]])[
        vapply(hit, function(i) if (length(i)) i[1] else NA_integer_, integer(1))]
      d <- merge(d, xy, by = c("lon", "lat"))
      ok <- paste(d$biome, d$year) %in% paste(coverage$biome, coverage$year)
      d$Mha[!ok] <- NA_real_
      d <- d[!is.na(d$Mha), setdiff(names(d), "biome")]
    }
  }
  tot <- d |> dplyr::group_by(lon, lat) |> dplyr::summarise(Mha = sum(Mha), .groups = "drop")

  has_gsp <- requireNamespace("ggspatial", quietly = TRUE)
  step <- function(v) { v <- sort(unique(round(v, 4))); if (length(v) < 2) NA_real_ else min(diff(v)) }
  if (is.null(res)) res <- c(step(d$lon), step(d$lat))
  if (anyNA(res)) res <- c(1, 1)
  if (length(res) == 1) res <- rep(res, 2)

  if (is.null(south_america) && requireNamespace("rnaturalearth", quietly = TRUE)) {
    south_america <- try(rnaturalearth::ne_countries(continent = "South America", returnclass = "sf"), silent = TRUE)
    if (inherits(south_america, "try-error")) south_america <- NULL
  }
  b_draw <- NULL
  if (!is.null(biomes)) {
    b <- sf::st_make_valid(sf::st_transform(biomes, 4326))
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
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(1, 0), legend.justification = c(1, 0))
  else ggplot2::theme(legend.position = c(1, 0), legend.justification = c(1, 0))
  th <- c(th_user, list(
    ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 20),
                   text = ggplot2::element_text(size = 20),
                   legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.8), colour = NA)),
    legend_in))
  bar <- ggplot2::guide_colourbar(barheight = grid::unit(4, "cm"))
  leg <- paste0(gsub(" ", "\n", label), "\n(Mha)")

  build <- function(data, lim, facet) {
    p <- ggplot2::ggplot(data, ggplot2::aes(lon, lat))
    for (l in base_layers) p <- p + l
    p <- p + ggplot2::geom_tile(ggplot2::aes(color = Mha, fill = Mha), width = res[1], height = res[2]) +
      ggplot2::scale_color_viridis_c(option = palette, limits = lim, oob = scales::squish, name = leg, guide = bar) +
      ggplot2::scale_fill_viridis_c(option = palette, limits = lim, oob = scales::squish, name = leg, guide = bar)
    for (l in top_border) p <- p + l
    if (facet) p <- p + ggplot2::facet_wrap(~ year, ncol = ncol)
    for (l in deco) p <- p + l
    p <- p + ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                               crs = sf::st_crs(4326), default_crs = sf::st_crs(4326)) +
      ggplot2::scale_x_continuous(breaks = seq(-70, -40, by = 10)) +
      ggplot2::scale_y_continuous(breaks = seq(-30, 0, by = 10)) +
      ggplot2::labs(x = "Longitude", y = "Latitude") + ggplot2::theme_bw()
    for (l in th) if (!is.null(l)) p <- p + l
    p
  }
  q <- function(v) as.numeric(stats::quantile(v[v > 0], q_max, na.rm = TRUE))
  lim_a <- if (!is.null(limits)) limits else c(0, q(d$Mha))
  maps <- list(annual = build(d, lim_a, TRUE))
  lim_t <- c(0, q(tot$Mha))
  if (total) maps$total <- build(tot, lim_t, FALSE)
  list(maps = maps, data = tibble::as_tibble(d), limits = list(annual = lim_a, total = lim_t))
}

#' Salva os mapas de map_disturbance_annual() em PNG: <prefix>_annual.png e <prefix>_total.png
save_disturbance_maps <- function(res, dir = ".", prefix = "disturbance", dpi = 300,
                                  panel = 6, ncol = 3) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  n <- length(unique(res$data$year)); nr <- ceiling(n / ncol)
  ggplot2::ggsave(file.path(dir, paste0(prefix, "_annual.png")), res$maps$annual,
                  width = panel * ncol, height = panel * nr, dpi = dpi, bg = "white", limitsize = FALSE)
  if (!is.null(res$maps$total))
    ggplot2::ggsave(file.path(dir, paste0(prefix, "_total.png")), res$maps$total,
                    width = panel * 1.4, height = panel * 1.2, dpi = dpi, bg = "white")
  invisible(dir)
}

