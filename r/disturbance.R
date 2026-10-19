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

