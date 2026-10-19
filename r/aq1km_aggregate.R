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
