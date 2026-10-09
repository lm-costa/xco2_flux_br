#' Extrai dados do OCO-2 MIP (NetCDF) e salva um .rds por submission x experiment
#'
#' @param nc_file       Caminho do arquivo .nc de entrada.
#' @param out_dir       Pasta de saída (criada se não existir).
#' @param geometry      Geometria para o recorte espacial: objeto sf/sfc ou
#'                      terra::SpatVector (qualquer CRS; é convertido para
#'                      EPSG:4326). Se NULL, não faz recorte espacial.
#'                      A geometria é validada (st_make_valid), dissolvida
#'                      (st_union) e recebe buffer antes do recorte.
#' @param buffer        Buffer em graus aplicado à geometria (padrão 0.05).
#'                      Use 0 ou NULL para desativar.
#' @param submissions   Vetor com nomes (ou índices) das submissions desejadas.
#'                      NULL = todas.
#' @param experiments   Vetor com nomes (ou índices) dos experimentos desejados.
#'                      NULL = todos.
#' @param drop_na       Remove linhas com xco2 simulado NA (fill values).
#' @param overwrite     Sobrescreve .rds já existentes?
#'
#' @return (invisível) vetor com os caminhos dos arquivos salvos.
mip_extractor <- function(nc_file,
                          out_dir,
                          geometry        = NULL,
                          buffer          = 0.05,
                          submissions     = NULL,
                          experiments     = NULL,
                          drop_na         = TRUE,
                          overwrite       = TRUE) {

  stopifnot(file.exists(nc_file))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  nc <- ncdf4::nc_open(nc_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  # ---- nomes de submissions / experimentos ---------------------------------
  all_subs <- trimws(as.character(ncdf4::ncvar_get(nc, "submissions")))
  all_exps <- trimws(as.character(ncdf4::ncvar_get(nc, "experiments")))

  pick <- function(sel, all, label) {
    if (is.null(sel)) return(seq_along(all))
    if (is.numeric(sel)) {
      if (any(sel < 1 | sel > length(all)))
        stop(label, ": índice fora do intervalo 1-", length(all))
      return(as.integer(sel))
    }
    idx <- match(trimws(sel), all)
    if (anyNA(idx))
      stop(label, " não encontrado(s): ", paste(sel[is.na(idx)], collapse = ", "),
           "\nDisponíveis: ", paste(all, collapse = ", "))
    idx
  }
  i_subs <- pick(submissions, all_subs, "submissions")
  i_exps <- pick(experiments, all_exps, "experiments")

  # ---- recorte espacial (feito uma única vez) ------------------------------
  lon <- as.numeric(ncdf4::ncvar_get(nc, "longitude"))
  lat <- as.numeric(ncdf4::ncvar_get(nc, "latitude"))
  keep <- !is.na(lon) & !is.na(lat)

  if (!is.null(geometry)) {
    if (inherits(geometry, "SpatVector")) geometry <- sf::st_as_sf(geometry)

    # buffer em graus exige s2 desligado (com s2 ligado a distância seria em metros)
    s2_old <- sf::sf_use_s2(FALSE)
    on.exit(sf::sf_use_s2(s2_old), add = TRUE)

    geom <- sf::st_geometry(geometry)
    geom <- sf::st_make_valid(geom)            # corrige geometrias inválidas
    geom <- sf::st_transform(geom, 4326)       # CRS dos dados (graus)
    geom <- sf::st_union(geom)                 # dissolve (multipolygon / várias feições)
    if (!is.null(buffer) && buffer > 0) {
      geom <- sf::st_buffer(geom, dist = buffer)  # buffer em graus
    }
    geom <- sf::st_make_valid(geom)            # garante validade após o buffer

    # 1) pré-filtro rápido por bbox
    bb <- sf::st_bbox(geom)
    keep <- keep &
      lon >= bb["xmin"] & lon <= bb["xmax"] &
      lat >= bb["ymin"] & lat <= bb["ymax"]

    # 2) interseção exata com a geometria apenas nos pontos pré-filtrados
    cand <- which(keep)
    keep[] <- FALSE
    if (length(cand) > 0) {
      pts <- sf::st_as_sf(
        data.frame(lon = lon[cand], lat = lat[cand]),
        coords = c("lon", "lat"), crs = 4326
      )
      inside <- lengths(sf::st_intersects(pts, geom)) > 0
      keep[cand[inside]] <- TRUE
    }
  }


  idx <- which(keep)
  if (length(idx) == 0) stop("Nenhuma observação dentro do recorte/filtros informados.")

  # ---- variáveis comuns (lidas uma vez e já subsetadas) --------------------
  get1 <- function(v) as.numeric(ncdf4::ncvar_get(nc, v))[idx]

  time <- as.numeric(ncdf4::ncvar_get(nc, "time"))[idx]
  base <- data.frame(
    lon         = lon[idx],
    lat         = lat[idx],
    time        = time,
    datetime    = as.POSIXct(time, origin = "1970-01-01", tz = "UTC"),
    xco2        = get1('xco2'),
    uncertanty  = get1("xco2_uncertainty"),
    model_error = get1("model_error"),
    altitude    = get1("altitude"),
    psurf       = get1("psurf"),
    airmass     = get1("airmass"),
    tcwv        = get1("tcwv")
  )

  # ---- loop: um arquivo por experimento x submission -----------------------
  clean <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  saved <- character(0)

  for (ie in i_exps) {
    for (is in i_subs) {

      out_file <- file.path(
        out_dir,
        paste0(clean(all_subs[is]), "__", clean(all_exps[ie]), ".rds")
      )
      if (file.exists(out_file) && !overwrite) {
        message("Já existe, pulando: ", basename(out_file))
        next
      }

      # dim de simulated_values: [experiments, submissions, obs]
      # (chunk [1,1,N] -> leitura de uma combinação por vez, sem estourar RAM)
      sim <- as.numeric(ncdf4::ncvar_get(
        nc, "simulated_values",
        start = c(ie, is, 1), count = c(1, 1, -1)
      ))[idx]

      dft <- tibble::as_tibble(base)
      dft$submission <- all_subs[is]
      dft$experiment <- all_exps[ie]
      dft$xco2_sim   <- sim   # valor simulado pelo modelo (xco2 = observado OCO-2)

      first <- c("lon", "lat", "time", "datetime", "submission",
                 "experiment", "xco2_sim", "xco2")
      dft <- dft[, c(first, setdiff(names(dft), first))]

      if (drop_na) dft <- dft[!is.na(dft$xco2_sim), ]

      if (nrow(dft) == 0) {
        warning("Sem dados para ", all_subs[is], " / ", all_exps[ie], "; não salvo.")
        next
      }

      saveRDS(dft, out_file)
      saved <- c(saved, out_file)
      message(sprintf("Salvo: %s (%s linhas)", basename(out_file),
                      format(nrow(dft), big.mark = ".")))
    }
  }

  invisible(saved)
}

# ---------------------------------------------------------------------------
# Exemplos de uso
# ---------------------------------------------------------------------------
# library(sf)
#
# # Recorte retangular equivalente ao filtro anterior (lon -75/-35, lat -35/5)
# bbox_br <- st_as_sfc(st_bbox(
#   c(xmin = -75, xmax = -35, ymin = -35, ymax = 5), crs = 4326
# ))
#
# # Todas as submissions e experimentos:
# mip_extractor("data/OCO2.nc", out_dir = "output/rds", geometry = bbox_br)
#
# # Selecionando algumas (por nome ou índice) e só dados assimilados:
# mip_extractor(
#   "data/OCO2.nc", "output/rds",
#   geometry        = geobr::read_state(code_state = "SP") ,
#   submissions     = c(1, 3),
#   experiments     = "IS",
#   assimilate_flag = 1
# )
#
# # Ler um resultado:
# dados_1 <- readRDS("output/rds/<submission>__<experiment>.rds")
