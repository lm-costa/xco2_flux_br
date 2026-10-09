#' Extrai fluxos do OCO-2 MIP (gridded fluxes, .nc4) e salva um .rds por arquivo
#'
#' Lê com ncdf4 os arquivos "<Submission>_gridded_fluxes_<Experimento>.nc4"
#' (variáveis [longitude, latitude, n_months]), faz recorte espacial opcional,
#' converte o tempo para data/hora e salva em formato longo (uma linha por
#' célula x mês) em .rds.
#'
#' @param in_dir       Pasta com os .nc4 de fluxo.
#' @param out_dir      Pasta de saída (criada se não existir).
#' @param geometry     sf/sfc ou terra::SpatVector (qualquer CRS). NULL = sem recorte.
#'                     É validada (st_make_valid), dissolvida (st_union) e recebe
#'                     buffer; mantém as células cujo centro cai dentro dela.
#' @param buffer       Buffer em graus (padrão 0.05). 0 ou NULL desativa.
#' @param submissions  Nomes das submissions (ex.: c("CT", "EnsMean")). NULL = todas
#'                     (inclui EnsMean e EnsStd, se existirem na pasta).
#' @param experiments  Nomes dos experimentos (ex.: c("IS", "LNLGIS", "Prior")).
#'                     NULL = todos.
#' @param variables    Variáveis de fluxo a extrair.
#' @param drop_na      Remove linhas em que todas as variáveis extraídas são NA.
#' @param overwrite    Sobrescreve .rds já existentes?
#'
#' @return (invisível) vetor com os caminhos dos arquivos salvos.
flux_extractor <- function(in_dir,
                           out_dir,
                           geometry    = NULL,
                           buffer      = 0.05,
                           submissions = NULL,
                           experiments = NULL,
                           variables   = c("land_flux", "ocean_flux",
                                           "fossil_flux", "net_flux"),
                           drop_na     = TRUE,
                           overwrite   = TRUE) {

  stopifnot(dir.exists(in_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  clean <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

  # ---- inventário: <Submission>_gridded_fluxes_<Experimento>.nc4 -----------
  pat   <- "^(.*)_gridded_fluxes_(.*)\\.nc4?$"
  files <- list.files(in_dir, pattern = pat, full.names = TRUE)
  if (length(files) == 0) stop("Nenhum arquivo '*_gridded_fluxes_*.nc4' em ", in_dir)
  bn    <- basename(files)
  sub_f <- sub(pat, "\\1", bn)
  exp_f <- sub(pat, "\\2", bn)

  pick <- function(sel, all, label) {
    if (is.null(sel)) return(rep(TRUE, length(all)))
    miss <- setdiff(sel, all)
    if (length(miss) > 0)
      stop(label, " não encontrado(s): ", paste(miss, collapse = ", "),
           "\nDisponíveis: ", paste(sort(unique(all)), collapse = ", "))
    all %in% sel
  }
  use   <- pick(submissions, sub_f, "submissions") &
    pick(experiments, exp_f, "experiments")
  files <- files[use]; sub_f <- sub_f[use]; exp_f <- exp_f[use]
  if (length(files) == 0) stop("Nenhum arquivo para a combinação submissions x experiments pedida.")

  # ---- geometria: make_valid -> CRS 4326 -> dissolve -> buffer -> make_valid
  geom <- NULL
  if (!is.null(geometry)) {
    if (inherits(geometry, "SpatVector")) geometry <- sf::st_as_sf(geometry)
    s2_old <- sf::sf_use_s2(FALSE)   # buffer em graus exige s2 desligado
    on.exit(sf::sf_use_s2(s2_old), add = TRUE)

    geom <- sf::st_geometry(geometry)
    geom <- sf::st_make_valid(geom)
    geom <- sf::st_transform(geom, 4326)
    geom <- sf::st_union(geom)
    if (!is.null(buffer) && buffer > 0) geom <- sf::st_buffer(geom, dist = buffer)
    geom <- sf::st_make_valid(geom)
    bb   <- sf::st_bbox(geom)
  }

  saved <- character(0)

  for (k in seq_along(files)) {

    out_file <- file.path(out_dir, paste0(clean(sub_f[k]), "__", clean(exp_f[k]), ".rds"))
    if (file.exists(out_file) && !overwrite) {
      message("Já existe, pulando: ", basename(out_file)); next
    }

    nc <- ncdf4::nc_open(files[k])

    # ---- grade --------------------------------------------------------------
    lon <- as.numeric(ncdf4::ncvar_get(nc, "longitude"))
    lat <- as.numeric(ncdf4::ncvar_get(nc, "latitude"))
    lon[lon > 180] <- lon[lon > 180] - 360       # garante -180..180
    nlon <- length(lon); nlat <- length(lat)

    # ---- máscara de células dentro da geometria ----------------------------
    mask <- matrix(TRUE, nlon, nlat)
    if (!is.null(geom)) {
      mask[] <- FALSE
      ii <- which(lon >= bb["xmin"] & lon <= bb["xmax"])
      jj <- which(lat >= bb["ymin"] & lat <= bb["ymax"])
      if (length(ii) > 0 && length(jj) > 0) {
        g   <- expand.grid(i = ii, j = jj)
        pts <- sf::st_as_sf(data.frame(lon = lon[g$i], lat = lat[g$j]),
                            coords = c("lon", "lat"), crs = 4326)
        inside <- lengths(sf::st_intersects(pts, geom)) > 0
        mask[cbind(g$i[inside], g$j[inside])] <- TRUE
      }
    }
    if (!any(mask)) {
      ncdf4::nc_close(nc)
      warning("Nenhuma célula dentro do recorte em ", basename(files[k]), "; não salvo.")
      next
    }

    # janela mínima que contém o recorte (lê só esse bloco do arquivo)
    ix <- range(which(rowSums(mask) > 0))
    iy <- range(which(colSums(mask) > 0))
    nx <- diff(ix) + 1; ny <- diff(iy) + 1
    m_sub <- mask[ix[1]:ix[2], iy[1]:iy[2], drop = FALSE]
    lon_s <- lon[ix[1]:ix[2]]
    lat_s <- lat[iy[1]:iy[2]]

    # ---- tempo --------------------------------------------------------------
    time <- as.numeric(ncdf4::ncvar_get(nc, "time"))
    nt   <- length(time)
    datetime <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")

    # índices (ordem de armazenamento: lon varia mais rápido, depois lat, depois tempo)
    ncell    <- nx * ny
    keep_all <- which(rep(as.vector(m_sub), times = nt))
    cell     <- rep(seq_len(ncell), times = nt)[keep_all]
    t_idx    <- rep(seq_len(nt), each = ncell)[keep_all]
    lon_cell <- rep(lon_s, times = ny)
    lat_cell <- rep(lat_s, each  = nx)

    # ---- variáveis de fluxo -------------------------------------------------
    vars <- intersect(variables, names(nc$var))
    if (length(vars) == 0) {
      ncdf4::nc_close(nc)
      stop("Nenhuma das variáveis pedidas existe em ", basename(files[k]),
           ". Disponíveis: ", paste(names(nc$var), collapse = ", "))
    }

    out <- data.frame(
      lon        = lon_cell[cell],
      lat        = lat_cell[cell],
      time       = time[t_idx],
      datetime   = datetime[t_idx],
      year       = as.integer(format(datetime[t_idx], "%Y")),
      month      = as.integer(format(datetime[t_idx], "%m")),
      submission = sub_f[k],
      experiment = exp_f[k]
    )
    for (v in vars) {
      a <- ncdf4::ncvar_get(nc, v,
                            start = c(ix[1], iy[1], 1),
                            count = c(nx, ny, -1))      # fill values -> NA
      out[[v]] <- as.numeric(a)[keep_all]
    }
    ncdf4::nc_close(nc)

    out <- tibble::as_tibble(out)
    if (drop_na) out <- out[rowSums(!is.na(out[vars])) > 0, ]

    if (nrow(out) == 0) {
      warning("Sem dados válidos em ", basename(files[k]), "; não salvo.")
      next
    }

    saveRDS(out, out_file)
    saved <- c(saved, out_file)
    message(sprintf("Salvo: %s (%s linhas)", basename(out_file),
                    format(nrow(out), big.mark = ".")))
  }

  invisible(saved)
}

# ---------------------------------------------------------------------------
# Exemplos de uso
# ---------------------------------------------------------------------------
# pasta <- list.dirs("data-raw/")[3]
#
# # Só o ensemble médio do experimento LNLGIS, recortado em SP:
# flux_extractor(
#   in_dir      = pasta,
#   out_dir     = "teste_flux/",
#   geometry    = geobr::read_state(year = 2025, code_state = "SP"),
#   submissions = "EnsMean",
#   experiments = "LNLGIS"
# )
#
# # Todas as submissions e experimentos, só net_flux:
# flux_extractor(pasta, "teste_flux/", geometry = geobr::read_state(code_state = "SP"),
#                variables = "net_flux")
