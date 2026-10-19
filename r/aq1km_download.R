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
