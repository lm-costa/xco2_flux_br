#' Carrega os pacotes do pipeline e instala os que faltarem
#'
#' @param extra     Pacotes adicionais (vetor de nomes) a carregar/instalar junto.
#' @param groups    Grupos de pacotes: "core" (agregação, comparação e gráficos),
#'                  "maps" (mapas e biomas), "io" (leitura de NetCDF/Excel, extratores).
#'                  Padrão: todos.
#' @param install   Se TRUE (padrão), instala o que faltar; se FALSE, só avisa.
#' @param load      Se TRUE (padrão), carrega (library) os pacotes; se FALSE, só garante
#'                  que estão instalados. Os scripts usam pkg::funcao(), então carregar
#'                  não é obrigatório, mas as funções do dplyr/ggplot2 no seu código
#'                  (mutate, ggplot...) precisam de library().
#' @param repos     Repositório CRAN.
#' @param github    Pacotes do GitHub, como c("owner/repo", ...) (usa remotes).
#' @param quiet     Se TRUE, não imprime mensagens de progresso.
#'
#' @return (invisível) data.frame com package, group, installed, loaded, version.
setup_packages <- function(extra   = NULL,
                           groups  = c("core", "maps", "io"),
                           install = TRUE,
                           load    = TRUE,
                           repos   = getOption("repos")[["CRAN"]],
                           github  = NULL,
                           quiet   = FALSE) {
  groups <- match.arg(groups, c("core", "maps", "io"), several.ok = TRUE)
  if (is.null(repos) || is.na(repos) || repos %in% c("@CRAN@", ""))
    repos <- "https://cloud.r-project.org"

  pkgs <- list(
    core = c("dplyr", "tibble", "ggplot2", "scales", "tidyr", "lubridate"),
    maps = c("sf", "geobr", "ggspatial", "rnaturalearth", "rnaturalearthdata"),
    io   = c("ncdf4", "readxl", "writexl")
  )
  tab <- do.call(rbind, lapply(groups, function(g)
    data.frame(package = pkgs[[g]], group = g, stringsAsFactors = FALSE)))
  if (length(extra))
    tab <- rbind(tab, data.frame(package = extra, group = "extra", stringsAsFactors = FALSE))
  tab <- tab[!duplicated(tab$package), ]
  msg <- function(...) if (!quiet) message(...)

  have <- function(p) requireNamespace(p, quietly = TRUE)
  missing <- tab$package[!vapply(tab$package, have, logical(1))]

  if (length(missing)) {
    if (!install) {
      warning("Pacotes ausentes (install = FALSE): ", paste(missing, collapse = ", "), call. = FALSE)
    } else {
      msg("Instalando: ", paste(missing, collapse = ", "))
      for (p in missing) {
        ok <- tryCatch({
          utils::install.packages(p, repos = repos, quiet = quiet)
          have(p)
        }, error = function(e) FALSE, warning = function(w) have(p))
        if (!ok) warning("Não consegui instalar '", p, "'. Instale manualmente ",
                         "(sf/geobr podem precisar de bibliotecas do sistema: GDAL, GEOS, PROJ).",
                         call. = FALSE)
      }
    }
  }

  if (length(github)) {
    if (!have("remotes") && install) utils::install.packages("remotes", repos = repos, quiet = quiet)
    for (r in github) {
      nm <- sub(".*/", "", sub("@.*$", "", r))
      if (!have(nm)) {
        if (install && have("remotes")) {
          msg("Instalando do GitHub: ", r)
          tryCatch(remotes::install_github(r, quiet = quiet, upgrade = "never"),
                   error = function(e) warning("Falha em ", r, ": ", conditionMessage(e), call. = FALSE))
        } else warning("Pacote do GitHub ausente: ", r, call. = FALSE)
      }
      tab <- rbind(tab, data.frame(package = nm, group = "github", stringsAsFactors = FALSE))
    }
  }

  tab$installed <- vapply(tab$package, have, logical(1))
  tab$loaded <- FALSE
  if (load) {
    for (i in which(tab$installed)) {
      tab$loaded[i] <- suppressPackageStartupMessages(
        tryCatch(require(tab$package[i], character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE),
                 error = function(e) FALSE))
    }
  }
  tab$version <- vapply(seq_len(nrow(tab)), function(i)
    if (tab$installed[i]) as.character(utils::packageVersion(tab$package[i])) else NA_character_, "")

  bad <- tab$package[!tab$installed]
  if (length(bad)) msg("Ainda ausentes: ", paste(bad, collapse = ", "))
  else msg("Pacotes prontos: ", sum(tab$installed), ".")
  invisible(tab)
}

# Exemplo:
# source("setup_packages.R")
# setup_packages()                              # tudo
# setup_packages(groups = c("core", "maps"))    # sem leitura de NetCDF/Excel
# setup_packages(extra = c("purrr", "readr"))
# setup_packages(install = FALSE)               # só checa
