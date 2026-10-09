#' Ensemble das submissions do OCO-2 MIP, por experimento
#'
#' Lê os .rds gerados por mip_extractor() (padrão "Submission__Experimento.rds"),
#' agrupa as submissions de cada experimento pelas coordenadas/tempo da
#' observação (lon, lat, time) e calcula estatísticas do ensemble sobre o
#' xco2 simulado (xco2_sim). Salva um .rds por experimento.
#'
#' @param in_dir       Pasta onde estão os .rds do mip_extractor().
#' @param experiments  Vetor com os experimentos desejados (ex.: c("IS", "LNLG"))
#'                     ou "ALL" para todos, exceto "Prior".
#' @param out_dir      Pasta de saída. Padrão: subpasta "ensemble" dentro de in_dir.
#' @param min_sub      Nº mínimo de submissions com valor válido para manter a
#'                     observação no ensemble (padrão 1).
#' @param overwrite    Sobrescreve arquivos existentes?
#'
#' @return (invisível) vetor com os caminhos dos arquivos salvos.
mip_ensemble <- function(in_dir,
                         experiments = "ALL",
                         out_dir     = file.path(in_dir, "ensemble"),
                         min_sub     = 1,
                         overwrite   = TRUE) {

  stopifnot(dir.exists(in_dir))
  clean <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

  # ---- inventário dos arquivos: Submission__Experimento.rds ----------------
  files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
  nm    <- tools::file_path_sans_ext(basename(files))
  ok    <- grepl("__", nm, fixed = TRUE)
  files <- files[ok]
  nm    <- nm[ok]
  if (length(files) == 0) stop("Nenhum arquivo 'Submission__Experimento.rds' em ", in_dir)

  exp_f <- sub("^(.*?)__(.*)$", "\\2", nm, perl = TRUE)

  # ---- seleção dos experimentos --------------------------------------------
  avail <- sort(unique(exp_f))
  if (length(experiments) == 1 && toupper(experiments) == "ALL") {
    sel <- avail[!grepl("^prior$", avail, ignore.case = TRUE)]
  } else {
    sel  <- clean(trimws(experiments))
    miss <- setdiff(sel, avail)
    if (length(miss) > 0)
      stop("Experimento(s) não encontrado(s): ", paste(miss, collapse = ", "),
           "\nDisponíveis: ", paste(avail, collapse = ", "))
  }
  if (length(sel) == 0) stop("Nenhum experimento selecionado.")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # colunas que são iguais em todas as submissions (vêm da base observada)
  static_cols <- c("datetime", "xco2", "uncertanty", "model_error", "altitude",
                   "psurf", "airmass", "tcwv", "assimilate_flag")
  saved <- character(0)

  for (e in sel) {

    out_file <- file.path(out_dir, paste0("ensemble__", e, ".rds"))
    if (file.exists(out_file) && !overwrite) {
      message("Já existe, pulando: ", basename(out_file))
      next
    }

    f <- files[exp_f == e]
    message(sprintf("Experimento %s: %d submission(s)", e, length(f)))

    dat <- dplyr::bind_rows(lapply(f, readRDS))
    keep_static <- intersect(static_cols, names(dat))

    ens <- dat |>
      dplyr::group_by(lon, lat, time) |>
      dplyr::summarise(
        experiment      = dplyr::first(experiment),
        n_sub           = sum(!is.na(xco2_sim)),
        xco2_ens_mean   = mean(xco2_sim,   na.rm = TRUE),
        xco2_ens_median = stats::median(xco2_sim, na.rm = TRUE),
        xco2_ens_se     = stats::sd(xco2_sim,     na.rm = TRUE)/sqrt(n_sub),
        xco2_ens_min    = min(xco2_sim,    na.rm = TRUE),
        xco2_ens_max    = max(xco2_sim,    na.rm = TRUE),
        dplyr::across(dplyr::all_of(keep_static), dplyr::first),
        .groups = "drop"
      ) |>
      dplyr::filter(n_sub >= min_sub) |>
      dplyr::relocate(lon, lat, time, dplyr::any_of("datetime"), experiment, n_sub)

    rm(dat); invisible(gc())

    if (nrow(ens) == 0) {
      warning("Sem observações válidas para o experimento ", e, "; não salvo.")
      next
    }

    saveRDS(ens, out_file)
    saved <- c(saved, out_file)
    message(sprintf("Salvo: %s (%s linhas)", basename(out_file),
                    format(nrow(ens), big.mark = ".")))
  }

  invisible(saved)
}

# ---------------------------------------------------------------------------
# Exemplos de uso
# ---------------------------------------------------------------------------
# # Todos os experimentos, exceto Prior:
# mip_ensemble(in_dir = "teste/", experiments = "ALL")
#
# # Combinação específica de experimentos (um ensemble por experimento):
# mip_ensemble("teste/", experiments = c("IS", "LNLG"), out_dir = "teste/ens")
#
# # Exigir pelo menos 5 submissions por observação:
# mip_ensemble("teste/", "ALL", min_sub = 5)
