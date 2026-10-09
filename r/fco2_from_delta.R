#' Converte o delta de XCO2 (ppm) em fluxo por sondagem
#'
#' Hipótese: o excesso de CO2 na coluna (delta = xco2_obs - xco2_sim_IS) se
#' acumulou ao longo de um tempo de residência `tau_days` sobre a sondagem.
#'
#'   n_dry  = (psurf*100/g - tcwv) / M_dry          [mol ar seco m-2]
#'   dM     = delta * 1e-6 * n_dry * M_C            [g C m-2]
#'   fco2   = dM / tau_days                         [g C m-2 d-1]
#'
#' @param df         data.frame com as colunas de delta, psurf e tcwv.
#' @param delta      Nome da coluna do delta (ppm). Padrão "delta".
#' @param psurf      Nome da coluna de pressão de superfície (hPa).
#' @param tcwv       Nome da coluna de vapor d'água total da coluna (kg m-2).
#' @param sigma_delta Nome da coluna de incerteza do delta (ppm), ex.: "total_erro".
#'                   NULL = não propaga incerteza.
#' @param tau_days   Tempo de residência em dias (padrão 1).
#' @param sigma_tau_days Incerteza (1 desvio padrão) do tau, em dias. NULL = tau exato.
#' @param method     "massa" (balanço de massa com psurf/tcwv), "SPT" ou "ambos".
#' @param spt_height,spt_mw,spt_vmol,spt_scale Constantes do SPT:
#'                   spt = spt_height * delta * (spt_mw/spt_vmol) * spt_scale
#'                   (padrão: 10000 * delta * (44/25) * 1e-3). Gera `spt`
#'                   (e `sigma_spt` se sigma_delta for informado).
#'
#' @return df com as colunas novas: k_gC_ppm (g C m-2 por ppm), dM_gC_m2,
#'         fco2 (g C m-2 d-1) e, se sigma_delta informado, sigma_f.
fco2_from_delta <- function(df,
                            delta          = "delta",
                            psurf          = "psurf",
                            tcwv           = "tcwv",
                            sigma_delta    = NULL,
                            tau_days       = 1,
                            sigma_tau_days = NULL,
                            method         = c("ambos", "massa", "SPT"),
                            spt_height     = 10000,
                            spt_mw         = 44,
                            spt_vmol       = 25,
                            spt_scale      = 1e-3) {

  method <- match.arg(method)

  need <- c(delta, if (method != "SPT") c(psurf, tcwv), sigma_delta)
  if (!all(need %in% names(df)))
    stop("Colunas ausentes em df: ", paste(setdiff(need, names(df)), collapse = ", "))

  if (method %in% c("massa", "ambos")) {

    g     <- 9.80665    # m s-2
    M_dry <- 0.028964   # kg mol-1
    M_C   <- 12.011     # g C mol-1

    # massa de ar seco por m2: massa total (psurf/g) menos a massa de água (tcwv)
    n_dry <- (df[[psurf]] * 100 / g - df[[tcwv]]) / M_dry      # mol m-2

    # g C m-2 por ppm de delta
    df$k_gC_ppm <- 1e-6 * n_dry * M_C

    df$dM_gC_m2 <- df[[delta]] * df$k_gC_ppm                   # g C m-2
    df$fco2     <- df$dM_gC_m2 / tau_days                      # g C m-2 d-1

    if (!is.null(sigma_delta)) {
      # σ_f² = (k σ_delta / tau)² + (fco2 σ_tau / tau)²
      var_f <- (df$k_gC_ppm * df[[sigma_delta]] / tau_days)^2
      if (!is.null(sigma_tau_days)) {
        var_f <- var_f + (df$fco2 * sigma_tau_days / tau_days)^2
      }
      df$sigma_f <- sqrt(var_f)
    }
  }   # fim do método "massa"

  if (method %in% c("SPT", "ambos")) {
    # SPT = 10000 * delta * (44/25) * 1e-3
    #   delta (ppm) * 44/25 -> mg CO2 m-3 (44 g/mol; ~25 L/mol de ar)
    #   * 1e-3              -> g CO2 m-3
    #   * 10000             -> g CO2 m-2 (camada de spt_height m, densidade constante)
    # Não tem base de tempo: é massa por área (g CO2 m-2).
    fator    <- spt_height * (spt_mw / spt_vmol) * spt_scale
    df$spt   <- df[[delta]] * fator
    if (!is.null(sigma_delta)) df$sigma_spt <- df[[sigma_delta]] * fator
  }

  df
}

# ---------------------------------------------------------------------------
# Exemplo de uso
# ---------------------------------------------------------------------------
# df_n <- fco2_from_delta(df_n, sigma_delta = "total_erro", tau_days = 1)
#
# # Sensibilidade ao tau:
# for (tau in c(0.5, 1, 2, 3)) {
#   tmp <- fco2_from_delta(df_n, sigma_delta = "total_erro", tau_days = tau)
#   cat("tau =", tau, "d | fco2 médio =", mean(tmp$fco2, na.rm = TRUE), "\n")
# }
#
# # Depois, a agregação na célula do fluxo:
# base <- agg_xco2_flux(df_n, df_flux, vars = "fco2", sigma = "sigma_f",
#                       min_n = 5, min_days = 2)
