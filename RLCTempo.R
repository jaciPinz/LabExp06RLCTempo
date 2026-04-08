# ==============================================================================
# Analisi transitorio RLC — tensione ai capi dell'induttore
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(quantities)
  library(nlstools)
  library(cowplot)
  library(gridExtra)
  library(grid)
  library(xtable)
  library(onls)
  library(MASS)
  library(minpack.lm)
  library(rlang)
  library(nortest)
  library(propagate)
  library(BSDA)
  library(deming)
  library(patchwork)
  library(readr)
})


# ==============================================================================
# 1. INCERTEZZE STRUMENTALI
# ==============================================================================

#' Assegna l'incertezza di misura ai dati dell'oscilloscopio.
#'
#' @param data  Vettore numerico o data.frame.
#' @param type  "volt" oppure "time".
#' @param scale Divisione di fondo scala (default 1/25).
#' @param col   Nome della colonna (richiesto se data è un data.frame).
#'
#' @details
#'   - Tensione: sigma = sqrt(sigmaR^2 + sigmaK^2),
#'               con sigmaR = 4% della scala e sigmaK = 1.2% del valore letto.
#'   - Tempo:    sigma = sigmaR costante (solo errore di lettura).
#'
#' @return Lo stesso tipo di `data` con la colonna/vettore di incertezza.
assignUncertainty <- function(data,
                              type,
                              scale = 1 / 25,
                              col = NULL) {
  sigmaR <- 0.04 * scale
  
  if (type == "volt") {
    values  <- if (is.data.frame(data))
      data[[col]]
    else
      data
    sigmaK  <- 0.012 * abs(values)
    uncert  <- sqrt(sigmaR^2 + sigmaK^2)
    
    if (is.data.frame(data)) {
      data[[paste0(col, "Unc")]] <- uncert
      return(data)
    }
    return(uncert)
    
  } else if (type == "time") {
    uncert <- rep(sigmaR, if (is.data.frame(data))
      nrow(data)
      else
        length(data))
    
    if (is.data.frame(data)) {
      data[[paste0(col, "Unc")]] <- uncert
      return(data)
    }
    return(uncert)
    
  } else {
    stop("'type' deve essere \"volt\" o \"time\".")
  }
}


# ==============================================================================
# 2. IMPORTAZIONE E PREPARAZIONE DEI DATI
# ==============================================================================

# Riga da cui iniziano i dati veri (le precedenti sono metadati del CSV)
START_ROW <- 176

raw <- read_csv("Desktop/ALL0000/F0000CH2.CSV", col_names = FALSE)

data <- data.frame(time = as.numeric(raw$X4[START_ROW:nrow(raw)]), volt = as.numeric(raw$X5[START_ROW:nrow(raw)]))
rm(raw)

# Aggiunta delle incertezze (scale ricavate dalla configurazione dell'oscilloscopio)
data <- assignUncertainty(data,
                          type = "time",
                          scale = 0.001,
                          col = "time")
data <- assignUncertainty(data,
                          type = "volt",
                          scale = 0.5,
                          col = "volt")


# ==============================================================================
# 3. MODELLI E FIT NON LINEARE
# ==============================================================================

#' Modello completo di VL(t): include t0 (istante di innesco) e v0 (baseline).
#'
#' Per t <= t0 il segnale vale v0; per t > t0 è un oscillatore smorzato.
modelFull <- function(t, A, omega0, tau, v0, t0) {
  Omega <- sqrt(pmax(1e-12, omega0^2 - 1 / tau^2))
  dt    <- t - t0
  
  Vosc <- A * exp(-dt / tau) *
    (cos(Omega * dt) + sin(Omega * dt) / (Omega * tau))
  
  return(Vosc + v0)
}

#' Modello ridotto di VL(t): usato dopo aver traslato t0->0 e rimosso v0.
modelReduced <- function(t, A, omega0, tau, phi) {
  Omega <- sqrt(pmax(1e-12, omega0^2 - 1 / tau^2))
  A * exp(-t / tau) * (cos(Omega * t + phi) + sin(Omega * t + phi) / (Omega *
                                                                        tau))
}

# ------------------------------------------------------------------------------

#' Stampa i risultati di un fit nls in modo compatto.
#'
#' @param fit  Oggetto nls.
#' @param step 1 = fit completo, 2 = fit ridotto (stampa anche Omega e f_d).
printFitResults <- function(fit, step) {
  params <- coef(fit)
  errs   <- sqrt(diag(vcov(fit)))
  resid  <- residuals(fit)
  w      <- fit$weights %||% rep(1, length(resid))   # %||% è "se NULL usa default"
  
  chi2 <- sum(w * resid^2)
  dof  <- length(resid) - length(params)
  
  cat(sprintf("  %-8s = %+.4e  ±  %.2e\n", names(params), params, errs), sep = "")
  cat(sprintf("  chi2 / dof  =  %.2f / %d  =  %.3f\n", chi2, dof, chi2 /
                dof))
  
  if (step == 2) {
    Omega <- sqrt(max(0, params["omega0"]^2 - 1 / params["tau"]^2))
    cat(sprintf("\n  Omega  = %.4e  rad/s\n", Omega))
    cat(sprintf("  f_d    = %.4e  Hz\n", Omega / (2 * pi)))
    cat(sprintf("  tau    = %.4e  s\n", params["tau"]))
  }
}

# ------------------------------------------------------------------------------

#' Fit del transitorio VL(t) con strategia a due passi.
#'
#' Passo 1 — fit completo: stima t0 e v0 sull'intero segnale.
#' Passo 2 — fit ridotto: taglia i dati a t > t0, rimuove la baseline e
#'           riadatta il modello senza t0/v0 per ottenere stime più precise
#'           di omega0, tau, A, phi.
#'
#' Se fixT0 = TRUE si salta il passo 1 (utile quando t0 = 0 è già noto).
#'
#' @return Lista invisibile con fitFull (solo se !fixT0), fitReduced, dataCut.
fitVL <- function(data,
                  omega0Init = 1,
                  tauInit = 1,
                  AInit = 1,
                  v0Init = 0,
                  t0Init = 0) {
  nlsControl <- nls.control(maxiter = 500, tol = 1e-8)
  
  # --- Fit completo ---
  fitFull <- nlsLM(
    volt ~ modelFull(time, A, omega0, tau, v0, t0),
    data    = data,
    start   = list(
      A = AInit,
      omega0 = omega0Init,
      tau = tauInit,
      v0 = v0Init,
      t0 = t0Init
    ),
    weights = 1 / data$voltUnc^2,
    control = nlsControl
  )
  
  cat("====== Fit completo ======\n")
  printFitResults(fitFull, step = 2)  # puoi lasciare step=2 per stampare Omega
  
  return(invisible(list(fitFull = fitFull)))
}

# ==============================================================================
# 4. STIMA AUTOMATICA DEI PARAMETRI INIZIALI
# ==============================================================================

#' Stima grezza dei parametri iniziali cercando i picchi del segnale.
#'
#' Algoritmo:
#'   1. Individua massimi e minimi locali con un intorno di kNeighbors punti.
#'   2. Filtra i picchi troppo piccoli rispetto al primo (soglia percentuale).
#'   3. Stima tau dal rapporto di ampiezza tra picchi consecutivi,
#'      omega0 dalla frequenza pseudo-oscillatoria e A dal primo picco.
#'
#' @param data              data.frame con colonne "time" e "volt".
#' @param kNeighbors        Ampiezza dell'intorno locale per la ricerca picchi.
#' @param sogliaPercentuale Frazione minima rispetto al primo picco (0–1).
#'
#' @return Lista con aInit, omega0Init, tauInit, phiInit, v0Init, t0Init, indici.
estimateTransientInitialsBrutal <- function(data,
                                            kNeighbors = 25,
                                            sogliaPercentuale = 0.1) {
  tempo    <- as.numeric(data$time)
  tensione <- as.numeric(data$volt)
  n        <- length(tensione)
  
  indiceMassimi <- integer(0)
  indiceMinimi  <- integer(0)
  
  # Ricerca massimi/minimi con gestione plateau
  i <- kNeighbors + 1
  while (i <= n - kNeighbors) {
    j <- i
    while (j < n - kNeighbors &&
           tensione[j + 1] == tensione[i])
      j <- j + 1
    
    val    <- tensione[i]
    sinSx  <- tensione[(i - kNeighbors):(i - 1)]
    sinDx  <- tensione[(j + 1):(j + kNeighbors)]
    idxMid <- round((i + j) / 2)
    
    if (all(val >= sinSx) && all(val >= sinDx)) {
      if (length(indiceMassimi) == 0 ||
          idxMid > max(indiceMassimi) + kNeighbors)
        indiceMassimi <- c(indiceMassimi, idxMid)
    }
    
    if (all(val <= sinSx) && all(val <= sinDx)) {
      if (length(indiceMinimi) == 0 ||
          idxMid > max(indiceMinimi) + kNeighbors)
        indiceMinimi <- c(indiceMinimi, idxMid)
    }
    
    i <- j + 1
  }
  
  # Filtro ampiezza sui massimi
  if (length(indiceMassimi) > 0) {
    soglia <- abs(tensione[indiceMassimi[1]]) * sogliaPercentuale
    indiceMassimi <- indiceMassimi[abs(tensione[indiceMassimi]) > soglia]
  }
  
  # Filtro ampiezza sui minimi
  if (length(indiceMinimi) > 0) {
    soglia <- abs(tensione[indiceMinimi[1]]) * sogliaPercentuale
    indiceMinimi <- indiceMinimi[abs(tensione[indiceMinimi]) > soglia]
  }
  
  # Stima baseline
  if (length(indiceMassimi) >= 1 && length(indiceMinimi) >= 1) {
    stimaV0 <- (tensione[indiceMassimi[1]] + tensione[indiceMinimi[1]]) / 2
  } else {
    stimaV0 <- mean(tensione[1:min(20, n)])
  }
  
  tensioneCentrata <- tensione - stimaV0
  
  # Estrae T, tau, omegaPseudo e omega0 da coppie consecutive
  stimaDaEstremi <- function(tempo, tensioneCentrata, indiceEstremi) {
    indiceEstremi <- sort(indiceEstremi)
    if (length(indiceEstremi) < 2)
      return(NULL)
    
    t <- tempo[indiceEstremi]
    v <- tensioneCentrata[indiceEstremi]
    
    T_i <- diff(t)
    v1  <- abs(v[-length(v)])
    v2  <- abs(v[-1])
    
    ok <- is.finite(T_i) & T_i > 0 &
      is.finite(v1) & is.finite(v2) &
      v1 > 0 & v2 > 0 & v1 > v2
    
    if (!any(ok))
      return(NULL)
    
    T_i <- T_i[ok]
    v1  <- v1[ok]
    v2  <- v2[ok]
    
    omegaPseudo_i <- 2 * pi / T_i
    tau_i <- T_i / log(v1 / v2)
    omega0_i <- sqrt(omegaPseudo_i^2 + 1 / tau_i^2)
    
    list(
      T = T_i,
      tau = tau_i,
      omegaPseudo = omegaPseudo_i,
      omega0 = omega0_i
    )
  }
  
  resMax <- stimaDaEstremi(tempo, tensioneCentrata, indiceMassimi)
  resMin <- stimaDaEstremi(tempo, tensioneCentrata, indiceMinimi)
  
  T_all <- c()
  tau_all <- c()
  omegaPseudo_all <- c()
  omega0_all <- c()
  
  if (!is.null(resMax)) {
    T_all <- c(T_all, resMax$T)
    tau_all <- c(tau_all, resMax$tau)
    omegaPseudo_all <- c(omegaPseudo_all, resMax$omegaPseudo)
    omega0_all <- c(omega0_all, resMax$omega0)
  }
  
  if (!is.null(resMin)) {
    T_all <- c(T_all, resMin$T)
    tau_all <- c(tau_all, resMin$tau)
    omegaPseudo_all <- c(omegaPseudo_all, resMin$omegaPseudo)
    omega0_all <- c(omega0_all, resMin$omega0)
  }
  
  if (length(tau_all) > 0) {
    periodoMedio <- mean(T_all)
    stimaTau     <- mean(tau_all)
    omegaPseudo  <- mean(omegaPseudo_all)
    stimaOmega0  <- mean(omega0_all)
  } else {
    periodoMedio <- NA_real_
    stimaTau     <- (max(tempo) - min(tempo)) / 5
    omegaPseudo  <- NA_real_
    stimaOmega0  <- 1e6
  }
  
  # Ampiezza iniziale: primo estremo valido, altrimenti massimo assoluto
  if (length(indiceMassimi) >= 1) {
    stimaA <- abs(tensioneCentrata[indiceMassimi[1]])
  } else if (length(indiceMinimi) >= 1) {
    stimaA <- abs(tensioneCentrata[indiceMinimi[1]])
  } else {
    stimaA <- max(abs(tensioneCentrata))
  }
  
  # --- Stima robusta di t0 basata sulla derivata iniziale ---
  
  nInit <- min(40, length(tensione))
  
  dVdt <- diff(tensione) / diff(tempo)
  tempoMid <- tempo[-1]
  
  sigmaDer <- sd(dVdt[1:(nInit - 1)], na.rm = TRUE)
  threshold <- 5 * sigmaDer
  
  idx <- which(abs(dVdt) > threshold)
  
  if (length(idx) > 0) {
    stimaT0 <- tempoMid[idx[1]]
  } else {
    stimaT0 <- median(tempo[1:nInit])
  }
  
  cat(sprintf(
    "Picchi validi trovati: max=%d, min=%d\n",
    length(indiceMassimi),
    length(indiceMinimi)
  ))
  
  list(
    aInit      = stimaA,
    omega0Init = stimaOmega0,
    tauInit    = stimaTau,
    phiInit    = 0,
    v0Init     = stimaV0,
    t0Init     = stimaT0,
    periodoMedio = periodoMedio,
    indiciMax   = indiceMassimi,
    indiciMin   = indiceMinimi
  )
}

# ==============================================================================
# 5. ESECUZIONE E VISUALIZZAZIONE
# ==============================================================================

initialGuesses <- estimateTransientInitialsBrutal(data, kNeighbors = 50)

# Subset con solo i punti corrispondenti ai picchi identificati
dataPMax <- data[initialGuesses$indiciMax, ]
dataPMin <- data[initialGuesses$indiciMin, ]

# Plot diagnostico: dati completi + massimi + minimi trovati
plot(
  data$time,
  data$volt,
  col  = "red",
  pch = 16,
  cex = 0.5,
  xlab = "Tempo (s)",
  ylab = "Tensione (V)",
  main = "Segnale VL(t) e picchi identificati"
)
points(
  dataPMax$time,
  dataPMax$volt,
  col = "blue",
  pch = 16,
  cex = 0.8
)
points(
  dataPMin$time,
  dataPMin$volt,
  col = "green",
  pch = 16,
  cex = 0.8
)
legend(
  "topright",
  legend = c("Dati grezzi", "Massimi", "Minimi"),
  col    = c("red", "blue", "green"),
  pch = 16
)

# ==============================================================================
# 6. FIT NON LINEARE
# ==============================================================================

fitResults <- fitVL(
  data       = data,
  AInit      = initialGuesses$aInit,
  omega0Init = initialGuesses$omega0Init,
  tauInit    = initialGuesses$tauInit,
  v0Init     = initialGuesses$v0Init,
  t0Init     = initialGuesses$t0Init
)

# ==============================================================================
# 7. PLOT DIAGNOSTICO CON RESIDUI REALI (VOLT) E BARRE D'ERRORE
# ==============================================================================

# 1. Recupero parametri
params <- coef(fitResults$fitFull)

# 2. Calcolo fitValues e residui REALI (non pesati)
data$fitValues <- modelFull(data$time, params["A"], params["omega0"], params["tau"], params["v0"], params["t0"])

# Calcolo differenza semplice in Volt
data$residReal <- data$volt - data$fitValues

# 3. Griglia per la curva continua
tGrid <- seq(min(data$time), max(data$time), length.out = 2000)
fitDF <- data.frame(time = tGrid,
                    volt = modelFull(tGrid, params["A"], params["omega0"], params["tau"], params["v0"], params["t0"]))

# --- PLOT SEGNALE CON BARRE D'ERRORE ---
p_signal <- ggplot(data, aes(x = time, y = volt)) +
  geom_errorbar(
    aes(ymin = volt - voltUnc, ymax = volt + voltUnc),
    width = 0,
    color = "red",
    alpha = 0.2
  ) +
  geom_point(color = "red",
             size = 0.5,
             alpha = 0.4) +
  geom_line(
    data = fitDF,
    aes(x = time, y = volt),
    color = "black",
    linewidth = 0.8
  ) +
  geom_point(data = data[initialGuesses$indiciMax, ],
             color = "blue",
             size = 2) +
  geom_point(data = data[initialGuesses$indiciMin, ],
             color = "green",
             size = 2) +
  labs(title = "Transitorio RLC — Fit e Incertezze", x = "Tempo (s)", y = "Tensione (V)") +
  theme_minimal()

# --- PLOT RESIDUI REALI (In Volt) ---
p_resid <- ggplot(data, aes(x = time, y = residReal)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black") +
  geom_errorbar(
    aes(ymin = residReal - voltUnc, ymax = residReal + voltUnc),
    width = 0,
    color = "purple",
    alpha = 0.2
  ) +
  geom_point(color = "purple",
             size = 0.7,
             alpha = 0.5) +
  labs(title = "Residui", x = "Tempo (s)", y = "Residui (V)") + # Unità di misura esplicita
  theme_minimal()

(p_signal / p_resid)

# ==============================================================================
# 8. MAPPATURA MANUALE DEL CHI² — INCERTEZZE SUI PARAMETRI
# ==============================================================================
#
# Strategia (analoga al codice Python di riferimento):
#
#   1. Estrai i parametri best fit da fitResults$fitFull e i loro errori
#      dalla matrice di covarianza (usati solo per definire il range di scansione).
#
#   2. Costruisci una griglia 3D di (A, omega0, tau), tenendo v0 e t0 fissi
#      al best fit. Per ogni punto della griglia calcola il chi².
#      → mappa3D[iA, iOm, iTau]
#
#   3. Profila:
#        - 1D: per ogni valore di un parametro, prendi il min sugli altri due
#              (questo è il "profilo likelihood" discretizzato).
#        - 2D: per ogni coppia, prendi il min sul terzo asse
#              (usato per le mappe contour).
#
#   4. Incertezza 1-sigma su ogni parametro = distanza dal best fit al punto
#      in cui il profilo 1D raggiunge chi2_min + 1  (Δchi²= 1).
#
#   I livelli di riferimento per k parametri di interesse sono:
#     - Δchi² = 1.00  → 1-sigma per 1 parametro (68.3 %)
#     - Δchi² = 2.30  → 1-sigma congiunta per 2 parametri (68.3 % in 2D)
#     - Δchi² = 3.84  → 2-sigma per 1 parametro (95 %)
#     - Δchi² = 6.18  → 2-sigma congiunta per 2 parametri (95 % in 2D)
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 8.1  Estrazione dei parametri best fit e del chi2 minimo
# ------------------------------------------------------------------------------

params_bf <- coef(fitResults$fitFull)
errs_bf   <- sqrt(diag(vcov(fitResults$fitFull)))

A_bf   <- params_bf["A"]
om_bf  <- params_bf["omega0"]
tau_bf <- params_bf["tau"]
v0_bf  <- params_bf["v0"]
t0_bf  <- params_bf["t0"]

# pesi usati nel fit
w <- 1 / data$voltUnc^2

# chi2 minimo (dal fit NLS)
chi2_min_fit <- sum(w * residuals(fitResults$fitFull)^2)
cat(sprintf("chi2_min (NLS)    = %.4f\n", chi2_min_fit))
cat(sprintf("dof               = %d\n", nrow(data) - length(params_bf)))


# ------------------------------------------------------------------------------
# 8.2  Definizione della griglia di scansione
# ------------------------------------------------------------------------------

N_SIGMA <- 20L    # ampiezza del range in unità di sigma (da covarianza)
N_STEP  <- 120L   # punti per asse  (80^3 = 512 000 valutazioni del modello)
#
# Aumenta N_STEP per una mappatura più fine, a scapito del tempo di calcolo.
# Con N_STEP=80 e la funzione vettorizzata sotto, il calcolo richiede ~10–30 s
# su un laptop moderno.

A_grid   <- seq(A_bf   - N_SIGMA * errs_bf["A"], A_bf   + N_SIGMA * errs_bf["A"], length.out = N_STEP)

om_grid  <- seq(om_bf  - N_SIGMA * errs_bf["omega0"], om_bf  + N_SIGMA * errs_bf["omega0"], length.out = N_STEP)

tau_grid <- seq(tau_bf - N_SIGMA * errs_bf["tau"],
                tau_bf + N_SIGMA * errs_bf["tau"],
                length.out = N_STEP)


# ------------------------------------------------------------------------------
# 8.3  Calcolo della mappa chi² 3D
# ------------------------------------------------------------------------------
#
# Per ogni (A, omega0, tau) nella griglia calcola:
#
#   chi2(A, omega0, tau) = sum_i  [ (V_i - V_fit(t_i)) / sigma_i ]^2
#
# con v0 e t0 fissati al best fit.
#
# Implementazione vettorizzata tramite outer product implicito:
#   - pre-calcola il modello su (N_STEP × N_dati) matrici
#   - sfrutta il riciclo di R per evitare loop espliciti
#
# NOTA: v0 e t0 sono fissati perché la profilazione completa su 5 parametri
# richiederebbe 80^5 ≈ 3×10^9 valutazioni. Una scansione 3D (come nel codice
# Python allegato) è il compromesso standard per esperimenti didattici.

cat("Calcolo mappa chi2 3D in corso...\n")
t_start <- proc.time()

# Tensioni osservate e pesi come vettori (evita accesso ripetuto al data.frame)
V_obs   <- data$volt
t_obs   <- data$time
sigma_i <- data$voltUnc

# Funzione chi2 per singolo punto della griglia (v0, t0 fissi)
chi2_point <- function(A, omega0, tau) {
  Omega  <- sqrt(max(1e-20, omega0^2 - 1 / tau^2))
  dt     <- t_obs - t0_bf
  V_fit  <- A * exp(-dt / tau) *
    (cos(Omega * dt) + sin(Omega * dt) / (Omega * tau)) + v0_bf
  sum(((V_obs - V_fit) / sigma_i)^2)
}

# Costruisce la griglia come indici e chiama chi2_point su ogni riga
grid_idx <- expand.grid(iA   = seq_len(N_STEP),
                        iOm  = seq_len(N_STEP),
                        iTau = seq_len(N_STEP))

chi2_vec <- mapply(
  function(iA, iOm, iTau)
    chi2_point(A_grid[iA], om_grid[iOm], tau_grid[iTau]),
  grid_idx$iA,
  grid_idx$iOm,
  grid_idx$iTau
)

# Riorganizza in array 3D: dimensione [iA, iOm, iTau]
mappa3D <- array(chi2_vec, dim = c(N_STEP, N_STEP, N_STEP))

t_elapsed <- (proc.time() - t_start)["elapsed"]
cat(sprintf("Completato in %.1f s.\n", t_elapsed))

# Minimo della mappa e sua posizione
chi2_map_min <- min(mappa3D)
idx_min      <- which(mappa3D == chi2_map_min, arr.ind = TRUE)[1L, ]

cat(sprintf("\nchi2_min (mappa)  = %.4f\n", chi2_map_min))
cat(sprintf(
  "  A_best     = %+.4e   [indice griglia %d / %d]\n",
  A_grid[idx_min[1]],
  idx_min[1],
  N_STEP
))
cat(sprintf(
  "  omega0_best= %+.4e   [indice griglia %d / %d]\n",
  om_grid[idx_min[2]],
  idx_min[2],
  N_STEP
))
cat(sprintf(
  "  tau_best   = %+.4e   [indice griglia %d / %d]\n",
  tau_grid[idx_min[3]],
  idx_min[3],
  N_STEP
))


# ------------------------------------------------------------------------------
# 8.4  Profilazione 1D — per ogni parametro, minimo sugli altri due
# ------------------------------------------------------------------------------

# apply(mappa3D, MARGIN, FUN) applica FUN lungo il margine indicato.
# apply(arr, 1, min) → vettore: per ogni iA, min su (iOm, iTau)
prof_A   <- apply(mappa3D, 1L, min)   # profilo su A
prof_om  <- apply(mappa3D, 2L, min)   # profilo su omega0
prof_tau <- apply(mappa3D, 3L, min)   # profilo su tau


# ------------------------------------------------------------------------------
# 8.5  Profilazione 2D — per ogni coppia, minimo sul terzo parametro
# ------------------------------------------------------------------------------

# apply(arr, c(1,2), min) → matrice [iA, iOm]: min su iTau
prof2D_A_om  <- apply(mappa3D, c(1L, 2L), min)   # (A, omega0),  marg. tau
prof2D_A_tau <- apply(mappa3D, c(1L, 3L), min)   # (A, tau),     marg. omega0
prof2D_om_tau <- apply(mappa3D, c(2L, 3L), min)  # (omega0, tau), marg. A


# ------------------------------------------------------------------------------
# 8.6  Calcolo delle incertezze: criterio Δchi² = 1
# ------------------------------------------------------------------------------

#' Trova i limiti ±1-sigma sul profilo 1D del chi².
#'
#' Algoritmo:
#'   1. Localizza il minimo del profilo.
#'   2. A sinistra del minimo: trova l'indice dove il profilo è più vicino
#'      a chi2_min + delta (interpolazione lineare opzionale).
#'   3. Idem a destra.
#'
#' @param prof      Vettore del profilo 1D del chi².
#' @param grid      Vettore delle coordinate del parametro (stessa lunghezza).
#' @param chi2_0    Valore minimo del chi² (può differire leggermente da
#'                  min(prof) se la griglia è rada).
#' @param delta     Soglia in unità di chi² (default 1 → 1-sigma).
#'
#' @return Lista: par_best, par_sx, par_dx, err_sx, err_dx.
findUncertainty <- function(prof, grid, chi2_0, delta = 1.0) {
  level    <- chi2_0 + delta
  idx_best <- which.min(prof)
  
  # --- lato sinistro ---
  if (idx_best > 1L) {
    left_vals <- prof[seq_len(idx_best)]
    idx_sx    <- which.min(abs(left_vals - level))
    # interpolazione lineare tra i due punti che straddlano il livello
    if (idx_sx < idx_best && idx_sx > 1L) {
      d1 <- abs(left_vals[idx_sx - 1L] - level)
      d2 <- abs(left_vals[idx_sx]      - level)
      par_sx <- (grid[idx_sx - 1L] * d2 + grid[idx_sx] * d1) / (d1 + d2)
    } else {
      par_sx <- grid[idx_sx]
    }
  } else {
    warning("Il minimo è al bordo sinistro della griglia — aumenta N_SIGMA.")
    par_sx <- grid[1L]
  }
  
  # --- lato destro ---
  if (idx_best < length(prof)) {
    right_vals <- prof[idx_best:length(prof)]
    idx_dx_rel <- which.min(abs(right_vals - level))
    idx_dx     <- idx_dx_rel + idx_best - 1L
    if (idx_dx > idx_best && idx_dx < length(prof)) {
      d1 <- abs(prof[idx_dx]     - level)
      d2 <- abs(prof[idx_dx + 1L] - level)
      par_dx <- (grid[idx_dx] * d2 + grid[idx_dx + 1L] * d1) / (d1 + d2)
    } else {
      par_dx <- grid[idx_dx]
    }
  } else {
    warning("Il minimo è al bordo destro della griglia — aumenta N_SIGMA.")
    par_dx <- grid[length(grid)]
  }
  
  list(
    par_best = grid[idx_best],
    par_sx   = par_sx,
    par_dx   = par_dx,
    err_sx   = grid[idx_best] - par_sx,
    err_dx   = par_dx - grid[idx_best]
  )
}

unc_A   <- findUncertainty(prof_A, A_grid, chi2_map_min)
unc_om  <- findUncertainty(prof_om, om_grid, chi2_map_min)
unc_tau <- findUncertainty(prof_tau, tau_grid, chi2_map_min)

# Omega e frequenza smorzata derivate dal best fit della mappa
Omega_d <- sqrt(max(0, unc_om$par_best^2 - 1 / unc_tau$par_best^2))

cat("\n====== Incertezze da mappatura chi² (Δchi² = 1) ======\n")
cat(sprintf(
  "  A      = %+.4e  - %.2e  + %.2e\n",
  unc_A$par_best,
  unc_A$err_sx,
  unc_A$err_dx
))
cat(
  sprintf(
    "  omega0 = %+.4e  - %.2e  + %.2e  rad/s\n",
    unc_om$par_best,
    unc_om$err_sx,
    unc_om$err_dx
  )
)
cat(
  sprintf(
    "  tau    = %+.4e  - %.2e  + %.2e  s\n",
    unc_tau$par_best,
    unc_tau$err_sx,
    unc_tau$err_dx
  )
)
cat(sprintf("  Omega_d = %.4e  rad/s\n", Omega_d))
cat(sprintf("  f_d     = %.4e  Hz\n", Omega_d / (2 * pi)))

# Confronto con gli errori dalla matrice di covarianza del fit NLS
cat("\n====== Confronto con errori dalla matrice di covarianza (NLS) ======\n")
cat(sprintf("  sigma(A)      (cov) = %.2e\n", errs_bf["A"]))
cat(sprintf("  sigma(omega0) (cov) = %.2e  rad/s\n", errs_bf["omega0"]))
cat(sprintf("  sigma(tau)    (cov) = %.2e  s\n", errs_bf["tau"]))


# ==============================================================================
# 9. GRAFICI DIAGNOSTICI CON ggplot2
# ==============================================================================
#
# Si producono tre tipi di figura:
#
#  (a) Profili 1D del chi² per A, omega0, tau — con linee Δchi² = 1, 2.3, 3.84.
#  (b) Mappe 2D profilate (contour filled) per le coppie:
#        (omega0, tau), (A, tau), (A, omega0)
#      con i profili 1D a margine (layout "corner plot" minimale).
#  (c) Plot riassuntivo con le bande di incertezza sovrapposte al segnale.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# Livelli di soglia del chi²
LVL <- list(
  sigma1_1D = chi2_map_min + 1.00,
  # 1-sigma,  1 param  (68.3 %)
  sigma1_2D = chi2_map_min + 2.30,
  # 1-sigma,  2 param  (68.3 % in 2D)
  sigma2_1D = chi2_map_min + 3.84,
  # 2-sigma,  1 param  (95 %)
  sigma2_2D = chi2_map_min + 6.18    # 2-sigma,  2 param  (95 % in 2D)
)

# Palette colori per le bande dei contour (da scuro a chiaro)
FILL_COLS  <- c("#2166AC", "#74ADD1", "#FEE08B", "#D73027")
LINE_COLOR <- "grey20"


# ------------------------------------------------------------------------------
# 9.1  Funzione ausiliaria: profilo 1D → ggplot
# ------------------------------------------------------------------------------

#' Costruisce il plot del profilo 1D del chi².
#'
#' @param grid    Vettore coordinate parametro.
#' @param prof    Vettore chi² profilato.
#' @param unc     Lista restituita da findUncertainty().
#' @param xlab    Etichetta asse x (stringa o expression()).
#' @param title   Titolo del plot.
#'
#' @return Oggetto ggplot.
plotProfile1D <- function(grid, prof, unc, xlab, title) {
  df <- data.frame(par = grid, chi2 = prof)
  
  # banda ombreggiata tra i limiti ±1-sigma
  df_band <- df |> filter(par >= unc$par_sx & par <= unc$par_dx)
  
  ggplot(df, aes(x = par, y = chi2)) +
    
    # banda ombreggiata entro Δchi²=1
    geom_ribbon(
      data  = df_band,
      aes(ymin = chi2_map_min, ymax = chi2),
      fill  = "#74ADD1",
      alpha = 0.25
    ) +
    
    # profilo
    geom_line(color = "steelblue", linewidth = 1) +
    
    # linee orizzontali di riferimento
    geom_hline(
      yintercept = LVL$sigma1_1D,
      linetype = "dashed",
      color = "red",
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = LVL$sigma1_2D,
      linetype = "dotted",
      color = "orange",
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = LVL$sigma2_1D,
      linetype = "dotdash",
      color = "purple",
      linewidth = 0.8
    ) +
    
    # linee verticali: best fit e limiti 1-sigma
    geom_vline(
      xintercept = unc$par_best,
      color = "grey30",
      linetype = "solid",
      linewidth = 0.6
    ) +
    geom_vline(
      xintercept = unc$par_sx,
      color = "red",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_vline(
      xintercept = unc$par_dx,
      color = "red",
      linetype = "dashed",
      linewidth = 0.6
    ) +
    
    # etichette Δchi² sui livelli
    annotate(
      "text",
      x = max(grid),
      y = LVL$sigma1_1D,
      label = expression(Delta * chi^2 == 1),
      hjust = 1.1,
      vjust = -0.3,
      size = 3,
      color = "red"
    ) +
    annotate(
      "text",
      x = max(grid),
      y = LVL$sigma2_1D,
      label = expression(Delta * chi^2 == 3.84),
      hjust = 1.1,
      vjust = -0.3,
      size = 3,
      color = "purple"
    ) +
    
    # valore best fit annotato
    annotate(
      "text",
      x = unc$par_best,
      y = chi2_map_min,
      label = sprintf("%.3e", unc$par_best),
      hjust = -0.15,
      vjust = 1.5,
      size = 2.8,
      color = "grey30"
    ) +
    
    labs(title = title,
         x = xlab,
         y = expression(chi^2)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

p_profA <- plotProfile1D(
  A_grid,
  prof_A,
  unc_A,
  xlab  = "A",
  title = expression("Profilo " * chi^2 * " — ampiezza A")
)

p_profOm <- plotProfile1D(
  om_grid,
  prof_om,
  unc_om,
  xlab  = expression(omega[0] * " (rad/s)"),
  title = expression("Profilo " * chi^2 * " — " * omega[0])
)

p_profTau <- plotProfile1D(
  tau_grid,
  prof_tau,
  unc_tau,
  xlab  = expression(tau * " (s)"),
  title = expression("Profilo " * chi^2 * " — " * tau)
)

# Composizione con patchwork
(p_profA | p_profOm | p_profTau) +
  plot_annotation(
    title   = expression("Profili 1D del " * chi^2 * " — transitorio RLC"),
    subtitle = sprintf(
      "chi²_min = %.2f  |  dof = %d  |  griglia %d×%d×%d",
      chi2_map_min,
      nrow(data) - 5L,
      N_STEP,
      N_STEP,
      N_STEP
    ),
    theme   = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  )


# ------------------------------------------------------------------------------
# 9.2  Funzione ausiliaria: mappa 2D → ggplot
# ------------------------------------------------------------------------------

#' Costruisce il plot contour del profilo 2D del chi² con i margini 1D.
#'
#' Layout: mappa centrale + profilo 1D sull'asse y (a sinistra) +
#'         profilo 1D sull'asse x (in basso).
#'
#' @param mat2D   Matrice [ix, iy] del chi² profilato (NON trasposta).
#' @param x_grid  Coordinate dell'asse x.
#' @param y_grid  Coordinate dell'asse y.
#' @param x_unc   Lista incertezza sull'asse x (da findUncertainty).
#' @param y_unc   Lista incertezza sull'asse y (da findUncertainty).
#' @param xlab    Etichetta asse x.
#' @param ylab    Etichetta asse y.
#' @param title   Titolo del pannello.
plotProfile2D <- function(mat2D,
                          x_grid,
                          y_grid,
                          x_unc,
                          y_unc,
                          xlab,
                          ylab,
                          title) {
  # --- 1. Preparazione dati per la mappa ---
  df2D <- expand.grid(x = x_grid, y = y_grid) |>
    mutate(chi2  = as.vector(mat2D), dchi2 = chi2 - chi2_map_min)
  
  # --- 2. Mappa centrale ---
  p_map <- ggplot(df2D, aes(x = x, y = y)) +
    geom_raster(aes(fill = chi2), interpolate = TRUE) +
    
    geom_contour(
      aes(z = chi2),
      bins = 15,
      color = "black",
      alpha = 0.3,
      linewidth = 0.2
    ) +
    
    geom_contour(
      aes(z = dchi2),
      breaks    = c(1.00, 2.30, 3.84, 6.18),
      color     = "black",
      linetype  = "dashed",
      linewidth = 0.6
    ) +
    
    # mirino 1σ
    geom_vline(
      xintercept = c(x_unc$par_sx, x_unc$par_dx),
      color = "grey20",
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.6
    ) +
    geom_hline(
      yintercept = c(y_unc$par_sx, y_unc$par_dx),
      color = "grey20",
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.6
    ) +
    
    geom_point(
      aes(x = x_unc$par_best, y = y_unc$par_best),
      shape = 16,
      size = 2.5,
      color = "black"
    ) +
    
    scale_fill_viridis_c(option = "plasma", name = expression(chi^2)) +
    labs(x = xlab, y = ylab) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "right",
      panel.grid      = element_blank(),
      axis.title      = element_text(size = 9),
      panel.border    = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      )
    )
  
  # --- 3. Profilo 1D asse Y (SINISTRA) ---
  df_py_side <- data.frame(par  = y_grid, chi2 = apply(mat2D, 2L, min)) |>
    arrange(par)
  
  p_side <- ggplot(df_py_side, aes(x = par, y = chi2)) +
    geom_line(color = "steelblue", linewidth = 0.9) +
    
    # linee coerenti con la mappa
    geom_vline(
      xintercept = c(y_unc$par_sx, y_unc$par_dx),
      color = "grey20",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = LVL$sigma1_1D,
      color = "red",
      linetype = "dotted"
    ) +
    
    coord_flip() +
    
    labs(x = expression(chi^2), y = NULL) +
    scale_x_continuous(n.breaks = 3) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), axis.text.y      = element_blank())
  
  # --- 4. Profilo 1D asse X (SOTTO) ---
  df_px_bottom <- data.frame(par  = x_grid, chi2 = apply(mat2D, 1L, min)) |>
    arrange(par)
  
  p_bottom <- ggplot(df_px_bottom, aes(x = par, y = chi2)) +
    geom_line(color = "steelblue", linewidth = 0.9) +
    
    geom_vline(
      xintercept = c(x_unc$par_sx, x_unc$par_dx),
      color = "grey20",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = LVL$sigma1_1D,
      color = "red",
      linetype = "dotted"
    ) +
    
    labs(x = NULL, y = expression(chi^2)) +
    scale_y_continuous(n.breaks = 3) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(), axis.text.x      = element_blank())
  
  # --- 5. Layout finale ---
  p_empty <- ggplot() + theme_void()
  
  (p_side | p_map) / (p_empty | p_bottom) +
    plot_layout(widths = c(1, 4), heights = c(4, 1)) +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(face = "bold", size = 11)))
}

# Mappa (omega0, tau) — profilata su A
p2D_om_tau <- plotProfile2D(
  mat2D  = prof2D_om_tau,
  # t() per avere [iOm, iTau] → x=om, y=tau
  x_grid = om_grid,
  y_grid = tau_grid,
  x_unc  = unc_om,
  y_unc  = unc_tau,
  xlab   = expression(omega[0] * " (rad/s)"),
  ylab   = expression(tau * " (s)"),
  title  = expression("Profilo 2D " * chi^2 * ": " * omega[0] * " vs " * tau * " (marg. su A)")
)
print(p2D_om_tau)

# Mappa (A, tau) — profilata su omega0
p2D_A_tau <- plotProfile2D(
  mat2D  = prof2D_A_tau,
  # [iA, iTau] → x=A, y=tau
  x_grid = A_grid,
  y_grid = tau_grid,
  x_unc  = unc_A,
  y_unc  = unc_tau,
  xlab   = "A",
  ylab   = expression(tau * " (s)"),
  title  = expression("Profilo 2D " * chi^2 * ": A vs " * tau * " (marg. su " * omega[0] * ")")
)
print(p2D_A_tau)

# Mappa (A, omega0) — profilata su tau
p2D_A_om <- plotProfile2D(
  mat2D  = prof2D_A_om,
  # [iA, iOm] → x=A, y=om
  x_grid = A_grid,
  y_grid = om_grid,
  x_unc  = unc_A,
  y_unc  = unc_om,
  xlab   = "A",
  ylab   = expression(omega[0] * " (rad/s)"),
  title  = expression("Profilo 2D " * chi^2 * ": A vs " * omega[0] * " (marg. su " * tau * ")")
)
print(p2D_A_om)


# ------------------------------------------------------------------------------
# 9.3  Plot riassuntivo: segnale + bande di incertezza
# ------------------------------------------------------------------------------
#
# Sovrappone al segnale e al fit best fit un ventaglio di curve ottenute
# variando ciascun parametro di ±1-sigma (dalla mappatura del chi²).
# Questo visualizza graficamente l'effetto delle incertezze sul segnale.

t_fine  <- seq(min(data$time), max(data$time), length.out = 2000)

# Curva best fit dalla mappa (usa i valori al minimo della griglia)
V_bf_map <- modelFull(
  t_fine,
  A      = unc_A$par_best,
  omega0 = unc_om$par_best,
  tau    = unc_tau$par_best,
  v0     = v0_bf,
  t0     = t0_bf
)
df_fit_map <- data.frame(time = t_fine, volt = V_bf_map, tipo = "Best fit (mappa)")

# Genera le curve ai limiti ±1-sigma per omega0 e tau (i parametri fisici chiave)
param_variants <- list(
  list(
    A = unc_A$par_best,
    om = unc_om$par_sx,
    tau = unc_tau$par_best,
    label = "omega0 - 1σ"
  ),
  list(
    A = unc_A$par_best,
    om = unc_om$par_dx,
    tau = unc_tau$par_best,
    label = "omega0 + 1σ"
  ),
  list(
    A = unc_A$par_best,
    om = unc_om$par_best,
    tau = unc_tau$par_sx,
    label = "tau - 1σ"
  ),
  list(
    A = unc_A$par_best,
    om = unc_om$par_best,
    tau = unc_tau$par_dx,
    label = "tau + 1σ"
  )
)

df_variants <- do.call(rbind, lapply(param_variants, function(v) {
  data.frame(
    time = t_fine,
    volt = modelFull(t_fine, v$A, v$om, v$tau, v0_bf, t0_bf),
    tipo = v$label
  )
}))

# Colori per le varianti
variant_cols <- c(
  "omega0 - 1σ" = "#2166AC",
  "omega0 + 1σ" = "#4DAC26",
  "tau - 1σ"    = "#D01C8B",
  "tau + 1σ"    = "#F1B6DA"
)

p_summary <- ggplot() +
  # dati con barre di errore
  geom_errorbar(
    data = data,
    aes(
      x = time,
      ymin = volt - voltUnc,
      ymax = volt + voltUnc
    ),
    width = 0,
    color = "red",
    alpha = 0.15
  ) +
  geom_point(
    data = data,
    aes(x = time, y = volt),
    color = "red",
    size = 0.4,
    alpha = 0.5
  ) +
  # varianti ±1-sigma
  geom_line(
    data = df_variants,
    aes(
      x = time,
      y = volt,
      color = tipo,
      linetype = tipo
    ),
    linewidth = 0.7,
    alpha = 0.85
  ) +
  # best fit
  geom_line(
    data = df_fit_map,
    aes(x = time, y = volt),
    color = "black",
    linewidth = 1,
    linetype = "solid"
  ) +
  scale_color_manual(values = variant_cols) +
  scale_linetype_manual(
    values = c(
      "omega0 - 1σ" = "dashed",
      "omega0 + 1σ" = "dashed",
      "tau - 1σ"    = "dotted",
      "tau + 1σ"    = "dotted"
    )
  ) +
  labs(
    title    = expression("Transitorio RLC — best fit e bande ±1" * sigma * " (mappa " * chi^2 * ")"),
    subtitle = sprintf(
      "A = %.3e  |  ω₀ = %.4e rad/s  |  τ = %.4e s",
      unc_A$par_best,
      unc_om$par_best,
      unc_tau$par_best
    ),
    x     = "Tempo (s)",
    y     = "Tensione (V)",
    color = "Variante",
    linetype = "Variante"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

print(p_summary)


# ------------------------------------------------------------------------------
# 9.4  Tabella riassuntiva delle incertezze (console)
# ------------------------------------------------------------------------------

cat("\n")
cat(strrep("=", 70), "\n")
cat("  RISULTATI FINALI — Transitorio RLC (mappatura chi²)\n")
cat(strrep("=", 70), "\n")
cat(sprintf(
  "  %-8s = %+.4e  - %.2e  + %.2e\n",
  "A",
  unc_A$par_best,
  unc_A$err_sx,
  unc_A$err_dx
))
cat(
  sprintf(
    "  %-8s = %+.4e  - %.2e  + %.2e  rad/s\n",
    "omega0",
    unc_om$par_best,
    unc_om$err_sx,
    unc_om$err_dx
  )
)
cat(
  sprintf(
    "  %-8s = %+.4e  - %.2e  + %.2e  s\n",
    "tau",
    unc_tau$par_best,
    unc_tau$err_sx,
    unc_tau$err_dx
  )
)
cat(sprintf("  %-8s = %+.4e  rad/s  (derivata dal best fit)\n", "Omega_d", Omega_d))
cat(sprintf("  %-8s = %+.4e  Hz\n", "f_d", Omega_d / (2 * pi)))
cat(sprintf("\n  chi2_min (mappa) = %.4f\n", chi2_map_min))
cat(sprintf("  chi2_min (NLS)   = %.4f\n", chi2_min_fit))
cat(sprintf(
  "  Discrepanza      = %.4f  (attesa ~0 se la griglia è sufficientemente fine)\n",
  abs(chi2_map_min - chi2_min_fit)
))
cat(strrep("=", 70), "\n")
