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
assignUncertainty <- function(data, type, scale = 1/25, col = NULL) {
  
  sigmaR <- 0.04 * scale
  
  if (type == "volt") {
    values  <- if (is.data.frame(data)) data[[col]] else data
    sigmaK  <- 0.012 * abs(values)
    uncert  <- sqrt(sigmaR^2 + sigmaK^2)
    
    if (is.data.frame(data)) {
      data[[paste0(col, "Unc")]] <- uncert
      return(data)
    }
    return(uncert)
    
  } else if (type == "time") {
    uncert <- rep(sigmaR, if (is.data.frame(data)) nrow(data) else length(data))
    
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

data <- data.frame(
  time = as.numeric(raw$X4[START_ROW:nrow(raw)]),
  volt = as.numeric(raw$X5[START_ROW:nrow(raw)])
)
rm(raw)

# Aggiunta delle incertezze (scale ricavate dalla configurazione dell'oscilloscopio)
data <- assignUncertainty(data, type = "time", scale = 0.001, col = "time")
data <- assignUncertainty(data, type = "volt", scale = 0.5,   col = "volt")


# ==============================================================================
# 3. MODELLI E FIT NON LINEARE
# ==============================================================================

#' Modello completo di VL(t): include t0 (istante di innesco) e v0 (baseline).
#'
#' Per t <= t0 il segnale vale v0; per t > t0 è un oscillatore smorzato.
modelFull <- function(t, A, omega0, tau, phi, v0, t0) {
  Omega  <- sqrt(pmax(1e-12, omega0^2 - 1/tau^2))
  dt     <- t - t0
  smooth <- 1 / (1 + exp(-50 * dt))   # 50 = ripidezza (regolabile)
  A * exp(-pmax(0, dt)/tau) *
    (cos(Omega * pmax(0, dt) + phi) +
       sin(Omega * pmax(0, dt) + phi) / (Omega * tau)) * smooth + v0
}

#' Modello ridotto di VL(t): usato dopo aver traslato t0->0 e rimosso v0.
modelReduced <- function(t, A, omega0, tau, phi) {
  Omega <- sqrt(pmax(1e-12, omega0^2 - 1/tau^2))
  A * exp(-t/tau) * (cos(Omega*t + phi) + sin(Omega*t + phi) / (Omega*tau))
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
  
  cat(sprintf("  %-8s = %+.4e  ±  %.2e\n",
              names(params), params, errs), sep = "")
  cat(sprintf("  chi2 / dof  =  %.2f / %d  =  %.3f\n",
              chi2, dof, chi2/dof))
  
  if (step == 2) {
    Omega <- sqrt(max(0, params["omega0"]^2 - 1/params["tau"]^2))
    cat(sprintf("\n  Omega  = %.4e  rad/s\n", Omega))
    cat(sprintf("  f_d    = %.4e  Hz\n",    Omega / (2*pi)))
    cat(sprintf("  tau    = %.4e  s\n",     params["tau"]))
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
                  omega0Init = 1, tauInit = 1, AInit = 1,
                  phiInit = 0,   v0Init  = 0, t0Init = 0,
                  fixT0 = FALSE) {
  
  nlsControl <- nls.control(maxiter = 500, tol = 1e-8)
  
  if (!fixT0) {
    
    # --- Passo 1: stima t0 e v0 ---
    fitFull <- nls(
      volt ~ modelFull(time, A, omega0, tau, phi, v0, t0),
      data    = data,
      start   = list(A = AInit, omega0 = omega0Init, tau = tauInit,
                     phi = phiInit, v0 = v0Init, t0 = t0Init),
      weights = 1 / data$voltUnc^2,
      control = nlsControl
    )
    
    t0BF <- coef(fitFull)["t0"]
    v0BF <- coef(fitFull)["v0"]
    
    cat("====== PASSO 1: fit completo ======\n")
    printFitResults(fitFull, step = 1)
    
    # Rimozione baseline e traslazione temporale
    dataCut       <- data[data$time > t0BF, ]
    dataCut$volt  <- dataCut$volt - v0BF
    dataCut$time  <- dataCut$time - t0BF
    
    # --- Passo 2: fit ridotto ---
    fitReduced <- nls(
      volt ~ modelReduced(time, A, omega0, tau, phi),
      data    = dataCut,
      start   = list(A = AInit, omega0 = omega0Init,
                     tau = tauInit, phi = phiInit),
      weights = 1 / dataCut$voltUnc^2,
      control = nlsControl
    )
    
    cat("\n====== PASSO 2: fit ridotto (t > t0, baseline rimossa) ======\n")
    printFitResults(fitReduced, step = 2)
    
    return(invisible(list(fitFull = fitFull, fitReduced = fitReduced, dataCut = dataCut)))
    
  } else {
    
    # --- Solo fit ridotto (t0 già noto e rimosso) ---
    fitReduced <- nls(
      volt ~ modelReduced(time, A, omega0, tau, phi),
      data    = data,
      start   = list(A = AInit, omega0 = omega0Init,
                     tau = tauInit, phi = phiInit),
      weights = 1 / data$voltUnc^2,
      control = nlsControl
    )
    
    cat("====== Fit ridotto ======\n")
    printFitResults(fitReduced, step = 2)
    
    return(invisible(list(fitReduced = fitReduced, dataCut = data)))
  }
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
estimateTransientInitialsBrutal <- function(data, kNeighbors = 25,
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
    while (j < n - kNeighbors && tensione[j + 1] == tensione[i]) j <- j + 1
    
    val    <- tensione[i]
    sinSx  <- tensione[(i - kNeighbors):(i - 1)]
    sinDx  <- tensione[(j + 1):(j + kNeighbors)]
    idxMid <- round((i + j) / 2)
    
    if (all(val >= sinSx) && all(val >= sinDx)) {
      if (length(indiceMassimi) == 0 || idxMid > max(indiceMassimi) + kNeighbors)
        indiceMassimi <- c(indiceMassimi, idxMid)
    }
    
    if (all(val <= sinSx) && all(val <= sinDx)) {
      if (length(indiceMinimi) == 0 || idxMid > max(indiceMinimi) + kNeighbors)
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
    if (length(indiceEstremi) < 2) return(NULL)
    
    t <- tempo[indiceEstremi]
    v <- tensioneCentrata[indiceEstremi]
    
    T_i <- diff(t)
    v1  <- abs(v[-length(v)])
    v2  <- abs(v[-1])
    
    ok <- is.finite(T_i) & T_i > 0 &
      is.finite(v1) & is.finite(v2) &
      v1 > 0 & v2 > 0 & v1 > v2
    
    if (!any(ok)) return(NULL)
    
    T_i <- T_i[ok]
    v1  <- v1[ok]
    v2  <- v2[ok]
    
    omegaPseudo_i <- 2 * pi / T_i
    tau_i <- T_i / log(v1 / v2)
    omega0_i <- sqrt(omegaPseudo_i^2 + 1 / tau_i^2)
    
    list(T = T_i, tau = tau_i, omegaPseudo = omegaPseudo_i, omega0 = omega0_i)
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
  
  sigmaDer <- sd(dVdt[1:(nInit-1)], na.rm = TRUE)
  threshold <- 5 * sigmaDer
  
  idx <- which(abs(dVdt) > threshold)
  
  if (length(idx) > 0) {
    stimaT0 <- tempoMid[idx[1]]
  } else {
    stimaT0 <- median(tempo[1:nInit])
  }
  
  cat(sprintf("Picchi validi trovati: max=%d, min=%d\n",
              length(indiceMassimi), length(indiceMinimi)))
  
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
  data$time, data$volt,
  col  = "red", pch = 16, cex = 0.5,
  xlab = "Tempo (s)", ylab = "Tensione (V)",
  main = "Segnale VL(t) e picchi identificati"
)
points(dataPMax$time, dataPMax$volt, col = "blue",  pch = 16, cex = 0.8)
points(dataPMin$time, dataPMin$volt, col = "green", pch = 16, cex = 0.8)
legend("topright",
       legend = c("Dati grezzi", "Massimi", "Minimi"),
       col    = c("red", "blue", "green"), pch = 16)

# ==============================================================================
# 6. FIT NON LINEARE
# ==============================================================================

fitResults <- fitVL(
  data       = data,
  AInit      = initialGuesses$aInit,
  omega0Init = initialGuesses$omega0Init,
  tauInit    = initialGuesses$tauInit,
  phiInit    = initialGuesses$phiInit,
  v0Init     = initialGuesses$v0Init,
  t0Init     = initialGuesses$t0Init,
  fixT0      = FALSE
)
