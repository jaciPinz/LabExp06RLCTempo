# ==============================================================================
# Analisi transitorio RLC — tensione ai capi dell'induttore
# Codice riordinato, camelCase, ottimizzato per calcolo fine del χ²
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
  library(viridis)
  library(metR)
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
assignUncertainty <- function(data, type, scale = 1 / 25, col = NULL) {
  sigmaR <- 0.04 * scale
  
  if (type == "volt") {
    values <- if (is.data.frame(data)) data[[col]] else data
    sigmaK <- 0.012 * abs(values)
    uncert <- sqrt(sigmaR^2 + sigmaK^2)
    
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

START_ROW <- 176  # riga da cui iniziano i dati veri

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
# 3. MODELLI
# ==============================================================================

#' Modello completo di VL(t): include t0 (istante di innesco) e v0 (baseline).
modelFull <- function(t, A, omega0, tau, v0, t0) {
  Omega <- sqrt(pmax(1e-12, omega0^2 - 1 / tau^2))
  dt    <- t - t0
  Vosc  <- A * exp(-dt / tau) * (cos(Omega * dt) + sin(Omega * dt) / (Omega * tau))
  return(Vosc + v0)
}

#' Modello ridotto di VL(t): dopo traslazione t0->0 e rimozione di v0.
modelReduced <- function(t, A, omega0, tau) {
  Omega <- sqrt(pmax(0, omega0^2 - 1 / tau^2))
  A * exp(-t / tau) * (cos(Omega * t) + sin(Omega * t) / (Omega * tau))
}

# ==============================================================================
# 4. FUNZIONI DI SUPPORTO PER FIT E STAMPA
# ==============================================================================

#' Stampa i risultati di un fit nls in modo compatto.
printFitResults <- function(fit, step = 2) {
  params <- coef(fit)
  errs   <- sqrt(diag(vcov(fit)))
  resid  <- residuals(fit)
  w      <- fit$weights %||% rep(1, length(resid))
  
  chi2 <- sum(w * resid^2)
  dof  <- length(resid) - length(params)
  
  for (i in seq_along(params)) {
    cat(sprintf("  %-8s = %+.4e  ±  %.2e\n", names(params)[i], params[i], errs[i]))
  }
  cat(sprintf("  chi2 / dof  =  %.2f / %d  =  %.3f\n", chi2, dof, chi2 / dof))
  
  if (step == 2) {
    Omega <- sqrt(max(0, params["omega0"]^2 - 1 / params["tau"]^2))
    cat(sprintf("\n  Omega  = %.4e  rad/s\n", Omega))
    cat(sprintf("  f_d    = %.4e  Hz\n",    Omega / (2 * pi)))
    cat(sprintf("  tau    = %.4e  s\n",     params["tau"]))
  }
}

# ==============================================================================
# 5. STIMA AUTOMATICA DEI PARAMETRI INIZIALI
# ==============================================================================

#' Stima grezza dei parametri iniziali cercando i picchi del segnale.
estimateTransientInitialsBrutal <- function(data, kNeighbors = 25, sogliaPercentuale = 0.1) {
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
  
  # Filtro ampiezza sui massimi e minimi
  if (length(indiceMassimi) > 0) {
    soglia <- abs(tensione[indiceMassimi[1]]) * sogliaPercentuale
    indiceMassimi <- indiceMassimi[abs(tensione[indiceMassimi]) > soglia]
  }
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
  
  T_all <- tau_all <- omegaPseudo_all <- omega0_all <- numeric(0)
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
  
  # Ampiezza iniziale
  if (length(indiceMassimi) >= 1) {
    stimaA <- abs(tensioneCentrata[indiceMassimi[1]])
  } else if (length(indiceMinimi) >= 1) {
    stimaA <- abs(tensioneCentrata[indiceMinimi[1]])
  } else {
    stimaA <- max(abs(tensioneCentrata))
  }
  
  # Stima robusta di t0 basata sulla derivata iniziale
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
  
  cat(sprintf("Picchi validi trovati: max=%d, min=%d\n",
              length(indiceMassimi), length(indiceMinimi)))
  
  list(
    aInit       = stimaA,
    omega0Init  = stimaOmega0,
    tauInit     = stimaTau,
    phiInit     = 0,
    v0Init      = stimaV0,
    t0Init      = stimaT0,
    periodoMedio = periodoMedio,
    indiciMax    = indiceMassimi,
    indiciMin    = indiceMinimi
  )
}

# ==============================================================================
# 6. FIT NON LINEARE (MODELLO COMPLETO)
# ==============================================================================

fitVL <- function(data, AInit, omega0Init, tauInit, v0Init, t0Init) {
  nlsControl <- nls.control(maxiter = 500, tol = 1e-8)
  
  fitFull <- nlsLM(
    volt ~ modelFull(time, A, omega0, tau, v0, t0),
    data    = data,
    start   = list(A = AInit, omega0 = omega0Init, tau = tauInit, v0 = v0Init, t0 = t0Init),
    weights = 1 / data$voltUnc^2,
    control = nlsControl
  )
  
  cat("====== Fit completo ======\n")
  printFitResults(fitFull, step = 2)
  
  invisible(list(fitFull = fitFull))
}

# ==============================================================================
# 7. ESECUZIONE FIT E PLOT DIAGNOSTICO INIZIALE
# ==============================================================================

initialGuesses <- estimateTransientInitialsBrutal(data, kNeighbors = 50)

# Plot diagnostico: dati completi + massimi + minimi trovati
plot(data$time, data$volt,
     col  = "red", pch = 16, cex = 0.5,
     xlab = "Tempo (s)", ylab = "Tensione (V)",
     main = "Segnale VL(t) e picchi identificati")
points(data$time[initialGuesses$indiciMax], data$volt[initialGuesses$indiciMax],
       col = "blue", pch = 16, cex = 0.8)
points(data$time[initialGuesses$indiciMin], data$volt[initialGuesses$indiciMin],
       col = "green", pch = 16, cex = 0.8)
legend("topright", legend = c("Dati grezzi", "Massimi", "Minimi"),
       col = c("red", "blue", "green"), pch = 16)

fitResults <- fitVL(
  data       = data,
  AInit      = initialGuesses$aInit,
  omega0Init = initialGuesses$omega0Init,
  tauInit    = initialGuesses$tauInit,
  v0Init     = initialGuesses$v0Init,
  t0Init     = initialGuesses$t0Init
)

# ==============================================================================
# 8. PLOT DIAGNOSTICO CON RESIDUI REALI E BARRE D'ERRORE
# ==============================================================================

params <- coef(fitResults$fitFull)
data$fitValues <- modelFull(data$time, params["A"], params["omega0"], params["tau"],
                            params["v0"], params["t0"])
data$residReal <- data$volt - data$fitValues

tGrid <- seq(min(data$time), max(data$time), length.out = 2000)
fitDF <- data.frame(time = tGrid,
                    volt = modelFull(tGrid, params["A"], params["omega0"], params["tau"],
                                     params["v0"], params["t0"]))

pSignal <- ggplot(data, aes(x = time, y = volt)) +
  geom_errorbar(aes(ymin = volt - voltUnc, ymax = volt + voltUnc),
                width = 0, color = "red", alpha = 0.2) +
  geom_point(color = "red", size = 0.5, alpha = 0.4) +
  geom_line(data = fitDF, aes(x = time, y = volt), color = "black", linewidth = 0.8) +
  geom_point(data = data[initialGuesses$indiciMax, ], color = "blue", size = 2) +
  geom_point(data = data[initialGuesses$indiciMin, ], color = "green", size = 2) +
  labs(title = "Transitorio RLC — Fit e Incertezze", x = "Tempo (s)", y = "Tensione (V)") +
  theme_minimal()

pResid <- ggplot(data, aes(x = time, y = residReal)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = residReal - voltUnc, ymax = residReal + voltUnc),
                width = 0, color = "purple", alpha = 0.2) +
  geom_point(color = "purple", size = 0.7, alpha = 0.5) +
  labs(title = "Residui", x = "Tempo (s)", y = "Residui (V)") +
  theme_minimal()

print(pSignal / pResid)

# ==============================================================================
# 9. FIT RIDOTTO E MAPPATURA COMPLETA DEL χ² (3D)
# ==============================================================================

# --- 9.1 Fit ridotto (solo A, omega0, tau) ---
pFull <- coef(fitResults$fitFull)
t0Est <- pFull["t0"]
v0Est <- pFull["v0"]

dataCut <- data[data$time > t0Est, ]
dataCut$timeShifted <- dataCut$time - t0Est
dataCut$voltShifted <- dataCut$volt - v0Est

fitRed <- nlsLM(
  voltShifted ~ modelReduced(timeShifted, A, omega0, tau),
  data    = dataCut,
  start   = list(A = pFull["A"], omega0 = pFull["omega0"], tau = pFull["tau"]),
  weights = 1 / dataCut$voltUnc^2,
  control = nls.control(maxiter = 500, tol = 1e-8)
)

cat("\n====== Fit ridotto (dopo rimozione t0 e v0) ======\n")
printFitResults(fitRed, step = 2)

pRed <- coef(fitRed)
eRed <- sqrt(diag(vcov(fitRed)))
Abf   <- pRed[1]
om0bf <- pRed[2]
taubf <- pRed[3]
eAbf   <- eRed[1]
eOm0bf <- eRed[2]
eTaubf <- eRed[3]

# --- 9.2 Preparazione dati e griglie ---
x <- dataCut$timeShifted
y <- dataCut$voltShifted
w <- 1 / dataCut$voltUnc^2

NSI <- 20
nStep <- 150   # finezza della griglia (puoi aumentare)

om0Lim <- c(om0bf - NSI * eOm0bf, om0bf + NSI * eOm0bf)
tauLim <- c(taubf - NSI * eTaubf, taubf + NSI * eTaubf)
ALim   <- c(Abf   - NSI * eAbf,   Abf   + NSI * eAbf)

om0Grid <- seq(om0Lim[1], om0Lim[2], length.out = nStep)
tauGrid <- seq(tauLim[1], tauLim[2], length.out = nStep)
AGrid   <- seq(ALim[1],   ALim[2],   length.out = nStep)

# --- 9.3 Calcolo vettoriale della mappa χ² 3D (ottimizzato) ---
# Per ogni coppia (omega0, tau) calcoliamo il vettore f e il χ² minimo rispetto ad A.
# Inoltre, per una griglia completa 3D calcoliamo il χ² per tutti i valori di A.

gridOmTau <- expand.grid(omega0 = om0Grid, tau = tauGrid)
nOmTau <- nrow(gridOmTau)

# Matrice exp(-x/tau) per ogni tau
expMat <- outer(x, tauGrid, function(x, tau) exp(-x / tau))

cat("Calcolo mappa χ² 3D (ottimizzata) ...\n")
pb <- txtProgressBar(min = 0, max = nOmTau, style = 3)

# Pre-allocazione
chi2MinA <- numeric(nOmTau)        # χ² minimo rispetto ad A
AottVec  <- numeric(nOmTau)        # A ottimale corrispondente

for (k in seq_len(nOmTau)) {
  om0 <- gridOmTau$omega0[k]
  tau <- gridOmTau$tau[k]
  Omega <- sqrt(max(0, om0^2 - 1/tau^2))
  
  idxTau <- which(tauGrid == tau)[1]
  f <- expMat[, idxTau] * (cos(Omega * x) + sin(Omega * x) / (Omega * tau))
  
  wf <- w * f
  sumWf2 <- sum(wf * f)
  Aott <- sum(wf * y) / sumWf2
  
  AottVec[k] <- Aott
  chi2MinA[k] <- sum(w * (y - Aott * f)^2)
  
  setTxtProgressBar(pb, k)
}
close(pb)

gridOmTau$Aott <- AottVec
gridOmTau$chi2 <- chi2MinA

# --- 9.4 Griglia 3D completa (A, omega0, tau) ---
# Per ogni punto della griglia 3D calcoliamo il χ² esatto.
# Sfruttiamo il fatto che per fissati (omega0, tau) il χ² è una parabola in A:
# χ²(A) = χ²_min + sumWf2 * (A - A_ott)^2
# Questo evita di ricalcolare il modello per ogni A.

cat("\nCostruzione griglia 3D completa (A, omega0, tau) ...\n")
grid3dFull <- expand.grid(A = AGrid, omega0 = om0Grid, tau = tauGrid)
n3d <- nrow(grid3dFull)

# Per ogni punto della griglia 3d dobbiamo trovare il corrispondente Aott e sumWf2
# Creiamo un indice che mappa (omega0, tau) alla riga di gridOmTau
gridOmTau$idx <- 1:nOmTau
grid3dFull <- merge(grid3dFull, gridOmTau[, c("omega0", "tau", "Aott", "chi2", "idx")],
                    by = c("omega0", "tau"), all.x = TRUE)

# Calcoliamo sumWf2 per ogni coppia (omega0, tau)
# sumWf2 = sum( w * f^2 )  --> lo ricalcoliamo una volta per ogni coppia
cat("Calcolo di sumWf2 per ogni (omega0, tau) ...\n")
pb <- txtProgressBar(min = 0, max = nOmTau, style = 3)
sumWf2Vec <- numeric(nOmTau)
for (k in seq_len(nOmTau)) {
  om0 <- gridOmTau$omega0[k]
  tau <- gridOmTau$tau[k]
  Omega <- sqrt(max(0, om0^2 - 1/tau^2))
  idxTau <- which(tauGrid == tau)[1]
  f <- expMat[, idxTau] * (cos(Omega * x) + sin(Omega * x) / (Omega * tau))
  wf <- w * f
  sumWf2Vec[k] <- sum(wf * f)
  setTxtProgressBar(pb, k)
}
close(pb)
gridOmTau$sumWf2 <- sumWf2Vec

# Merge di sumWf2 nel grid3dFull
grid3dFull <- merge(grid3dFull, gridOmTau[, c("omega0", "tau", "sumWf2")],
                    by = c("omega0", "tau"), all.x = TRUE)

# Ora calcoliamo χ² = chi2_min + sumWf2 * (A - A_ott)^2
grid3dFull$chi2 <- grid3dFull$chi2 + grid3dFull$sumWf2 * (grid3dFull$A - grid3dFull$Aott)^2

# Rimuoviamo colonne di servizio
grid3dFull$idx <- NULL
grid3dFull$sumWf2 <- NULL

# --- 9.5 Profilazione 1D per A, omega0, tau ---
profile1d <- function(df, paramName) {
  prof <- aggregate(chi2 ~ get(paramName), data = df, FUN = min)
  names(prof) <- c("paramValue", "chi2")
  chi2min <- min(prof$chi2)
  target <- chi2min + 1.0
  idx <- which(prof$chi2 <= target)
  if (length(idx) > 0) {
    lower <- prof$paramValue[min(idx)]
    upper <- prof$paramValue[max(idx)]
  } else {
    lower <- upper <- NA
  }
  list(profile = prof, chi2min = chi2min, lower = lower, upper = upper)
}

profOm0 <- profile1d(grid3dFull, "omega0")
profTau <- profile1d(grid3dFull, "tau")
profA   <- profile1d(grid3dFull, "A")

cat("\n========== RISULTATI MAPPA χ² (modello ridotto) ==========\n")
cat(sprintf("A      = %.4e  + %.1e / - %.1e\n",
            Abf, profA$upper - Abf, Abf - profA$lower))
cat(sprintf("omega0 = %.4e  + %.1e / - %.1e\n",
            om0bf, profOm0$upper - om0bf, om0bf - profOm0$lower))
cat(sprintf("tau    = %.4e  + %.1e / - %.1e\n",
            taubf, profTau$upper - taubf, taubf - profTau$lower))
cat(sprintf("χ²_min = %.3f\n", profOm0$chi2min))
cat("==========================================================\n")

# Salviamo la griglia 3D completa e i profili per le sezioni successive
# (in una sessione interattiva questi oggetti restano nell'ambiente)

# ==============================================================================
# 10. VISUALIZZAZIONE CONTOUR E PROFILI (MODELLO RIDOTTO) – TUTTE LE COPPIE
# ==============================================================================

# Funzione per generare il layout combinato per una data coppia di parametri
plotPairReduced <- function(grid3d, paramX, paramY, 
                            profX, profY, 
                            xlab, ylab, title = NULL) {
  
  # Aggrega la griglia 3D sui due parametri (minimizza sul terzo)
  grid2d <- aggregate(chi2 ~ get(paramX) + get(paramY), data = grid3d, FUN = min)
  names(grid2d)[1:2] <- c("x", "y")
  
  chi2Min2d <- min(grid2d$chi2)
  contourLevels <- c(chi2Min2d + 0.0001, chi2Min2d + 1, 
                     chi2Min2d + 2.3, chi2Min2d + 3.8)
  
  xLim <- range(grid2d$x)
  yLim <- range(grid2d$y)
  
  # Mappa centrale
  pMap <- ggplot(grid2d, aes(x = x, y = y, z = chi2)) +
    geom_contour_filled(bins = 100) +
    geom_contour(breaks = contourLevels, color = "black", alpha = 0.5, linetype = "dotted") +
    scale_fill_viridis_d(option = "plasma", direction = -1, name = expression(chi^2)) +
    geom_vline(xintercept = c(profX$lower, profX$upper), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = c(profY$lower, profY$upper), linetype = "dashed", color = "gray40") +
    coord_cartesian(xlim = xLim, ylim = yLim, expand = FALSE) +
    labs(x = xlab, y = ylab) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
  
  # Profilo laterale sinistro (parametro Y)
  pSide <- ggplot(profY$profile, aes(x = paramValue, y = chi2)) +
    geom_line(color = "steelblue", linewidth = 0.9) +
    geom_vline(xintercept = c(profY$lower, profY$upper), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = chi2Min2d + 1, linetype = "dotted", color = "red") +
    coord_flip(xlim = yLim, ylim = range(profY$profile$chi2)) +
    labs(x = ylab, y = expression(chi^2)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
  
  # Profilo inferiore (parametro X)
  pBottom <- ggplot(profX$profile, aes(x = paramValue, y = chi2)) +
    geom_line(color = "steelblue", linewidth = 0.9) +
    geom_vline(xintercept = c(profX$lower, profX$upper), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = chi2Min2d + 1, linetype = "dotted", color = "red") +
    coord_cartesian(xlim = xLim, ylim = range(profX$profile$chi2)) +
    labs(x = xlab, y = expression(chi^2)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
  
  pEmpty <- ggplot() + theme_void()
  
  layout <- (pSide | pMap) / (pEmpty | pBottom) +
    plot_layout(widths = c(1, 4), heights = c(4, 1))
  
  if (!is.null(title)) {
    layout <- layout + plot_annotation(title = title,
                                       theme = theme(plot.title = element_text(face = "bold", size = 11)))
  }
  return(layout)
}

# Definizione delle coppie da plottare
pairsReduced <- list(
  list(x = "omega0", y = "tau",   xProf = profOm0, yProf = profTau,
       xlab = expression(omega[0]~(rad/s)), ylab = expression(tau~(s)),
       title = "Modello ridotto – ω₀ vs τ"),
  list(x = "A",      y = "tau",   xProf = profA,   yProf = profTau,
       xlab = "A", ylab = expression(tau~(s)),
       title = "Modello ridotto – A vs τ"),
  list(x = "A",      y = "omega0", xProf = profA,   yProf = profOm0,
       xlab = "A", ylab = expression(omega[0]~(rad/s)),
       title = "Modello ridotto – A vs ω₀")
)

# Generazione e visualizzazione di tutti i layout
plotsReduced <- list()
for (pair in pairsReduced) {
  cat("Generazione plot:", pair$title, "\n")
  plotsReduced[[pair$title]] <- plotPairReduced(
    grid3d = grid3dFull,
    paramX = pair$x, paramY = pair$y,
    profX  = pair$xProf, profY = pair$yProf,
    xlab   = pair$xlab, ylab = pair$ylab,
    title  = pair$title
  )
  print(plotsReduced[[pair$title]])
}

# Salvataggio opzionale
# ggsave("chi2_reduced_omega0_tau.png", plotsReduced[[1]], width = 8, height = 6)
# ggsave("chi2_reduced_A_tau.png",     plotsReduced[[2]], width = 8, height = 6)
# ggsave("chi2_reduced_A_omega0.png",  plotsReduced[[3]], width = 8, height = 6)

# ==============================================================================
# 11. PROFILAZIONE PER IL MODELLO COMPLETO (5 PARAMETRI) – TUTTE LE COPPIE
# ==============================================================================

# Funzione χ² per il modello completo
chi2Full <- function(par, data, fixPar = NULL, fixNames = NULL) {
  allPar <- numeric(5)
  names(allPar) <- c("A", "omega0", "tau", "v0", "t0")
  if (!is.null(fixPar)) allPar[fixNames] <- fixPar
  freeIdx <- setdiff(names(allPar), fixNames)
  allPar[freeIdx] <- par
  
  yFit <- modelFull(data$time, allPar["A"], allPar["omega0"], allPar["tau"],
                    allPar["v0"], allPar["t0"])
  w <- 1 / data$voltUnc^2
  sum(w * (data$volt - yFit)^2)
}

profile1dFull <- function(paramName, gridValues, data, bestFit) {
  fixNames <- paramName
  freeNames <- setdiff(names(bestFit), fixNames)
  start <- bestFit[freeNames]
  
  chi2Prof <- numeric(length(gridValues))
  for (i in seq_along(gridValues)) {
    fixVal <- gridValues[i]
    opt <- optim(par = start, fn = chi2Full, data = data,
                 fixPar = fixVal, fixNames = fixNames,
                 method = "BFGS", control = list(maxit = 1000, reltol = 1e-8))
    chi2Prof[i] <- opt$value
    start <- opt$par
  }
  data.frame(paramValue = gridValues, chi2 = chi2Prof)
}

profile2dFull <- function(paramNames, grid1, grid2, data, bestFit) {
  fixNames <- paramNames
  freeNames <- setdiff(names(bestFit), fixNames)
  start <- bestFit[freeNames]
  
  ng1 <- length(grid1)
  ng2 <- length(grid2)
  chi2Mat <- matrix(NA, nrow = ng1, ncol = ng2)
  
  pb <- txtProgressBar(min = 0, max = ng1 * ng2, style = 3)
  counter <- 0
  for (i in seq_along(grid1)) {
    for (j in seq_along(grid2)) {
      fixVal <- c(grid1[i], grid2[j])
      opt <- optim(par = start, fn = chi2Full, data = data,
                   fixPar = fixVal, fixNames = fixNames,
                   method = "BFGS", control = list(maxit = 1000, reltol = 1e-8))
      chi2Mat[i, j] <- opt$value
      start <- opt$par
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }
  close(pb)
  
  df <- expand.grid(p1 = grid1, p2 = grid2)
  names(df) <- paramNames
  df$chi2 <- as.vector(chi2Mat)
  return(df)
}

# Parametri best‑fit ed errori dal fit completo
bestFull <- coef(fitResults$fitFull)
names(bestFull) <- c("A", "omega0", "tau", "v0", "t0")
errFull <- sqrt(diag(vcov(fitResults$fitFull)))
names(errFull) <- names(bestFull)

# χ² minimo dal fit completo
wFull <- 1 / data$voltUnc^2
residFull <- residuals(fitResults$fitFull)
chi2MinFull <- sum(wFull * residFull^2)

NSI_full <- 20
nStep <- 150   # risoluzione per le mappe 2D (aumentabile)

cat("\n====== PROFILAZIONE MODELLO COMPLETO (5 PARAMETRI) ======\n")
cat("χ² minimo (fit completo):", chi2MinFull, "\n")
cat("Intervalli di scansione: ±", NSI_full, "σ\n")

# ---------- Profilazione 1D per tutti i parametri ----------
paramNames <- c("A", "omega0", "tau", "v0", "t0")
profiles1D <- list()

for (pname in paramNames) {
  cat("Profilo 1D per", pname, "...\n")
  gridVals <- seq(bestFull[pname] - NSI_full * errFull[pname],
                  bestFull[pname] + NSI_full * errFull[pname],
                  length.out = 50)
  dfProf <- profile1dFull(pname, gridVals, data, bestFull)
  chi2min <- min(dfProf$chi2)
  idx <- which(dfProf$chi2 <= chi2min + 1)
  lower <- if (length(idx) > 0) dfProf$paramValue[min(idx)] else NA
  upper <- if (length(idx) > 0) dfProf$paramValue[max(idx)] else NA
  
  profiles1D[[pname]] <- list(df = dfProf, lower = lower, upper = upper)
  
  cat(sprintf("  %s = %.4e  + %.1e / - %.1e\n",
              pname, bestFull[pname], upper - bestFull[pname], bestFull[pname] - lower))
}

# ---------- Funzione per generare layout da matrice 2D e profili ----------
plotPairFull <- function(mat2D, xGrid, yGrid,
                         profX, profY,
                         xlab, ylab, title = NULL) {
  
  df2D <- expand.grid(x = xGrid, y = yGrid)
  df2D$chi2 <- as.vector(mat2D)
  chi2Min2d <- min(df2D$chi2, na.rm = TRUE)
  contourLevels <- c(chi2Min2d + 0.0001, chi2Min2d + 1,
                     chi2Min2d + 2.3, chi2Min2d + 3.8)
  
  xLim <- range(xGrid)
  yLim <- range(yGrid)
  
  pMap <- ggplot(df2D, aes(x = x, y = y, z = chi2)) +
    geom_contour_filled(bins = 100) +
    geom_contour(breaks = contourLevels, color = "black", alpha = 0.5, linetype = "dotted") +
    scale_fill_viridis_d(option = "plasma", direction = -1, name = expression(chi^2)) +
    geom_vline(xintercept = c(profX$lower, profX$upper), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = c(profY$lower, profY$upper), linetype = "dashed", color = "gray40") +
    coord_cartesian(xlim = xLim, ylim = yLim, expand = FALSE) +
    labs(x = xlab, y = ylab) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
  
  pSide <- ggplot(profY$df, aes(x = paramValue, y = chi2)) +
    geom_line(color = "steelblue", linewidth = 0.9) +
    geom_vline(xintercept = c(profY$lower, profY$upper), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = chi2Min2d + 1, linetype = "dotted", color = "red") +
    coord_flip(xlim = yLim, ylim = range(profY$df$chi2)) +
    labs(x = ylab, y = expression(chi^2)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
  
  pBottom <- ggplot(profX$df, aes(x = paramValue, y = chi2)) +
    geom_line(color = "steelblue", linewidth = 0.9) +
    geom_vline(xintercept = c(profX$lower, profX$upper), linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = chi2Min2d + 1, linetype = "dotted", color = "red") +
    coord_cartesian(xlim = xLim, ylim = range(profX$df$chi2)) +
    labs(x = xlab, y = expression(chi^2)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(0, 0, 0, 0))
  
  pEmpty <- ggplot() + theme_void()
  
  layout <- (pSide | pMap) / (pEmpty | pBottom) +
    plot_layout(widths = c(1, 4), heights = c(4, 1))
  
  if (!is.null(title)) {
    layout <- layout + plot_annotation(title = title,
                                       theme = theme(plot.title = element_text(face = "bold", size = 11)))
  }
  return(layout)
}

# ---------- Profilazione 2D per le coppie di interesse ----------
pairsFull <- list(
  list(x = "omega0", y = "tau",   xlab = expression(omega[0]~(rad/s)), ylab = expression(tau~(s)),
       title = "Modello completo – ω₀ vs τ"),
  list(x = "A",      y = "tau",   xlab = "A", ylab = expression(tau~(s)),
       title = "Modello completo – A vs τ"),
  list(x = "A",      y = "omega0", xlab = "A", ylab = expression(omega[0]~(rad/s)),
       title = "Modello completo – A vs ω₀"),
  list(x = "v0",     y = "t0",    xlab = expression(v[0]~(V)), ylab = expression(t[0]~(s)),
       title = "Modello completo – v₀ vs t₀")
)

plotsFull <- list()

for (pair in pairsFull) {
  pX <- pair$x
  pY <- pair$y
  cat("Mappa 2D per", pX, "vs", pY, "...\n")
  
  xVals <- seq(bestFull[pX] - NSI_full * errFull[pX],
               bestFull[pX] + NSI_full * errFull[pX],
               length.out = nStep)
  yVals <- seq(bestFull[pY] - NSI_full * errFull[pY],
               bestFull[pY] + NSI_full * errFull[pY],
               length.out = nStep)
  
  mat2D <- profile2dFull(c(pX, pY), xVals, yVals, data, bestFull)
  # mat2D è un data.frame con colonne pX, pY, chi2
  # Convertiamo in matrice per la funzione di plot
  chi2Mat <- matrix(mat2D$chi2, nrow = length(xVals), ncol = length(yVals), byrow = FALSE)
  
  layoutObj <- plotPairFull(
    mat2D = chi2Mat,
    xGrid = xVals, yGrid = yVals,
    profX = profiles1D[[pX]], profY = profiles1D[[pY]],
    xlab = pair$xlab, ylab = pair$ylab,
    title = pair$title
  )
  
  plotsFull[[pair$title]] <- layoutObj
  print(layoutObj)
}

# Salvataggio opzionale
# for (i in seq_along(plotsFull)) {
#   ggsave(paste0("chi2_full_", names(plotsFull)[i], ".png"), 
#          plotsFull[[i]], width = 8, height = 6)
# }
