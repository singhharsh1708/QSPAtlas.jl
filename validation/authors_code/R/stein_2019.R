# Stein 2019: the authors' MLXTRAN model is printed only inside the supplement PDF
# (PSP4-8-285-s001.pdf, section "Model MLXTRAN file"); there is no .mlxtran or Monolix project file, so
# monolix2rx has nothing to read. Translated BY HAND, line by line, into rxode2 (algebraic model, no ODE).
# The script checks that the PDF still contains the block it translates (text extraction with pdftotext).
#
# Two notes on the printed code:
#  * "R1 = R0*exp(rho* t1)" is broken over three lines in the PDF text layer; read as one expression.
#  * the third branch is printed "elseif t < tmax" (lower case) while the parameter is "Tmax". Mlxtran
#    is case sensitive, so as printed tmax is undefined; read as Tmax (a transcription slip in the PDF).
#    Without comedication t1 = t2 = Tmax, so this branch is never reached and the slip has no effect.
#
# Parameters: atlas values (models/stein_2019/model.jl TYPICAL = Table 1 fixed effects).
# Run from QSPAtlas.jl: Rscript validation/authors_code/R/stein_2019.R
source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

pdf <- supplement("stein_2019", "PSP4-8-285-s001.pdf")
txt <- tryCatch(system2("pdftotext", c(shQuote(pdf), "-"), stdout = TRUE), error = function(e) character())
if (length(txt)) {
  need <- c("Model MLXTRAN file", "CmaxAdjust = R2*exp(F1*F2*rho*(Tmax -t2))",
            "ylin = AA*exp(-alpha*(t-Tmax)) + BB*exp(-beta*(t-Tmax))", "elseif t < tmax")
  found <- sapply(need, function(s) any(grepl(s, txt, fixed = TRUE)))
  print(found)
  if (!all(found)) warning("the PDF text does not contain every translated line; check the translation")
}

stein <- rxode2({
  # first comedication -> (t1, F1), second -> (t2, F2); 99999 = not given (as in the authors' data set)
  if (TOCI1T <= STERSTT) {
    t1 <- min(TOCI1T, Tmax)
    F1c <- Ftoci
    t2 <- min(STERSTT, Tmax)
    F2c <- Fster
  } else {
    t1 <- min(STERSTT, Tmax)
    F1c <- Fster
    t2 <- min(TOCI1T, Tmax)
    F2c <- Ftoci
  }
  rho <- log(foldx) / Tmax
  R0 <- Cmax / exp(rho * Tmax)
  R1 <- R0 * exp(rho * t1)
  R2 <- R1 * exp(F1c * rho * (t2 - t1))
  CmaxAdjust <- R2 * exp(F1c * F2c * rho * (Tmax - t2))
  AA <- (1 - FB) * CmaxAdjust
  BB <- FB * CmaxAdjust
  if (t < 0) {
    ylin <- R0
  } else if (t < t1) {
    ylin <- R0 * exp(rho * t)
  } else if (t < t2) {
    ylin <- R1 * exp(F1c * rho * (t - t1))
  } else if (t < Tmax) {          # printed "tmax"
    ylin <- R2 * exp(F1c * F2c * rho * (t - t2))
  } else {
    ylin <- AA * exp(-alpha * (t - Tmax)) + BB * exp(-beta * (t - Tmax))
  }
  logy <- log(ylin)
})

theta <- c(foldx = 3900, Tmax = 9.3, Cmax = 24000, Ftoci = 1.2, Fster = 1.0, alpha = 0.16, FB = 0.0079, beta = 0.0032)
grid <- sort(unique(c(seq(0, 30, by = 0.1), seq(31, 380, by = 1))))
cases <- list(typical = c(99999, 99999), toci4_ster6 = c(4, 6), ster2_toci6 = c(6, 2), toci12 = c(12, 99999))
res <- do.call(rbind, lapply(names(cases), function(nm) {
  ev <- data.frame(id = 1, time = grid, evid = 0, TOCI1T = cases[[nm]][1], STERSTT = cases[[nm]][2])
  s <- as.data.frame(rxSolve(stein, theta, ev, returnType = "data.frame"))
  data.frame(case = nm, time = s$time, var = "C", value = exp(s$logy))
}))
write_long(res, file.path(OUT, "rx_stein_2019.csv"))
