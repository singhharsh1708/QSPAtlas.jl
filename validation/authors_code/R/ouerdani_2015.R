# Ouerdani 2015: the authors' NONMEM control streams, translated with nonmem2rx.
#   s004.txt: preclinical (CAKI-2 mice), s001.txt: clinical (RCC patients). ADVAN9, two states CELLS, CAPACITY.
# The streams carry INITIAL estimates; the THETAs are set to the atlas values (models/ouerdani_2015/model.jl,
# Table 1), so the comparison tests equations and inputs, not values. The data files ($DATA) are not published:
# nonmem2rx runs without them (no output-file validation possible), and the data items the code reads are
# supplied as covariates: DOSE (> 0 = treated), AUC, and the baseline VB (mice) / BASE (patients).
# Inputs as in the atlas fixes (protocol.jl FIXES): P0 = 150 mm^3 and AUC = 1140.8 (mice), P0 = 66.8 mm and
# AUC = 771.6 (patients); NONMEM time T in days from the start of treatment.
# Run from QSPAtlas.jl: Rscript validation/authors_code/R/ouerdani_2015.R
source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))
suppressMessages(library(nonmem2rx))

# nonmem2rx refuses ADVAN9 ("Differential Algebra Equations are not supported"). These streams have no
# $EQUILIBRIUM/$AES block, so ADVAN9 is used here as a plain stiff ODE solver; the only edit is
# "ADVAN9" -> "ADVAN13" in a temporary copy (tempdir(), outside the repository).
translate <- function(file) {
  src <- readLines(supplement("ouerdani_2015", file), warn = FALSE)
  stopifnot(!any(grepl("^\\$(AES|EQUILIBRIUM)", src)), sum(grepl("ADVAN9", src)) == 1)
  tmp <- file.path(tempdir(), file)
  writeLines(sub("ADVAN9", "ADVAN13", src), tmp)
  ui <- suppressWarnings(nonmem2rx(tmp, save = FALSE, validate = FALSE))
  cat("\n==== nonmem2rx translation of", file, "====\n"); print(ui)
  ui
}
pre <- translate("PSP4-4-660-s004.txt")
cli <- translate("PSP4-4-660-s001.txt")

# THETA index -> atlas value (Table 1). THETA(2) V0 is 0 FIX and unused; THETA(10), THETA(11) only scale etas.
theta_pre <- c(`1` = 0.166, `3` = 543.0, `4` = 0.0183, `5` = 0.007, `6` = 0.332, `7` = 0.251, `8` = 0.0, `9` = 0.196)
theta_cli <- c(`1` = 0.0021, `3` = 329.0, `4` = 0.0392, `5` = 0.0023, `6` = 0.142, `7` = 0.0032, `8` = 0.125, `9` = 0.0153)

set_thetas <- function(ui, th) {
  nm <- ui$iniDf$name[!is.na(ui$iniDf$ntheta)]
  idx <- ui$iniDf$ntheta[!is.na(ui$iniDf$ntheta)]
  vals <- setNames(th[as.character(idx[idx %in% as.integer(names(th))])], nm[idx %in% as.integer(names(th))])
  cat("theta names set:", paste(names(vals), signif(vals, 4), sep = "=", collapse = ", "), "\n")
  vals
}

simulate <- function(ui, th, cov, grid, case) {
  m <- rxode2::zeroRe(ui)                      # typical values: all etas 0 (and no residual error)
  ev <- data.frame(id = 1, time = grid, evid = 0, amt = 0, cov)
  s <- as.data.frame(rxSolve(m, params = set_thetas(ui, th), events = ev, atol = 1e-12, rtol = 1e-10,
                             maxsteps = 1e6, returnType = "data.frame"))
  st <- names(s)[tolower(names(s)) %in% c("cells", "capacity")]
  data.frame(case = case, time = rep(s$time, 2), var = rep(c("P", "K"), each = nrow(s)),
             value = c(s[[st[tolower(st) == "cells"]]], s[[st[tolower(st) == "capacity"]]]))
}

mg <- seq(0, 75, by = 0.25)
pg <- seq(0, 2500, by = 5)
res <- rbind(
  simulate(pre, theta_pre, data.frame(DOSE = 100, AUC = 1140.8, VB = 150), mg, "mouse_fig4"),
  simulate(pre, theta_pre, data.frame(DOSE = 10, AUC = 220.2, VB = 150), mg, "mouse_auc220"),
  simulate(pre, theta_pre, data.frame(DOSE = 30, AUC = 656.8, VB = 150), mg, "mouse_auc657"),
  simulate(pre, theta_pre, data.frame(DOSE = 0, AUC = 0, VB = 150), mg, "mouse_control_auc0"),
  simulate(cli, theta_cli, data.frame(DOSE = 800, AUC = 771.6, BASE = 66.8), pg, "patient_fig4"))
write_long(res, file.path(OUT, "rx_ouerdani_2015.csv"))
