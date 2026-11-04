# Case 8 in mrgsolve (same model as pop_qspatlas.jl). Writes out/pop_mrg.csv, out/pop_mrg_eta.csv,
# out/pop_mrg_zero.csv. Run with SDKROOT set (see cases_mrgsolve.R).
suppressMessages(library(mrgsolve))
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
code <- '
$PARAM TVCL = 5, TVV = 50, KA = 1
$CMT DEPOT CENT
$OMEGA @block
0.09 0.024 0.04
$SIGMA 0.0225
$MAIN
double CL = TVCL * exp(ETA(1));
double V  = TVV  * exp(ETA(2));
$ODE
dxdt_DEPOT = -KA * DEPOT;
dxdt_CENT  = KA * DEPOT - CL / V * CENT;
$TABLE
capture IPRED = CENT / V;
capture DV = IPRED * (1 + EPS(1));
capture ETA1 = ETA(1);
capture ETA2 = ETA(2);
'
mod <- mcode("qsppop", code, quiet = TRUE)
mod <- update(mod, rtol = 1e-10, atol = 1e-12, maxsteps = 1e7)
obst <- c(0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 24.5, 26, 30, 36, 48, 50, 54, 60, 72)
e <- ev(amt = 100, cmt = 1, ii = 24, addl = 2, ID = 1:5000)
set.seed(20261104)
out <- as.data.frame(mrgsim(mod, events = e, tgrid = obst, obsonly = TRUE, recsort = 3, end = -1))
write.csv(out[, c("ID", "time", "IPRED", "DV")], file.path("out", "pop_mrg.csv"), row.names = FALSE)
eta <- out[!duplicated(out$ID), c("ID", "ETA1", "ETA2")]
write.csv(eta, file.path("out", "pop_mrg_eta.csv"), row.names = FALSE)
# a second, independent run: 4 x 5000 subjects with other seeds
big <- do.call(rbind, lapply(1:4, function(k) {
  set.seed(1000 + k)
  o <- as.data.frame(mrgsim(mod, events = e, tgrid = obst, obsonly = TRUE, recsort = 3, end = -1))
  o$ID <- paste0(k, "_", o$ID); o[, c("ID", "time", "IPRED", "DV")]
}))
write.csv(big, file.path("out", "pop_mrg_20k.csv"), row.names = FALSE)
grid <- seq(0, 72, by = 0.25)
z <- as.data.frame(mrgsim(zero_re(mod), events = ev(amt = 100, cmt = 1, ii = 24, addl = 2), tgrid = grid, obsonly = TRUE, recsort = 3, end = -1))
d <- as.data.frame(mrgsim(update(mod, omega = dmat(0, 0), sigma = dmat(0)), events = ev(amt = 100, cmt = 1, ii = 24, addl = 2),
                          tgrid = grid, obsonly = TRUE, recsort = 3, end = -1))
write.csv(data.frame(time = z$time, ZERO = z$IPRED, DET = d$IPRED), file.path("out", "pop_mrg_zero.csv"), row.names = FALSE)
cat("rows", nrow(out), "subjects", nrow(eta), "\n")
