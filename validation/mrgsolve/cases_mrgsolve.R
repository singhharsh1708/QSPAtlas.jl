# Deterministic cases 1-7 in mrgsolve. Writes out/mrg_<case>.csv (time, DEPOT, CENT, PERIPH, RESP).
# Run: SDKROOT=.../MacOSX26.5.sdk Rscript cases_mrgsolve.R   (the default SDK 27 breaks the linker)
suppressMessages(library(mrgsolve))
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
dir.create("out", showWarnings = FALSE)

code <- '
$PARAM CL = 2, V1 = 10, Q = 3, V2 = 20, KA = 1.1,
       F1 = 1, LAG1 = 0, F2 = 1, LAG2 = 0, DUR2 = 0,
       KIN = 10, KOUT = 0.1, IMAX = 0.8, IC50 = 2
$CMT DEPOT CENT PERIPH RESP
$MAIN
F_DEPOT = F1; ALAG_DEPOT = LAG1;
F_CENT = F2;  ALAG_CENT = LAG2;
if(DUR2 > 0) D_CENT = DUR2;
RESP_0 = KIN / KOUT;
$ODE
double CPi = CENT / V1;
dxdt_DEPOT  = -KA * DEPOT;
dxdt_CENT   = KA * DEPOT - (CL + Q) / V1 * CENT + Q / V2 * PERIPH;
dxdt_PERIPH = Q / V1 * CENT - Q / V2 * PERIPH;
dxdt_RESP   = KIN * (1 - IMAX * CPi / (IC50 + CPi)) - KOUT * RESP;
$TABLE capture CP = CENT / V1;
'
mod0 <- mcode("qspx", code, quiet = TRUE)
mod0 <- update(mod0, rtol = 1e-10, atol = 1e-12, maxsteps = 1e7, ss_rtol = 1e-10, ss_atol = 1e-12)

ONE <- list(Q = 0)            # one-compartment: periphery never filled
DEP <- 1; CEN <- 2

cases <- list(
  c1   = list(p = ONE, tend = 96,  e = ev(amt = 100, cmt = CEN, ii = 12, addl = 5)),
  c2a  = list(p = list(), tend = 192, e = ev(amt = 100, cmt = CEN, rate = 100, ii = 12, addl = 13)),
  c2b  = list(p = list(), tend = 192, e = ev(amt = 100, cmt = CEN, rate = 100, ii = 12, addl = 13) +
                                         ev(amt = 50, cmt = CEN, rate = 10, time = 0.5, ii = 12, addl = 13)),
  c2c  = list(p = list(), tend = 192, e = ev(amt = 150, cmt = CEN, rate = 10, ii = 12, addl = 13)),
  c3a  = list(p = c(ONE, F1 = 0.6, LAG1 = 0.75), tend = 168, e = ev(amt = 100, cmt = DEP, ii = 24, addl = 6)),
  c3b_addl = list(p = c(ONE, LAG1 = 13), tend = 110, e = ev(amt = 100, cmt = DEP, ii = 12, addl = 6)),   # crashes mrgsolve 2.0.1
  c4a  = list(p = c(ONE, F2 = 0.5), tend = 72, e = ev(amt = 100, cmt = CEN, rate = 20, ii = 24, addl = 2)),
  c4b  = list(p = c(ONE, F2 = 0.5, DUR2 = 5), tend = 72, e = ev(amt = 100, cmt = CEN, rate = -2, ii = 24, addl = 2)),
  c4c  = list(p = c(ONE, F2 = 0.5, LAG2 = 1.5), tend = 72, e = ev(amt = 100, cmt = CEN, rate = 20, ii = 24, addl = 2)),
  c5a  = list(p = ONE, tend = 48, e = ev(amt = 100, cmt = CEN, ii = 12, addl = 3, ss = 1)),
  c5b  = list(p = c(ONE, F1 = 0.6, LAG1 = 0.75), tend = 72, e = ev(amt = 100, cmt = DEP, ii = 24, addl = 2, ss = 1)),
  c5b2 = list(p = c(ONE, LAG1 = 13), tend = 48, e = ev(amt = 100, cmt = DEP, ii = 12, addl = 3, ss = 1)),
  c5c  = list(p = list(), tend = 48, e = ev(amt = 100, cmt = CEN, rate = 20, ii = 12, addl = 3, ss = 1)),
  c5d  = list(p = list(), tend = 48, e = ev(amt = 180, cmt = CEN, rate = 10, ii = 12, addl = 3, ss = 1)),
  c5e  = list(p = ONE, tend = 48, e = ev(amt = 100, cmt = CEN, ii = 12, addl = 3, ss = 1)),
  c5f  = list(p = c(ONE, KIN = 0.05, KOUT = 0.0005), tend = 48, ss_n = 20000,
              e = ev(amt = 100, cmt = CEN, ii = 12, addl = 3, ss = 1)),
  c5f500 = list(p = c(ONE, KIN = 0.05, KOUT = 0.0005), tend = 48, ss_n = 500,
              e = ev(amt = 100, cmt = CEN, ii = 12, addl = 3, ss = 1)),
  c5g  = list(p = ONE, tend = 60, e = ev(amt = 500, cmt = CEN) + ev(amt = 100, cmt = CEN, time = 24, ii = 12, addl = 2, ss = 1)),
  c6a  = list(p = ONE, tend = 24, e = ev(amt = 100, cmt = CEN, rate = 10) + ev(time = 4, evid = 3, cmt = CEN, amt = 0) +
                                      ev(amt = 50, cmt = CEN, time = 8)),
  c6b  = list(p = ONE, tend = 24, e = ev(amt = 100, cmt = CEN, rate = 10) + ev(time = 4, evid = 4, cmt = CEN, amt = 50)),
  c6c  = list(p = ONE, tend = 24, e = ev(amt = 100, cmt = CEN, rate = 10) + ev(time = 4, evid = 3, cmt = CEN, amt = 0) +
                                      ev(amt = 50, cmt = CEN, rate = 5, time = 6)),
  c6d  = list(p = ONE, tend = 24, e = ev(amt = 100, cmt = CEN) + ev(time = 12, evid = 3, cmt = CEN, amt = 0)),
  c6e  = list(p = c(ONE, LAG1 = 2), tend = 12, e = ev(amt = 100, cmt = DEP, time = 3) + ev(time = 4, evid = 3, cmt = DEP, amt = 0)),
  c6f  = list(p = ONE, tend = 24, e = ev(amt = 100, cmt = CEN, ii = 6, addl = 3) + ev(time = 8, evid = 3, cmt = CEN, amt = 0)),
  # workaround for the addl + lag > ii crash: the same doses written out one by one
  c3b  = list(p = c(ONE, LAG1 = 13), tend = 110, e = ev(amt = 100, cmt = DEP, time = 12 * (0:6))),
  c7a  = list(p = ONE, tend = 24, e = ev(amt = 100, cmt = CEN) + ev(amt = 50, cmt = CEN) +
                                      ev(amt = 100, cmt = CEN, time = 6) + ev(amt = 50, cmt = CEN, time = 6)),
  c7b  = list(p = ONE, tend = 24, e = ev(amt = 100, cmt = CEN, rate = 50) + ev(amt = 30, cmt = CEN, time = 2) +
                                      ev(amt = 60, cmt = CEN, rate = 20, time = 2))
)

for (nm in names(cases)) {
  cs <- cases[[nm]]
  mod <- update(mod0, param = cs$p)
  ss_n <- if (is.null(cs$ss_n)) 500 else cs$ss_n
  res <- tryCatch(withCallingHandlers(
    as.data.frame(mrgsim(mod, events = cs$e, tgrid = seq(0, cs$tend, by = 0.25), recsort = 3,
                         obsonly = TRUE, ss_n = ss_n)),
    warning = function(w) { message(nm, " WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(err) { message(nm, " ERROR: ", conditionMessage(err)); NULL })
  if (is.null(res)) next
  write.csv(res[, c("time", "DEPOT", "CENT", "PERIPH", "RESP", "CP")], file.path("out", paste0("mrg_", nm, ".csv")), row.names = FALSE)
  cat(nm, nrow(res), "rows\n")
}

# c5b2: SS=1 with lag 13 > ii 12 crashes mrgsolve; reference = 200 explicit doses q12h, then the
# window starting at the 201st dose (t = 2400), shifted to t = 0.
mod <- update(mod0, param = c(ONE, LAG1 = 13))
res <- as.data.frame(mrgsim(mod, events = ev(amt = 100, cmt = DEP, time = 12 * (0:203)),
                            tgrid = 2400 + seq(0, 48, by = 0.25), recsort = 3, obsonly = TRUE))
res$time <- res$time - 2400
write.csv(res[, c("time", "DEPOT", "CENT", "PERIPH", "RESP", "CP")], file.path("out", "mrg_c5b2.csv"), row.names = FALSE)
cat("c5b2 (200 explicit doses) ", nrow(res), "rows\n")

# c5h (mrgsolve only): where does the SS=1 iteration start? Slow turnover, ss_n = 500 (not
# converged), SS dose at t = 24 with and without a large bolus at t = 0 before it.
mod <- update(mod0, param = c(ONE, KIN = 0.05, KOUT = 0.0005))
ss24 <- ev(amt = 100, cmt = CEN, time = 24, ii = 12, ss = 1)
r1 <- suppressWarnings(as.data.frame(mrgsim(mod, events = ss24, tgrid = 24, recsort = 3, obsonly = TRUE, end = -1, ss_n = 500)))
r2 <- suppressWarnings(as.data.frame(mrgsim(mod, events = ev(amt = 3000, cmt = CEN) + ss24, tgrid = c(23.99, 24), recsort = 3, obsonly = TRUE, end = -1, ss_n = 500)))
r3 <- suppressWarnings(as.data.frame(mrgsim(mod, events = ss24, tgrid = 24, recsort = 3, obsonly = TRUE, end = -1, ss_n = 20000)))
cat(sprintf("c5h RESP at T=24 after SS: no prior dose %.6f; after 3000 bolus at t=0 (RESP just before T %.4f) %.6f; converged (ss_n 20000) %.6f\n",
            r1$RESP[1], r2$RESP[1], r2$RESP[2], r3$RESP[1]))
