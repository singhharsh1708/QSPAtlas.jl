# Dataset cases d1 to d7 in mrgsolve, from the same files as dataset_cases_qspatlas.jl
# (datasets/*.csv). Default record handling: the data set order is kept (recsort 1, the
# default, does not reorder data set records; it places ADDL doses before observations at the
# same time), nocb = TRUE unless the case says otherwise. Writes out/mrg_<case>.csv with one
# row per observation record (EVID 0), in record order.
# Run: SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk Rscript dataset_cases_mrgsolve.R
suppressMessages(library(mrgsolve))
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
dir.create("out", showWarnings = FALSE)

code <- '
$PARAM CL = 2, V1 = 10, Q = 3, V2 = 20, KA = 1.1,
       F1 = 1, LAG1 = 0, F2 = 1, LAG2 = 0,
       KIN = 10, KOUT = 0.1, IMAX = 0.8, IC50 = 2, WT = 70, RC = 0, DC = 0
$CMT DEPOT CENT PERIPH RESP
$MAIN
double CLi = CL * pow(WT / 70, 0.75);
double V1i = V1 * WT / 70;
F_DEPOT = F1; ALAG_DEPOT = LAG1;
F_CENT = F2;  ALAG_CENT = LAG2;
if(RC > 0) R_CENT = RC;
if(DC > 0) D_CENT = DC;
RESP_0 = KIN / KOUT;
$ODE
double CPi = CENT / V1i;
dxdt_DEPOT  = -KA * DEPOT;
dxdt_CENT   = KA * DEPOT - (CLi + Q) / V1i * CENT + Q / V2 * PERIPH;
dxdt_PERIPH = Q / V1i * CENT - Q / V2 * PERIPH;
dxdt_RESP   = KIN * (1 - IMAX * CPi / (IC50 + CPi)) - KOUT * RESP;
$TABLE capture CP = CENT / V1i;
'
mod0 <- mcode("qspxd", code, quiet = TRUE)
mod0 <- update(mod0, rtol = 1e-10, atol = 1e-12, maxsteps = 1e7, ss_rtol = 1e-10, ss_atol = 1e-12)

cases <- list(
  d1_order    = list(file = "d1_order.csv"),
  d2_cov_nocb = list(file = "d2_cov.csv", nocb = TRUE),
  d2_cov_locf = list(file = "d2_cov.csv", nocb = FALSE),
  d3_evid2    = list(file = "d3_evid2.csv"),
  d4_ss2      = list(file = "d4_ss2.csv"),
  d4_ss2_lin  = list(file = "d4_ss2.csv", p = list(IMAX = 0)),
  d5_ssinf    = list(file = "d5_ssinf.csv", p = list(RC = 30)),
  d6_rate     = list(file = "d6_rate.csv", p = list(F2 = 0.5, RC = 20, DC = 5)),
  d7_reset    = list(file = "d7_reset.csv")
)

for (nm in names(cases)) {
  cs <- cases[[nm]]
  df <- read.csv(file.path("datasets", cs$file), na.strings = c(".", ""))
  for (c in setdiff(names(df), "DV")) df[[c]][is.na(df[[c]])] <- 0
  mod <- if (is.null(cs$p)) mod0 else update(mod0, param = cs$p)
  nocb <- if (is.null(cs$nocb)) TRUE else cs$nocb
  res <- tryCatch(withCallingHandlers(
    as.data.frame(mrgsim(mod, data = df, carry_out = "EVID", nocb = nocb, ss_n = 20000)),
    warning = function(w) { message(nm, " WARNING: ", conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(err) { message(nm, " ERROR: ", conditionMessage(err)); NULL })
  if (is.null(res)) next
  res <- res[res$EVID == 0, ]
  stopifnot(nrow(res) == sum(df$EVID == 0))
  write.csv(res[, c("ID", "TIME", "DEPOT", "CENT", "PERIPH", "RESP", "CP")], file.path("out", paste0("mrg_", nm, ".csv")), row.names = FALSE)
  cat(nm, nrow(res), "observation records\n")
}
