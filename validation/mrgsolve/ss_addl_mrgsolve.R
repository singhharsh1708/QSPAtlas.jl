# Case d8_ss_addl in mrgsolve: an SS=1 record after a record with ADDL doses still to come
# (datasets/d8_ss_addl.csv). Subject 1: bolus 100 every 12 h (ADDL 5), oral SS=1 dose at 30 h;
# subject 2: the same without SS (control); subject 3: infusions of 10 h every 12 h (ADDL 3),
# SS=1 bolus at 18 h while the second infusion runs. Same model and settings as
# dataset_cases_mrgsolve.R. Writes out/mrg_d8_ss_addl.csv; compare with dataset_compare.R.
# If ss_addl_qspatlas.jl has run, also compares mrgsolve with QSPAtlas run with the pending ADDL
# doses kept as separate dose records (out/jl_ss_addl_kept.csv).
# Run: SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk Rscript ss_addl_mrgsolve.R
suppressMessages(library(mrgsolve))
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
dir.create("out", showWarnings = FALSE)

code <- '
$PARAM CL = 2, V1 = 10, Q = 3, V2 = 20, KA = 1.1,
       KIN = 10, KOUT = 0.1, IMAX = 0.8, IC50 = 2, WT = 70
$CMT DEPOT CENT PERIPH RESP
$MAIN
double CLi = CL * pow(WT / 70, 0.75);
double V1i = V1 * WT / 70;
RESP_0 = KIN / KOUT;
$ODE
double CPi = CENT / V1i;
dxdt_DEPOT  = -KA * DEPOT;
dxdt_CENT   = KA * DEPOT - (CLi + Q) / V1i * CENT + Q / V2 * PERIPH;
dxdt_PERIPH = Q / V1i * CENT - Q / V2 * PERIPH;
dxdt_RESP   = KIN * (1 - IMAX * CPi / (IC50 + CPi)) - KOUT * RESP;
$TABLE capture CP = CENT / V1i;
'
mod <- mcode("qspxd8", code, quiet = TRUE)
mod <- update(mod, rtol = 1e-10, atol = 1e-12, maxsteps = 1e7, ss_rtol = 1e-10, ss_atol = 1e-12)

df <- read.csv(file.path("datasets", "d8_ss_addl.csv"), na.strings = c(".", ""))
for (c in setdiff(names(df), "DV")) df[[c]][is.na(df[[c]])] <- 0
res <- as.data.frame(mrgsim(mod, data = df, carry_out = "EVID", ss_n = 20000))
res <- res[res$EVID == 0, ]
stopifnot(nrow(res) == sum(df$EVID == 0))
write.csv(res[, c("ID", "TIME", "DEPOT", "CENT", "PERIPH", "RESP", "CP")], file.path("out", "mrg_d8_ss_addl.csv"), row.names = FALSE)
print(res[, c("ID", "TIME", "CENT", "RESP")], digits = 10, row.names = FALSE)

fk <- file.path("out", "jl_ss_addl_kept.csv")
if (file.exists(fk)) {
  k <- read.csv(fk)
  cols <- c("DEPOT", "CENT", "PERIPH", "RESP", "CP")
  B <- as.matrix(res[, cols]); A <- as.matrix(k[, cols])
  stopifnot(nrow(A) == nrow(B), all(k$ID == res$ID), all(abs(k$TIME - res$TIME) < 1e-12))
  den <- pmax(abs(B), matrix(1e-6 * apply(abs(B), 2, max), nrow(B), length(cols), byrow = TRUE))
  cat("mrgsolve vs QSPAtlas with the ADDL doses kept as dose records: largest relative difference",
      format(max(abs(A - B) / den), digits = 3), "\n")
}
