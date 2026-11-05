# Dose control with random effects and between-occasion variability in mrgsolve, the same cases
# as iov_cases_qspatlas.jl:
#   F (depot) = TVF exp(ETAF) WT / 70, ALAG (depot) = TVLAG exp(ETALAG), CL = TVCL exp(kappa),
#   kappa = K1 on occasion 1 and K2 on occasion 2 (OCC data item, NONMEM-style indicator code).
#  1. Exact case: ETAF, ETALAG, K1, K2 supplied in the data set (datasets/iov_exact.csv); default
#     record handling (nocb = TRUE, data set order kept). Writes out/mrg_iov_exact.csv.
#  2. Distribution case: out/iov_pop.csv (written by the Julia script), ETAs drawn by mrgsolve:
#     ETA(1) = ETAF, ETA(2) = ETALAG, ETA(3) and ETA(4) = kappa on occasions 1 and 2 (independent,
#     same variance: between-occasion variability). Writes out/mrg_iov_pop.csv.
# Run (after the Julia script): SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX26.5.sdk Rscript iov_cases_mrgsolve.R
suppressMessages(library(mrgsolve))
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
dir.create("out", showWarnings = FALSE)

code_exact <- '
$PARAM TVCL = 2, V = 10, KA = 1.1, TVF = 0.6, TVLAG = 0.75, WT = 70, OCC = 1,
       ETAF = 0, ETALAG = 0, K1 = 0, K2 = 0
$CMT DEPOT CENT
$MAIN
double KAP = (OCC == 1) ? K1 : K2;
double CLi = TVCL * exp(KAP);
F_DEPOT = TVF * exp(ETAF) * WT / 70;
ALAG_DEPOT = TVLAG * exp(ETALAG);
$ODE
dxdt_DEPOT = -KA * DEPOT;
dxdt_CENT  = KA * DEPOT - CLi / V * CENT;
$TABLE capture CP = CENT / V;
'
mod <- mcode("qspiov", code_exact, quiet = TRUE)
mod <- update(mod, rtol = 1e-12, atol = 1e-14, maxsteps = 1e8)
df <- read.csv(file.path("datasets", "iov_exact.csv"), na.strings = c(".", ""))
for (c in setdiff(names(df), "DV")) df[[c]][is.na(df[[c]])] <- 0
res <- as.data.frame(mrgsim(mod, data = df, carry_out = "EVID"))
res <- res[res$EVID == 0, ]
stopifnot(nrow(res) == sum(df$EVID == 0))
write.csv(res[, c("ID", "TIME", "DEPOT", "CENT", "CP")], file.path("out", "mrg_iov_exact.csv"), row.names = FALSE)

code_pop <- '
$PARAM TVCL = 2, V = 10, KA = 1.1, TVF = 0.6, TVLAG = 0.75, WT = 70, OCC = 1
$OMEGA 0.09 0.16 0.0625 0.0625
$CMT DEPOT CENT
$MAIN
double KAP = (OCC == 1) ? ETA(3) : ETA(4);
double CLi = TVCL * exp(KAP);
F_DEPOT = TVF * exp(ETA(1)) * WT / 70;
ALAG_DEPOT = TVLAG * exp(ETA(2));
$ODE
dxdt_DEPOT = -KA * DEPOT;
dxdt_CENT  = KA * DEPOT - CLi / V * CENT;
$TABLE capture CP = CENT / V;
'
modp <- mcode("qspiovpop", code_pop, quiet = TRUE)
modp <- update(modp, rtol = 1e-10, atol = 1e-12, maxsteps = 1e8)
dp <- read.csv(file.path("out", "iov_pop.csv"), na.strings = c(".", ""))
for (c in names(dp)) dp[[c]][is.na(dp[[c]])] <- 0
set.seed(20261105)
rp <- as.data.frame(mrgsim(modp, data = dp, obsonly = TRUE))
write.csv(rp[, c("ID", "TIME", "CP")], file.path("out", "mrg_iov_pop.csv"), row.names = FALSE)
cat("mrgsolve", as.character(packageVersion("mrgsolve")), ": exact", nrow(res), "rows, distribution", nrow(rp), "rows\n")
