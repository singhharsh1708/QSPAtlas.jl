# Lindauer 2017: the authors' NONMEM control stream for the mouse PK/RO/tumour model (PSP4-6-11-s006.txt,
# ADVAN13, 11 states), translated with nonmem2rx.
# s006 is a first-estimate stream: THETA(20)-(22) and (26)-(31) are initial values (lnTmulti 1.4, lnEMAXTP 4.6,
# lnEC50TP -6.6, KdegPD1 0.02, L0 0.12, L1 170, lnSLtg -10, lngamma 0.74) and Kon/Koff carry the in-vitro
# figures 340.2 and 0.10584. They are set to the atlas mouse values (models/lindauer_2017/model.jl mouse_raw(),
# Table 1 / Table S2), on the scale each THETA has in the code (log for 20-22, 30, 31; uM for EC50TP).
# The mouse data file is not used. Dosing as in the code: AMT in ug into C1 (DEFDOSE), converted to nmol by
# F1 = 1000/MW. Design: 20 g mouse, 0, 0.1, 1.4, 5 mg/kg on days 0, 5, 9, 13, 17 (Methods, test data set),
# 30 days. Outputs: the $ERROR quantities CP [mg/L], ROb, ROt [%], TV [mm^3].
# Run from QSPAtlas.jl: Rscript validation/authors_code/R/lindauer_2017.R
source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))
suppressMessages(library(nonmem2rx))

ui <- suppressWarnings(nonmem2rx(supplement("lindauer_2017", "PSP4-6-11-s006.txt"), save = FALSE, validate = FALSE))
cat("\n==== nonmem2rx translation of PSP4-6-11-s006.txt ====\n"); print(ui)

atlas_theta <- c(`20` = log(4.32), `21` = log(94.7), `22` = log(1.46 / 1000),    # Tmulti, EMAXTP, EC50TP [nM -> uM]
                 `23` = 340.0, `24` = 0.106, `26` = 0.0194,                       # KonPD1 [1e6/M/h], KoffPD1, KdegPD1 [1/h]
                 `27` = 0.113, `28` = 187.0, `29` = 170.0,                         # L0 [1/day], L1 [mm3/day], W0 [mm3]
                 `30` = log(1.98e-5), `31` = log(2.28))                            # SLtg [1/(% day)], gamma
th <- ui$iniDf[!is.na(ui$iniDf$ntheta), c("ntheta", "name", "est")]
par <- setNames(th$est, th$name)
for (k in names(atlas_theta)) par[th$name[th$ntheta == as.integer(k)]] <- atlas_theta[[k]]
cat("thetas changed from the stream:\n")
print(data.frame(theta = th$ntheta, name = th$name, stream = signif(th$est, 5), used = signif(par, 5))[th$est != par, ], row.names = FALSE)

m <- rxode2::zeroRe(ui)
MW <- 149000
grid <- seq(0, 720, by = 2)
dose_h <- c(0, 5, 9, 13, 17) * 24
res <- do.call(rbind, lapply(c(0, 0.1, 1.4, 5), function(mgkg) {
  ug <- mgkg * 0.02 * 1000                                       # mg/kg * 0.02 kg = mg -> ug
  obs <- data.frame(id = 1, time = grid, evid = 0, amt = 0, cmt = 1)
  ev <- if (mgkg > 0) rbind(data.frame(id = 1, time = dose_h, evid = 1, amt = ug, cmt = 1), obs) else obs
  ev <- ev[order(ev$time, -ev$evid), ]
  # data items read by $PK/$ERROR: DOSE (only copied to DOSECA), VARID (picks IPRED; 1 = PK); AMT is rxode2 amt
  ev$DOSE <- mgkg; ev$VARID <- 1
  s <- as.data.frame(rxSolve(m, params = par, events = ev, atol = 1e-14, rtol = 1e-10, maxsteps = 1e7,
                             returnType = "data.frame"))
  case <- paste0("mouse_", sub(".", "p", sprintf("%.1f", mgkg), fixed = TRUE), "mgkg")
  # nonmem2rx lower-cases the $ERROR names: CP -> cp, ROb -> rob, ROt -> rot, TV -> tv
  out <- c(CP = "cp", ROb = "rob", ROt = "rot", TV = "tv")
  do.call(rbind, lapply(names(out), function(v)
    data.frame(case = case, time = s$time, var = v, value = s[[out[[v]]]])))
}))
write_long(res, file.path(OUT, "rx_lindauer_2017.csv"))
