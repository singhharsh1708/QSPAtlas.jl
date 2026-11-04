# Case 8: compare the population distributions of QSPAtlas and mrgsolve (different RNGs).
# Per time and quantile (type 7, as QSPAtlas `percentiles`): difference, bootstrap SE of each
# tool's quantile (B = 2000 resamples of subjects), z = diff / sqrt(se_jl^2 + se_mrg^2), and a
# 95% bootstrap interval of the difference. Writes out/compare_population.csv.
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
set.seed(1)
args <- commandArgs(TRUE)
tag <- if (length(args) && args[1] == "20k") "_20k" else ""
a <- read.csv(paste0("out/pop_jl", tag, ".csv")); b <- read.csv(paste0("out/pop_mrg", tag, ".csv"))
q <- c(0.05, 0.5, 0.95); B <- if (tag == "") 2000 else 500
rows <- list()
for (tt in sort(unique(a$time))) for (v in c("IPRED", "DV")) {
  x <- a[a$time == tt, v]; y <- b[b$time == tt, v]
  qx <- quantile(x, q, type = 7); qy <- quantile(y, q, type = 7)
  bx <- replicate(B, quantile(sample(x, replace = TRUE), q, type = 7))
  by <- replicate(B, quantile(sample(y, replace = TRUE), q, type = 7))
  sx <- apply(bx, 1, sd); sy <- apply(by, 1, sd)
  ci <- apply(bx - by, 1, quantile, c(0.025, 0.975))
  rows[[length(rows) + 1]] <- data.frame(time = tt, var = v, q = q, qspatlas = qx, mrgsolve = qy,
    diff = qx - qy, rel_diff = (qx - qy) / qy, se_qspatlas = sx, se_mrgsolve = sy,
    z = (qx - qy) / sqrt(sx^2 + sy^2), diff_lo = ci[1, ] , diff_hi = ci[2, ])
}
tab <- do.call(rbind, rows)
write.csv(tab, paste0("out/compare_population", tag, ".csv"), row.names = FALSE)
cat("quantile comparisons:", nrow(tab), "\n")
cat("max |z| =", round(max(abs(tab$z)), 2), " count |z| > 2:", sum(abs(tab$z) > 2), "(expected ~", round(0.0455 * nrow(tab), 1),
    ") count |z| > 3:", sum(abs(tab$z) > 3), "\n")
cat("max |rel diff| =", signif(max(abs(tab$rel_diff)), 3), "\n")
cat("chi-square sum z^2 =", round(sum(tab$z^2), 1), "on", nrow(tab), "(correlated across rows, indicative only)\n")
print(head(tab[order(-abs(tab$z)), c("time", "var", "q", "qspatlas", "mrgsolve", "rel_diff", "z")], 8), digits = 4, row.names = FALSE)
cat("sign of z: negative", sum(tab$z < 0), "positive", sum(tab$z > 0), "\n")
if (tag != "") quit(save = "no")
# random effects: sample covariance
ea <- read.csv("out/pop_jl_eta.csv"); eb <- read.csv("out/pop_mrg_eta.csv")
cat("\neta covariance QSPAtlas:\n"); print(round(cov(ea[, c("ETA1", "ETA2")]), 4))
cat("eta covariance mrgsolve:\n"); print(round(cov(eb[, c("ETA1", "ETA2")]), 4))
cat("target: 0.09 0.024 0.04; correlation QSPAtlas", round(cor(ea$ETA1, ea$ETA2), 3), "mrgsolve", round(cor(eb$ETA1, eb$ETA2), 3), "\n")
# residual: (DV/IPRED - 1) sd
ra <- a$DV / a$IPRED - 1; rb <- b$DV / b$IPRED - 1
cat("residual sd (DV/IPRED - 1): QSPAtlas", round(sd(ra), 4), "mrgsolve", round(sd(rb), 4), "target 0.15\n")
# zero-RE prediction
za <- read.csv("out/pop_jl_zero.csv"); zb <- read.csv("out/pop_mrg_zero.csv")
cat("\nzero-RE vs no-random-effect model: QSPAtlas max|diff| =", max(abs(za$ZERO - za$DET)),
    "; mrgsolve max|diff| =", max(abs(zb$ZERO - zb$DET)), "\n")
cat("zero-RE QSPAtlas vs mrgsolve: max|diff| =", signif(max(abs(za$ZERO - zb$ZERO)), 3),
    "max rel =", signif(max(abs(za$ZERO - zb$ZERO) / pmax(abs(zb$ZERO), 1e-12))[1], 3), "\n")
