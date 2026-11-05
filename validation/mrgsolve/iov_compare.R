# Compare the IOV / dose-control cases (iov_cases_qspatlas.jl, iov_cases_mrgsolve.R).
#  1. Exact case: largest absolute and relative differences per subject (one row per observation
#     record, same order).
#  2. Distribution case: 5th, 50th and 95th percentiles of CP per time in both tools, with the
#     difference in units of its bootstrap standard error (subjects resampled, 400 draws).
# Writes out/compare_iov.csv.
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
a <- read.csv(file.path("out", "jl_iov_exact.csv"))
b <- read.csv(file.path("out", "mrg_iov_exact.csv"))
stopifnot(nrow(a) == nrow(b), all(a$ID == b$ID), all(abs(a$TIME - b$TIME) < 1e-12))
cols <- c("DEPOT", "CENT", "CP")
rows <- list()
for (id in unique(a$ID)) {
  i <- a$ID == id
  d <- abs(as.matrix(a[i, cols]) - as.matrix(b[i, cols]))
  scale <- apply(abs(as.matrix(b[i, cols])), 2, max)
  den <- pmax(abs(as.matrix(b[i, cols])), matrix(1e-6 * scale, sum(i), length(cols), byrow = TRUE))
  rows[[length(rows) + 1]] <- data.frame(case = "exact", part = paste("ID", id), n = sum(i), maxabs = max(d), maxrel = max(d / den), z = NA)
}
ex <- do.call(rbind, rows)
print(ex, digits = 3, row.names = FALSE)

pa <- read.csv(file.path("out", "jl_iov_pop.csv"))
pb <- read.csv(file.path("out", "mrg_iov_pop.csv"))
qs <- c(0.05, 0.5, 0.95)
set.seed(1)
B <- 400
out <- list()
for (tt in sort(unique(pa$TIME))) {
  xa <- pa$CP[abs(pa$TIME - tt) < 1e-9]
  xb <- pb$CP[abs(pb$TIME - tt) < 1e-9]
  qa <- quantile(xa, qs, type = 7); qb <- quantile(xb, qs, type = 7)
  ba <- replicate(B, quantile(sample(xa, replace = TRUE), qs, type = 7))
  bb <- replicate(B, quantile(sample(xb, replace = TRUE), qs, type = 7))
  se <- sqrt(apply(ba, 1, var) + apply(bb, 1, var))
  # quantiles that are exactly 0 in both tools (before the lagged dose enters): no difference
  for (k in seq_along(qs)) {
    out[[length(out) + 1]] <- data.frame(case = "distribution", part = sprintf("t=%g q=%g", tt, qs[k]), n = length(xa),
                                         maxabs = abs(qa[k] - qb[k]), maxrel = ifelse(qa[k] == qb[k], 0, abs(qa[k] - qb[k]) / abs(qb[k])),
                                         z = ifelse(qa[k] == qb[k], 0, (qa[k] - qb[k]) / se[k]))
  }
}
di <- do.call(rbind, out)
cat("\ndistribution: ", nrow(di), " quantiles, largest |z| ", round(max(abs(di$z)), 2),
    ", largest relative difference ", signif(max(di$maxrel), 3), ", |z| > 2: ", sum(abs(di$z) > 2), "\n", sep = "")
print(di[order(-abs(di$z)), ][1:5, ], digits = 3, row.names = FALSE)
write.csv(rbind(ex, di), file.path("out", "compare_iov.csv"), row.names = FALSE)
