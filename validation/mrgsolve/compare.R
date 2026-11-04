# Compare out/jl_<case>.csv with out/mrg_<case>.csv. Writes out/compare_deterministic.csv and
# prints, for cases that differ, the first rows where they do.
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
cols <- c("DEPOT", "CENT", "PERIPH", "RESP")
jl <- sub("^jl_(.*)\\.csv$", "\\1", list.files("out", "^jl_c.*\\.csv$"))
mg <- sub("^mrg_(.*)\\.csv$", "\\1", list.files("out", "^mrg_c.*\\.csv$"))
cases <- sort(union(jl, mg))
rows <- list()
for (cs in cases) {
  fj <- file.path("out", paste0("jl_", cs, ".csv")); fm <- file.path("out", paste0("mrg_", cs, ".csv"))
  if (cs == "c5f500") fj <- file.path("out", "jl_c5f.csv")      # mrgsolve default ss_n vs converged
  if (!file.exists(fj) || !file.exists(fm)) {
    rows[[cs]] <- data.frame(case = cs, maxabs = NA, maxrel = NA, worst_col = NA, worst_time = NA,
                             note = paste("missing", if (!file.exists(fj)) "QSPAtlas" else "mrgsolve"))
    next
  }
  a <- read.csv(fj); b <- read.csv(fm)
  stopifnot(nrow(a) == nrow(b), all(abs(a$time - b$time) < 1e-12))
  d <- abs(as.matrix(a[, cols]) - as.matrix(b[, cols]))
  scale <- apply(abs(as.matrix(b[, cols])), 2, max)
  rel <- d / pmax(abs(as.matrix(b[, cols])), matrix(1e-6 * pmax(scale, 1e-300), nrow(d), length(cols), byrow = TRUE))
  w <- which(rel == max(rel), arr.ind = TRUE)[1, ]
  rows[[cs]] <- data.frame(case = cs, maxabs = max(d), maxrel = max(rel), worst_col = cols[w[2]],
                           worst_time = a$time[w[1]], note = "")
  if (max(rel) > 1e-6) {
    bad <- which(apply(rel, 1, max) > 1e-6)
    cat("\n==", cs, ": ", length(bad), "of", nrow(a), "grid points differ (rel > 1e-6); first ones:\n")
    i <- head(bad, 8)
    print(data.frame(time = a$time[i], jl_CENT = a$CENT[i], mrg_CENT = b$CENT[i], jl_DEPOT = a$DEPOT[i],
                     mrg_DEPOT = b$DEPOT[i], jl_RESP = a$RESP[i], mrg_RESP = b$RESP[i]), digits = 10)
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path("out", "compare_deterministic.csv"), row.names = FALSE)
print(tab, digits = 3, row.names = FALSE)
