# Compare out/jl_<case>.csv with out/mrg_<case>.csv for the dataset cases (one row per observation
# record, same order). Writes out/compare_datasets.csv. Case d5_ssinf is compared in two parts,
# because mrgsolve and NONMEM differ there (see the validation page): at the steady-state record
# time, and QSPAtlas subject 2 (steady state plus an explicit continuing infusion, the NONMEM
# way to keep it running) against mrgsolve subject 1 (mrgsolve keeps the infusion running).
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
cols <- c("DEPOT", "CENT", "PERIPH", "RESP", "CP")
relmax <- function(a, b, cols) {
  d <- abs(as.matrix(a[, cols]) - as.matrix(b[, cols]))
  scale <- apply(abs(as.matrix(b[, cols])), 2, max)
  den <- pmax(abs(as.matrix(b[, cols])), matrix(1e-6 * pmax(scale, 1e-300), nrow(d), length(cols), byrow = TRUE))
  list(maxabs = max(d), maxrel = max(d / den), rel = d / den)
}
cases <- sub("^jl_(.*)\\.csv$", "\\1", list.files("out", "^jl_d.*\\.csv$"))
rows <- list()
for (cs in cases) {
  a <- read.csv(file.path("out", paste0("jl_", cs, ".csv")))
  fm <- file.path("out", paste0("mrg_", cs, ".csv"))
  if (!file.exists(fm)) { rows[[cs]] <- data.frame(case = cs, part = "", n = nrow(a), maxabs = NA, maxrel = NA, note = "no mrgsolve output"); next }
  b <- read.csv(fm)
  stopifnot(nrow(a) == nrow(b), all(a$ID == b$ID), all(abs(a$TIME - b$TIME) < 1e-12))
  all <- seq_len(nrow(a))
  parts <- list(all = list(i = all, j = all, cols = cols))
  if (startsWith(cs, "d4_ss2")) parts$PK <- list(i = all, j = all, cols = c("DEPOT", "CENT", "PERIPH", "CP"))
  if (cs == "d2_cov_locf") parts$states <- list(i = all, j = all, cols = c("DEPOT", "CENT", "PERIPH", "RESP"))
  for (pn in names(parts)) {
    p <- parts[[pn]]
    r <- relmax(a[p$i, ], b[p$j, ], p$cols)
    rows[[paste(cs, pn)]] <- data.frame(case = cs, part = pn, n = length(p$i), maxabs = r$maxabs, maxrel = r$maxrel, note = "")
    if (r$maxrel > 1e-6 && pn == "all") {
      bad <- which(apply(r$rel, 1, max) > 1e-6)
      cat("\n==", cs, pn, ":", length(bad), "of", length(p$i), "records differ (rel > 1e-6)\n")
      print(cbind(a[p$i, ][bad, c("ID", "TIME", "CENT", "RESP", "CP")], mrg_CENT = b[p$j, ][bad, "CENT"],
                  mrg_RESP = b[p$j, ][bad, "RESP"], mrg_CP = b[p$j, ][bad, "CP"]), digits = 10, row.names = FALSE)
    }
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path("out", "compare_datasets.csv"), row.names = FALSE)
print(tab, digits = 3, row.names = FALSE)
