# Compare the authors'-code translations (out/rx_<id>.csv) with the atlas curves (out/atlas_<id>.csv).
# Writes out/summary_<id>.csv and summary.csv (the tables in README.md): per (case, variable) the largest
# absolute difference, the largest relative difference |rx - atlas| / max(|atlas|, floor) with
# floor = 1e-6 * max|atlas| over the curve, and the value scale of the curve.
# Run from QSPAtlas.jl: Rscript validation/authors_code/R/compare.R
f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
OUT <- file.path(dirname(normalizePath(f)), "..", "out")
SUM <- file.path(dirname(normalizePath(f)), "..", "summary")
dir.create(SUM, showWarnings = FALSE)

rd <- function(p) { d <- read.csv(p, stringsAsFactors = FALSE); d$key <- round(d$time, 6); d }

compare_pair <- function(rx, at, rx_case, at_case, vars, drop_times = numeric()) {
  out <- list()
  for (v in vars) {
    a <- at[at$case == at_case & at$var == v, ]
    r <- rx[rx$case == rx_case & rx$var == v, ]
    m <- merge(a[, c("key", "value")], r[, c("key", "value")], by = "key", suffixes = c(".atlas", ".rx"))
    m <- m[!(m$key %in% round(drop_times, 6)), ]
    if (!nrow(m)) stop("no common times for ", rx_case, " / ", at_case, " / ", v)
    scale <- max(abs(m$value.atlas))
    fl <- max(1e-6 * scale, 1e-300)
    d <- abs(m$value.rx - m$value.atlas)
    i <- which.max(d)
    out[[v]] <- data.frame(rx_case = rx_case, atlas_case = at_case, var = v, n = nrow(m),
      max_abs = signif(max(d), 4), at_time = m$key[i], max_rel = signif(max(d / pmax(abs(m$value.atlas), fl)), 4),
      max_atlas = signif(scale, 6))
  }
  do.call(rbind, out)
}

run <- function(id, pairs, vars, drop_times = numeric()) {
  rx <- rd(file.path(OUT, paste0("rx_", id, ".csv")))
  at <- rd(file.path(OUT, paste0("atlas_", id, ".csv")))
  s <- do.call(rbind, lapply(pairs, function(p) compare_pair(rx, at, p[1], p[2], vars, drop_times)))
  s <- cbind(model = id, s)
  write.csv(s, file.path(SUM, paste0("summary_", id, ".csv")), row.names = FALSE)
  print(s, row.names = FALSE)
  s
}

all <- list()
if (file.exists(file.path(OUT, "rx_ouerdani_2015.csv"))) {
  cs <- c("mouse_fig4", "mouse_auc220", "mouse_auc657", "patient_fig4", "mouse_control_auc0")
  all$ouerdani <- run("ouerdani_2015", lapply(cs, function(x) c(x, x)), c("P", "K"))
}
if (file.exists(file.path(OUT, "rx_hutchinson_2016.csv"))) {
  rx <- unique(read.csv(file.path(OUT, "rx_hutchinson_2016.csv"))$case)
  pairs <- lapply(rx, function(x) c(x, sub("_samenum$", "", sub("_m2rx$", "", x))))
  all$hutchinson <- run("hutchinson_2016", pairs, c("T", "V", "N"))
}
if (file.exists(file.path(OUT, "rx_stein_2019.csv"))) {
  cs <- c("typical", "toci4_ster6", "ster2_toci6", "toci12")
  all$stein <- run("stein_2019", lapply(cs, function(x) c(x, x)), "C")
}
if (file.exists(file.path(OUT, "rx_lindauer_2017.csv"))) {
  doses <- c("0p0", "0p1", "1p4", "5p0")
  pairs <- c(lapply(doses, function(d) c(paste0("mouse_", d, "mgkg"), paste0("mouse_", d, "mgkg_guard"))),
             lapply(doses, function(d) c(paste0("mouse_", d, "mgkg"), paste0("mouse_", d, "mgkg"))))
  # dose times (h) dropped: the value saved exactly at a bolus time is pre-dose in one tool, post-dose in the other
  all$lindauer <- run("lindauer_2017", pairs, c("CP", "ROb", "ROt", "TV"), drop_times = c(0, 120, 216, 312, 408))
}
s <- do.call(rbind, all)
write.csv(s, file.path(SUM, "summary.csv"), row.names = FALSE)
