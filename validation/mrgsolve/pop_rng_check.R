# Compare the QSPAtlas 10^6-draw sampling check (pop_rng_check.jl) with the 4e6-draw reference
# (pop_reference.R). SE scaled from the 5000-subject bootstrap SE: se * sqrt(5000 / n).
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
a <- read.csv("out/rng_check_jl.csv"); r <- read.csv("out/reference_population.csv")
m <- merge(a, r[, c("time", "var", "q", "reference", "se_qspatlas", "se_mrgsolve")], by = c("time", "var", "q"))
se <- sqrt(((m$se_qspatlas + m$se_mrgsolve) / 2)^2 * (5000 / 1e6 + 5000 / 4e6))
m$z <- (m$value - m$reference) / se
write.csv(m, "out/rng_check_compare.csv", row.names = FALSE)
cat(sprintf("QSPAtlas 1e6 draws vs 4e6 reference: max |rel diff| %.2e, mean z %+.2f, max |z| %.2f, |z|>2: %d of %d, sum z^2 %.1f\n",
  max(abs(m$value / m$reference - 1)), mean(m$z), max(abs(m$z)), sum(abs(m$z) > 2), nrow(m), sum(m$z^2)))
