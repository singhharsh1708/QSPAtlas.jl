# Case 8 arbiter: near-exact population quantiles from the closed-form one-compartment oral
# solution with 4e6 draws of eta (MASS::mvrnorm) and eps, no ODE solver. Each tool's quantiles are
# compared with it: z = (tool - reference) / bootstrap SE of the tool's quantile.
# Usage: Rscript pop_reference.R [20k]   (uses out/compare_population[_20k].csv for the SEs)
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
args <- commandArgs(TRUE); tag <- if (length(args) && args[1] == "20k") "_20k" else ""
set.seed(99)
N <- 4e6
Om <- matrix(c(0.09, 0.024, 0.024, 0.04), 2)
eta <- MASS::mvrnorm(N, c(0, 0), Om)
CL <- 5 * exp(eta[, 1]); V <- 50 * exp(eta[, 2]); KA <- 1; k <- CL / V
ipred <- function(t) {
  s <- 0
  for (td in c(0, 24, 48)) if (t >= td) s <- s + 100 * KA / (V * (KA - k)) * (exp(-k * (t - td)) - exp(-KA * (t - td)))
  s
}
cmp <- read.csv(paste0("out/compare_population", tag, ".csv"))
q <- c(0.05, 0.5, 0.95)
ref <- list()
for (tt in sort(unique(cmp$time))) {
  ip <- ipred(tt); dv <- ip * (1 + 0.15 * rnorm(N))
  ref[[length(ref) + 1]] <- data.frame(time = tt, var = c(rep("IPRED", 3), rep("DV", 3)), q = rep(q, 2),
                                       reference = c(quantile(ip, q, type = 7), quantile(dv, q, type = 7)))
}
ref <- do.call(rbind, ref)
m <- merge(cmp, ref, by = c("time", "var", "q"))
m$z_qspatlas <- (m$qspatlas - m$reference) / m$se_qspatlas
m$z_mrgsolve <- (m$mrgsolve - m$reference) / m$se_mrgsolve
write.csv(m, paste0("out/reference_population", tag, ".csv"), row.names = FALSE)
for (tool in c("qspatlas", "mrgsolve")) {
  z <- m[[paste0("z_", tool)]]
  cat(sprintf("%-9s vs reference: mean z %+.2f, max |z| %.2f, |z|>2: %d of %d, negative %d, sum z^2 %.1f\n",
              tool, mean(z), max(abs(z)), sum(abs(z) > 2), length(z), sum(z < 0), sum(z^2)))
}
