# Shared helpers for the authors'-code translations. Sourced by the per-model scripts.
# The authors' files are read from $QSPATLAS_PAPERS/<id>/supplement/ (default: ../papers next to the repo);
# they are parsed by the translators, never copied into the repository and never run as given.
suppressMessages(library(rxode2))

script_dir <- function() {
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
HERE <- script_dir()                                   # validation/authors_code/R
AC_DIR <- normalizePath(file.path(HERE, ".."))        # validation/authors_code
REPO <- normalizePath(file.path(AC_DIR, "..", ".."))   # QSPAtlas.jl
PAPERS <- Sys.getenv("QSPATLAS_PAPERS", file.path(REPO, "..", "papers"))
OUT <- file.path(AC_DIR, "out")
dir.create(OUT, showWarnings = FALSE)

supplement <- function(id, file) {
  p <- file.path(PAPERS, id, "supplement", file)
  if (!file.exists(p)) stop("authors' file not found: ", p, " (set QSPATLAS_PAPERS)")
  normalizePath(p)
}

# long format: case, time, var, value
write_long <- function(df, path) {
  write.csv(df, path, row.names = FALSE, quote = FALSE)
  cat("wrote", path, "(", nrow(df), "rows )\n")
}

to_long <- function(sol, case, vars) {
  sol <- as.data.frame(sol)
  do.call(rbind, lapply(vars, function(v)
    data.frame(case = case, time = sol$time, var = v, value = sol[[v]])))
}

SOLVE_TOL <- list(atol = 1e-12, rtol = 1e-10, maxsteps = 1e6)
