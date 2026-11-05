# Hutchinson 2016: the authors' MLXTRAN structural model (PSP4-5-636-s004.txt).
# Route 1 (translator): monolix2rx, which reads Monolix PROJECT files; the supplement has only the structural
#   model, so a minimal wrapper project and a header-converted copy of the model go to tempdir() (see below).
# Route 2 (by hand): the same model written directly in rxode2, line by line from s004 (rxode2 rejects T as a
#   state name in hand-written code, so the states are Tum (= T) and Ves (= V)).
# Both are simulated for the atlas Figure 2 cases and Figure 6a; the comparison is in compare.R.
#
# Dosing as in the code: iv(adm=2, cmt=1) / iv(adm=4, cmt=2) put a dose into the amount compartments Bev / Van,
# which have ddt = 0 (no elimination); the window factor is N = Nmax*(Bev+Van) + 1 and the vessel kill is
# (delta_V_bev*Bev + delta_V_van*Van)*V. The dose records of the authors' data file are not published; one
# dose of amount 1 at day 38 is used (so the kill is delta_V and N = Nmax + 1 in the window).
# Like-for-like: Nmax(code) = N(atlas) - 1, so the window factor is the same number; the "samenum" cases use
# Nmax(code) = the atlas number to show the size of the +1.
# Parameters: atlas values (models/hutchinson_2016/protocol.jl fig2_base, DV, NM, fig6_base; T0 = 80 for Fig 6a
# as in the atlas fix T0_fig6_from_figure).
# Run from QSPAtlas.jl: Rscript validation/authors_code/R/hutchinson_2016.R
source(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), "common.R"))

mlx <- supplement("hutchinson_2016", "PSP4-5-636-s004.txt")
src <- readLines(mlx, warn = FALSE)
stopifnot(any(grepl("N = Nmax*(Bev+Van) +1", src, fixed = TRUE)), any(grepl("V_0 = K*X", src, fixed = TRUE)),
          any(grepl("t_0 = 38", src, fixed = TRUE)))

# ---- route 1: monolix2rx on a generated wrapper project ----------------------------------------------
# s004 is a Monolix 4 style structural model ("INPUT: parameter = {...}, regressor = X", no [LONGITUDINAL]
# header); monolix2rx reads only the current syntax, so a temporary copy (tempdir(), outside the repository)
# gets the header lines "[LONGITUDINAL] / input = {...} / X = {use=regressor}" in place of the INPUT: block.
# Everything from "PK:" to the end is passed unchanged. The wrapper project declares the parameters as fixed
# population values without variability (normal distribution, so a value of 0 is allowed).
pnames <- c("alpha_V", "alpha_T", "K", "Tnorm1", "Tnorm2", "Nmax", "delta_V_bev", "delta_V_van")
m2rx_model <- function() {
  tmpd <- file.path(tempdir(), "hutchinson_m2rx"); dir.create(tmpd, showWarnings = FALSE)
  i1 <- grep("^INPUT:", src); i2 <- grep("^PK:", src)
  stopifnot(length(i1) == 1, length(i2) == 1, i1 < i2,
            any(grepl("parameter = {alpha_V, alpha_T, K, Tnorm1, Tnorm2, Nmax, delta_V_bev, delta_V_van}", src[i1:i2], fixed = TRUE)))
  writeLines(c("[LONGITUDINAL]", paste0("input = {", paste(pnames, collapse = ", "), ", X}"), "X = {use=regressor}", "",
               src[i2:length(src)]), file.path(tmpd, "model.txt"))
  write.csv(data.frame(ID = 1, TIME = c(38, 39), Y = c(70, 80), AMT = c(1, NA), ADM = c(2, NA), X = 70),
            file.path(tmpd, "data.csv"), row.names = FALSE, quote = FALSE, na = ".")
  txt <- c("<DATAFILE>", "", "[FILEINFO]", "file = 'data.csv'", "delimiter = comma",
           "header = {ID, TIME, Y, AMT, ADM, X}", "", "[CONTENT]",
           "ID = {use=identifier}", "TIME = {use=time}", "Y = {use=observation, name=y, type=continuous}",
           "AMT = {use=amount}", "ADM = {use=administration}", "X = {use=regressor}", "",
           "<MODEL>", "", "[INDIVIDUAL]", paste0("input = {", paste0(pnames, "_pop", collapse = ", "), "}"), "",
           "DEFINITION:", paste0(pnames, " = {distribution=normal, typical=", pnames, "_pop, no-variability}"), "",
           "[LONGITUDINAL]", "input = {a}", "", "file = 'model.txt'", "",
           "DEFINITION:", "y = {distribution=normal, prediction=T, errorModel=constant(a)}", "",
           "<PARAMETER>", paste0(pnames, "_pop = {value=1, method=FIXED}"), "a = {value=1, method=MLE}", "")
  writeLines(txt, file.path(tmpd, "project.mlxtran"))
  suppressMessages(library(monolix2rx))
  ui <- monolix2rx(file.path(tmpd, "project.mlxtran"), update = FALSE)
  cat("\n==== monolix2rx translation of s004 ====\n"); print(ui)
  ui
}
translated <- m2rx_model()

# ---- cases (atlas protocol.jl: fig2_base, DV, NM; fig6_base) ----------------------------------------
fig2 <- c(alpha_T = 0.1, alpha_V = 0.09, K = 2, X = 70, Tnorm1 = 52, Tnorm2 = 62)
fig6 <- c(alpha_T = 0.109, alpha_V = 0.119, K = 1.14, X = 80, Tnorm1 = 53.4, Tnorm2 = 59.2)
cases <- list(fig2_a = c(fig2, dose = 0, delta = 0, Nmax = 0))
for (d in c(4, 6, 8, 10, 12)) cases[[sprintf("fig2_c_d%02d", d)]] <- c(fig2, dose = 1, delta = d / 100, Nmax = 0)
for (n in c(2, 3, 6, 11, 21)) {
  cases[[sprintf("fig2_e_N%d", n)]] <- c(fig2, dose = 1, delta = 0.1, Nmax = n - 1)          # like-for-like
  cases[[sprintf("fig2_e_N%d_samenum", n)]] <- c(fig2, dose = 1, delta = 0.1, Nmax = n)      # code with the same number
}
cases$fig6_a <- c(fig6, dose = 1, delta = 0.115, Nmax = 6.7 - 1)                             # vanucizumab, Table 1
cases$fig6_a_samenum <- c(fig6, dose = 1, delta = 0.115, Nmax = 6.7)

grid <- seq(38, 71, by = 0.05)
# Monolix "t_0 = 38": the initial values T_0 = X, V_0 = K*X, Bev_0 = Van_0 = 0 hold at day 38. rxode2 (and the
# monolix2rx translation, which drops t_0) starts at time 0, so a reset event (evid 3) at day 38 puts the
# states back to their initial values there; the dose record at day 38 follows the reset (checked below).
# rxode2 drops a reset that is the first record, hence the leading observation at time 0 (discarded).
events <- function(cs) {
  ev <- data.frame(id = 1, time = c(0, 38), evid = c(0, 3), amt = 0, cmt = "Van")
  if (cs[["dose"]] > 0) ev <- rbind(data.frame(id = 1, time = 38, evid = 1, amt = cs[["dose"]], cmt = "Van"), ev)
  ev <- rbind(ev, data.frame(id = 1, time = grid, evid = 0, amt = 0, cmt = "Van"))
  ev$X <- cs[["X"]]
  ev
}
# The dose goes to Van (adm=4, cmt=2) for every treated case; Bev and Van enter the equations symmetrically,
# so this equals a Bev dose with delta_V_bev = delta. The window uses absolute days, as in the code.
run_case <- function(model, nm, cs, m2rx) {
  par <- c(alpha_T = cs[["alpha_T"]], alpha_V = cs[["alpha_V"]], K = cs[["K"]], Tnorm1 = cs[["Tnorm1"]],
           Tnorm2 = cs[["Tnorm2"]], Nmax = cs[["Nmax"]], delta_V_bev = 0, delta_V_van = cs[["delta"]])
  if (m2rx) names(par) <- paste0(names(par), "_pop")
  s <- as.data.frame(rxSolve(model, params = par, events = events(cs), atol = 1e-12, rtol = 1e-10,
                             maxsteps = 1e6, returnType = "data.frame"))
  Tn <- if (m2rx) "T" else "Tum"; Vn <- if (m2rx) "V" else "Ves"
  s <- s[s$time >= 38, ]
  stopifnot(abs(s[[Tn]][1] - cs[["X"]]) < 1e-9, abs(s[[Vn]][1] - cs[["K"]] * cs[["X"]]) < 1e-9,
            abs(s$Van[1] - cs[["dose"]]) < 1e-12)
  rbind(data.frame(case = nm, time = s$time, var = "T", value = s[[Tn]]),
        data.frame(case = nm, time = s$time, var = "V", value = s[[Vn]]),
        data.frame(case = nm, time = s$time, var = "N", value = s$N))
}
# ---- route 2: by hand ------------------------------------------------------------------------------
hand <- rxode2({
  Tum(0) <- X
  Ves(0) <- K * X
  if (t > Tnorm1 & t < Tnorm2) {
    N <- Nmax * (Bev + Van) + 1
  } else {
    N <- 1
  }
  d/dt(Bev) <- 0
  d/dt(Van) <- 0
  d/dt(Tum) <- alpha_T * Tum * (1 - Tum / (N * Ves))
  d/dt(Ves) <- alpha_V * Tum^(2/3) * Ves^(1/3) - (delta_V_bev * Bev + delta_V_van * Van) * Ves
})
res <- rbind(
  do.call(rbind, lapply(names(cases), function(nm) run_case(hand, nm, cases[[nm]], FALSE))),
  do.call(rbind, lapply(names(cases), function(nm) run_case(rxode2::zeroRe(translated), paste0(nm, "_m2rx"), cases[[nm]], TRUE))))
write_long(res, file.path(OUT, "rx_hutchinson_2016.csv"))
