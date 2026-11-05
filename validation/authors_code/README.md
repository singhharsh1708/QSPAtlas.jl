# Cross-check against the authors' own code

The atlas models (`models/<id>/`) were written from the papers. This folder checks four of them a third way:
the NONMEM and Monolix code that the authors published in their supplements is translated to
[rxode2](https://nlmixr2.github.io/rxode2/) (R) with an independent toolchain, run with the parameter values and
inputs the atlas uses, and compared point by point with the atlas curves. The parameter values are set to the
atlas values on purpose (the supplement code often carries initial estimates): the check is on equations,
structure, dosing and units, not on the values.

The authors' files are not in this repository (the supplements are not ours to redistribute). The scripts read
them from `$QSPATLAS_PAPERS/<id>/supplement/` (default `../papers` next to the repository), parse them with the
translators and never run them as given. Where a translator needed an edit, the edit is made on a temporary copy
in `tempdir()` and is stated below.

## Files

| file | what it does |
|---|---|
| `julia/atlas_curves.jl` | atlas curves: loads each `models/<id>/model.jl` and `protocol.jl`, runs the atlas functions with the atlas fixes (`FIXES`), writes `out/atlas_<id>.csv` |
| `R/common.R` | paths (`QSPATLAS_PAPERS`), output helpers |
| `R/ouerdani_2015.R` | nonmem2rx translation of `PSP4-4-660-s004.txt` (mice) and `s001.txt` (patients) |
| `R/lindauer_2017.R` | nonmem2rx translation of `PSP4-6-11-s006.txt` (mouse PK / RO / tumour) |
| `R/hutchinson_2016.R` | monolix2rx translation of `PSP4-5-636-s004.txt`, plus a hand translation of the same file |
| `R/stein_2019.R` | hand translation of the MLXTRAN block printed in `PSP4-8-285-s001.pdf` |
| `R/compare.R` | per case and variable: largest absolute and relative difference; writes `summary/` |
| `summary/summary_<id>.csv`, `summary/summary.csv` | the comparison tables (committed) |

Run from the repository root (outputs in `out/`, not committed):

```bash
julia --startup-file=no --project=. validation/authors_code/julia/atlas_curves.jl
for m in ouerdani_2015 lindauer_2017 hutchinson_2016 stein_2019; do Rscript validation/authors_code/R/$m.R; done
Rscript validation/authors_code/R/compare.R
```

Versions: R 4.6.2, rxode2 5.1.8, nonmem2rx 0.1.12, monolix2rx 0.0.6. With Homebrew R on macOS, rxode2 compiles
each model and may need `SDKROOT` set to an installed SDK (as for mrgsolve, see `validation/mrgsolve/`);
building RcppParallel from source needs `cmake` on the PATH.

Relative difference = |rxode2 - atlas| / max(|atlas|, 1e-6 x max|atlas| of the curve). rxode2 was run with
rtol 1e-10, atol 1e-12 (1e-14 for Lindauer); the atlas with its own settings.

## Results

| model | route | cases | max abs diff | max rel diff | verdict |
|---|---|---|---|---|---|
| Ouerdani 2015 | nonmem2rx (ADVAN9 renamed ADVAN13 in a temp copy) | Fig 4 mouse (AUC 1140.8) and patient (AUC 771.6); mouse AUC 220.2 and 656.8 | 9.1e-8 mm^3 (P), 5.3e-9 (K) | 5.2e-10 | agrees |
| Hutchinson 2016 | monolix2rx (header converted in a temp copy) and by hand; identical results | Fig 2a, 2c (5 doses), 2e (5 N values), Fig 6a | 1.7e-5 mm^3 | 7.6e-8 (window cases), 1.7e-10 (others) | agrees, with the N = Nmax + 1 reading below |
| Stein 2019 | by hand (code only printed in a PDF) | typical; tocilizumab day 4 + steroid day 6; steroid day 2 + tocilizumab day 6; tocilizumab day 12 | 1.2e-7 copies/ug | 2.2e-12 | agrees |
| Lindauer 2017 | nonmem2rx, no edit | mouse, 0 / 0.1 / 1.4 / 5 mg/kg on days 0, 5, 9, 13, 17, 30 days | CP 3.8e-8 mg/L, ROb 4.2e-8 %, ROt 7.1e-8 %, TV 5.1e-7 mm^3 | 1.2e-8 | agrees |

Per-curve numbers are in `summary/summary_<id>.csv` (columns `max_abs`, `at_time`, `max_rel`, `max_atlas`).

### Ouerdani 2015

- Translation: nonmem2rx refuses ADVAN9 ("Differential Algebra Equations are not supported"). Neither stream has an
  `$AES` or `$EQUILIBRIUM` block, so ADVAN9 is only a stiff ODE solver here; the script checks that and renames it to
  ADVAN13 in a temporary copy. Nothing else is changed. The data files are not published; the data items the code
  reads are given as covariates: `DOSE` (> 0 = treated), `AUC`, and the baseline `VB` (mice) / `BASE` (patients).
- Inputs: THETAs set to Table 1 (the streams carry initial values such as LAMBDA 0.144 / 0.0015), P0 150 mm^3 and
  66.8 mm, AUC 1140.8 and 771.6 ug h/mL (atlas fixes `p0`, `mouse_auc`, `patient_auc_mean`), etas zero, time in days.
- Agreement: all four treated curves agree to 5e-10 relative. The equations, the exponent n (1 in s004,
  `B*A(1)`; 0.5 in s001, `B*A(1)**0.5`), the AUC power laws and the exp(-delta T) clock from T = 0 match.
- Difference, explained: control animals. The code sets GAMMA, ALPHA and DELTA to 0 when `DOSE = 0`. The atlas has
  no dose switch; with AUC = 0 its gamma is 0 but alpha stays alpha0 because beta_alpha = 0 (AUC^0 = 1), so an
  atlas "control" still has the cytotoxic term (`mouse_control_auc0`: P differs by up to 216 mm^3, rel 1.18).
  No atlas target is a control simulation, so the grade is not affected; anyone simulating controls with the atlas
  model must switch the drug terms off as the code does.

### Hutchinson 2016

- Translation, route 1: monolix2rx reads Monolix project files, and s004 is a bare structural model in the older
  syntax (`INPUT:` block with `parameter = {...}` and `regressor = X`, no `[LONGITUDINAL]` header). The script
  replaces that block with `[LONGITUDINAL]` / `input = {...}` / `X = {use=regressor}` in a temporary copy (from `PK:`
  to the end the file is unchanged) and wraps it in a generated project with the parameters fixed. monolix2rx then
  translates the equations, but drops `t_0 = 38` (it warns: non-zero t_0 is unsupported). Route 2 is a hand
  translation of the same lines. Both routes give the same numbers to the last digit.
- Initial time: in the code the initial values T_0 = X, V_0 = K*X hold at t_0 = 38. rxode2 starts at time 0, so a
  reset event at day 38 puts the states back to their initial values there (checked in the script: T(38) = X,
  V(38) = K X). Without it the tumour grows from day 0 and is about 7 times too large at day 38.
- Dosing: `Bev` and `Van` are dose compartments with `ddt = 0`; the kill is `delta_V*dose` and the window factor
  `N = Nmax*dose + 1`. The dose records of the authors' data are not published. One dose of amount 1 at day 38 is
  used, so the kill equals delta_V and N = Nmax + 1 in the window.
- Like-for-like (code Nmax = atlas N - 1): all 12 cases agree. Fig 2a and 2c agree to 1.7e-10 relative. The window
  cases (2e, 6a) differ by up to 7.6e-8 relative at the window end (t = 62). This is the atlas's own solver
  tolerance (Tsit5, 1e-9): re-solving the atlas problem with Vern9 at 1e-13 moves it by 7.5e-8 and then agrees with
  rxode2 to 2.5e-10.
- Difference, explained: the +1. With the same number in both (code Nmax = atlas N, cases `*_samenum`), the code
  gives a window factor one larger. Effect on T: 26.7 mm^3 (19 %) for N = 2, 17.0 (10 %) for N = 3, 6.5 (3.2 %) for
  N = 6, 2.4 (1.1 %) for N = 11, 0.78 (0.3 %) for N = 21, and 2.1 mm^3 (1.4 %) for Fig 6a with Table 1 Nmax = 6.7
  (code N = 7.7). The atlas follows Eq. 5 (N = Nmax), which Fig 2e/f support (atlas report, Issues 2). Since the
  Table 1 estimates were fitted with this code, Nmax = 6.7 there means a factor of 7.7 if the data carried a single
  dose of amount 1; with repeated dose records `Bev`/`Van` would accumulate (no elimination) and both N and the kill
  would step up with each dose. This cannot be settled without the data file. The Figure 6 chemotherapy terms
  (Eqs 7-11) are not in the supplement code, so Fig 6b-d cannot be checked this way.

### Stein 2019

- Translation: the supplement has no `.mlxtran` or project file; the MLXTRAN model is printed as text in
  `PSP4-8-285-s001.pdf` ("Model MLXTRAN file"). Translated by hand into an algebraic rxode2 model; the script checks
  with `pdftotext` that the translated lines are still in the PDF. Two reading notes: `R1 = R0*exp(rho* t1)` is
  broken over three lines in the PDF text layer; and the third branch is printed `elseif t < tmax` (lower case)
  while the parameter is `Tmax`. Mlxtran is case sensitive, so as printed `tmax` is undefined; it is read as `Tmax`.
  Without comedication t1 = t2 = Tmax and that branch is never reached.
- Agreement: the atlas ODE (effector / memory compartments) and the code's closed form agree to 2.2e-12 relative for
  the typical curve and for three comedication orders, including the clipping of a comedication after Tmax
  (`toci12` equals `typical`) and the `TOCI1T <= STERSTT` ordering rule. No difference.

### Lindauer 2017

- Translation: nonmem2rx, no edit. The stream is the first-estimate run: THETA(20)-(22), (26)-(31) are initial
  values and Kon/Koff carry 340.2 and 0.10584. They are set to the atlas mouse values on the code's own scales
  (lnTmulti = log 4.32, lnEMAXTP = log 94.7, lnEC50TP = log(1.46e-3 uM), Kon 340, Koff 0.106, KdegPD1 0.0194,
  L0 0.113, L1 187, lnSLtg = log 1.98e-5, lngamma = log 2.28). Doses enter as in the code: AMT in ug into C1,
  converted to nmol by `F1 = 1000/MW`; 20 g mouse.
- Agreement: CP, ROb, ROt and TV agree to 1.2e-8 relative at all four doses, against the atlas run with the
  `IF(ROtumor>0.01)` branch on (`guard = true`) and also with it off (the atlas default): the branch changes TV by
  less than 1e-9 relative here. Values exactly at dose times are left out (pre- or post-dose depends on the tool).
- Notes, no numerical effect: the `$MODEL` compartment names in s006 do not match their use in `$DES` (A(3), named
  Ce_ub, is the vascular amount; A(5), named Cis, is Ce_b; and so on). Both translations follow the indices, so the
  labels do not matter. The human scenario of the atlas (Figure 4) is not in the code and is not checked here.

## Verdicts

- Ouerdani 2015: the atlas equations and inputs match the authors' NONMEM code for every graded curve. The code's
  control-arm switch (`DOSE = 0` turns all drug terms off) is not in the atlas; it matters only for control runs.
- Hutchinson 2016: the atlas matches the authors' MLXTRAN code once the code's window factor `Nmax*dose + 1` is
  mapped to the atlas N; with the same number the code's factor is one larger (up to 19 % on T at N = 2). The
  data coding that would decide which Nmax Table 1 reports is not published.
- Stein 2019: exact agreement with the printed MLXTRAN, including the comedication branches; one case slip
  (`tmax`) in the printed code, without effect.
- Lindauer 2017: the atlas mouse model matches the authors' NONMEM stream to solver precision.
