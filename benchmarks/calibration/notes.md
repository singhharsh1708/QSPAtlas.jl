# Recalibration of hutchinson_2016 against Figure 1d (PEtab.jl MLE and Turing.jl NUTS)

Notes written 2026-11-04 (IST). The plan below was written BEFORE any fit was run (digitised data frozen first,
`data/hutchinson_fig1d.csv`, 2026-11-04 about 10:00 IST).
Recalibration never changes the atlas grade (hutchinson_2016 stays C). Results are in `results/calibration.md`.

## Plan (written 2026-11-04 ~10:15 IST, before fitting)

Data: 3 group-mean curves (Control, Bevacizumab, Vanucizumab), 10 days each (38, 42, 45, 49, 52, 56, 59, 63, 67, 71),
digitised markers of Fig 1d. n = 10 mice per group, the curves are arithmetic means of individual volumes.
The day-38 values are the common starting size (about 70 mm^3, Methods) and are used only to set T0, not as observations,
so 27 observations enter the likelihood (9 per group).

Model reading (stated, not tested):
- Treatment: the paper assumes antiangiogenic drug is fully effective from the first dose at day 38 (Methods/Model
  development: "from the onset of treatment, the concentration ... is maintained at a high enough level to ensure maximum
  efficacy", no PK). So delta_V is constant from t = 38; the weekly arrows (38, 45, 52, 59, 66) do not enter the equations.
  Time starts at 38 with the treatment already on. The supplement MLXTRAN code has Bev/Van amounts with iv doses and
  ddt = 0 (so the amount would be a dose-coded constant, possibly accumulating with each dose; the data file is not
  available). The paper reading is used (indicator = 1 from day 38), i.e. N = Nmax in the window and vessel kill delta_V
  throughout. Table 1 delta_V and Nmax are defined for that reading.
- Control: Eqs 1-2 (delta_V = 0, N = 1). Treated: Eqs 3-5 with N = Nmax (Eq. 5 as printed) in (tnorm1, tnorm2).
- V0 = K T0 (typo fix, as in the hutchinson_2016 atlas model).
- Table 1 (final model population values): alpha_T 0.109, alpha_V 0.119, K 1.14, tnorm1 53.4, tnorm2 59.2, Nmax 6.7,
  delta_V(bev) 0.113, delta_V(van) 0.115, proportional error b 18.5 % (individual-level residual).
  alpha_T, alpha_V, K are shared by all three groups; Nmax, tnorm1, tnorm2 are shared by the two treated groups ("the same
  values of Nmax, tnorm1 and tnorm2 for both treatment groups", results section); delta_V differs per drug.

Estimated parameters (7): alpha_T, alpha_V, K, Nmax, delta_V_bev, delta_V_van, sigma (log-normal residual sd).
Fixed:
- T0 = 70 mm^3 for all groups (first observed means are 70, 70, 71 mm^3 within digitising error and the paper says
  mean 70 at day 38; estimating three T0 from three single points would add 3 parameters that only absorb the
  digitising noise at day 38). Sensitivity run with T0 estimated (one shared T0 and three group T0s) in PEtab only.
- tnorm1 = 53.4, tnorm2 = 59.2 (Table 1). The data have only 2 sampling days inside the window (56, 59) and
  the window edges are not identifiable from 9 mean points per group; the paper itself fixed their IIV at 0.1 and
  they only get a small RSE from 30 individual curves. Nmax and the window edges are also confounded (a long short-lived
  boost vs a short strong one), so only Nmax is estimated, given the edges.
- beta = 2/3, gamma = 1/3 (paper fixed these on geometric grounds).
No cytotoxic terms (Fig 6 not used).

Error model: log-normal, log y = log T(t) + sigma eps. (The paper used a proportional error on the individual level, b =
18.5 %; for means of 10 mice the residual is expected to be smaller.)

PEtab.jl: MTK system with parameters alpha_T, alpha_V, K, T0, delta_V, Nmax_g, Nw (window factor, event-driven);
events t == tnorm1 -> Nw = Nmax_g, t == tnorm2 -> Nw = 1; PEtabConditions for control / bev / van (t0 = 38),
log10-scale estimation, bounds alpha_T, alpha_V in [1e-3, 5], K in [0.05, 50], Nmax in [1, 200], delta_V in [1e-3, 3],
sigma in [1e-3, 3]. Multistart 50 starts (Latin hypercube) with Fides (BFGS/ Hessian as supported), best fit and spread of
the top starts, profile likelihood (own implementation: fix one parameter on a grid, re-optimise others) if time.
MLE uses no prior.

Turing.jl: @model with the same ODE solved inside (remake of the MTK problem, one solve per group, the window as
two segments), AD ForwardDiff, solver Tsit5 or Rodas5P with tight tolerances. Priors: LogNormal centred (median) at the
Table 1 values with sd 1 on the log scale (Nmax truncated at 1; sigma LogNormal(log 0.2, 1)). NUTS target acceptance 0.8,
4 chains x 1000 (after 500 warm-up) if time, else fewer. Priors centred on the published values are an informative choice
that favours the published values; the sd of 1 is wide (factor e = 2.7 one sd), and prior sensitivity is checked with a
second set of priors that do not use the published values (generic priors) if time allows.

Report: table Table 1 vs MLE vs posterior median [90 % CI], R-hat, ESS, pairwise correlations, plots with the published
parameter simulation, PEtab fit and Turing posterior predictive band.

Pre-registered success criterion for "published agrees with the data": the Table 1 value of a parameter lies inside the
90 % posterior interval and the Table 1 parameter vector has a likelihood within the 95 % joint LR region of the MLE
(2 [nllh(published) - nllh(MLE)] < chi2_7,0.95 = 14.07). Also reported: Table 1 simulation vs the data directly.

## Log
- 10:00 Fig 1d markers digitised (`digitise_fig1d.jl`); the overlay on the paper's figure (kept locally only) was checked by eye.
  Hidden markers: day 38 Control and Bevacizumab (under the Vanucizumab triangle, set to 70), Control days 42 and 45
  (partly under green markers, read by eye 85 and 95).
- 10:17-10:35 PEtab: events trigger needs t0 = 0 (SBMLImporter adds a tstop at min(event times)/2, error "tstop behind current time" with t0 = 38), so PEtab uses time since day 38. MLE found by 50/50 starts at nllh 103.39; Nmax runs to its bound (raised from 200 to 1e4), flat likelihood. Profiles first done with Optim Fminbox (poorly converged, jagged), redone with remake + Fides (smooth).
- 10:25-10:30 Turing 4x1000: Rhat <= 1.007, ESS >= 700, 1 divergence.
- 10:36 generic-prior Turing run started (results in `results/turing_*_generic.csv`); the T0-estimated sensitivity was not run.
- Env additions: Fides (brings .CondaPkg in benchmarks/calibration), ForwardDiff, MCMCChains, StatsBase, CSV, ComponentArrays, DynamicPPL.
