# Recalibration of hutchinson_2016 against its own data (Figure 1d)

Recalibration does not change the atlas grade (hutchinson_2016 stays C). Written 2026-11-04, phase 1.

## Question

Do the published population parameters (Hutchinson et al. 2016, Table 1, final model) agree with the paper's own tumour data, and can they be identified from it?

## Data

Figure 1d of the paper (`PSP4-5-636-g001.jpg`) shows the mean tumour volume of three groups of 10 mice (Control = omalizumab, Bevacizumab, Vanucizumab) on days 38, 42, 45, 49, 52, 56, 59, 63, 67, 71. The marker centres were digitised (not the connecting lines) with the atlas digitiser helpers: colour classification, bounding-box midpoint of each marker, calibration on the axis ticks (1 px = about 5.7 mm^3, estimated error about 1 px). Day 38 Control and Bevacizumab markers lie under the Vanucizumab triangle and were set to 70 mm^3 (Methods: mean size 70 at day 38). Control days 42 and 45 are partly covered by the green markers and were read by eye (85 and 95, error about 6 mm^3). Overlay for checking: `hutchinson_fig1d_overlay.png`. Data: `../data/hutchinson_fig1d.csv`, script `../digitise_fig1d.jl`.

## Model and what was estimated

Eqs 1-5 with the K typo fix (V0 = K T0). Treatment reading: the paper assumes full drug effect from the first dose at day 38 (no PK), so delta_V is constant from day 38 and the weekly dosing arrows (38, 45, 52, 59, 66) do not enter the equations. The supplement MLXTRAN code has Bev/Van dose amounts with zero dynamics and the dose coding is not available, so the code reading could not be checked. Control: delta_V = 0, N = 1. Treated: N = Nmax in (53.4, 59.2).

Estimated (7): alpha_T, alpha_V, K (shared by all groups), Nmax (shared by the two treated groups), delta_V for bevacizumab and vanucizumab, and sigma (log-normal residual: log y ~ Normal(log T, sigma)).
Fixed: T0 = 70 mm^3 for all groups (the day-38 means are 70, 70, 71; estimating them would only fit digitising noise; the day-38 points are not used in the likelihood, 27 observations remain), window edges tnorm1 = 53.4 and tnorm2 = 59.2 (Table 1; only two sampling days, 56 and 59, fall inside the window), beta = 2/3, gamma = 1/3. The plan was written before fitting (`../notes.md`).

## PEtab.jl (maximum likelihood)

50 Latin-hypercube starts with Fides (BFGS), log10 scale, no prior. All 50 starts reach the same optimum (nllh 103.39), so there is no sign of local minima. The window is implemented as two PEtab time events, with time counted from day 38 (an SBMLImporter limitation, see `../notes.md`). PEtab and a plain `remake` simulation agree to 1e-8.

Profile likelihoods (`fig_profiles.png`): each parameter fixed on a grid and the rest re-fitted with Fides.

## Turing.jl (NUTS)

Priors: log-normal with median at the Table 1 value and sd 1 on the log scale (Nmax truncated at 1; sigma median 0.2). 4 chains x 1000 samples after 500 warm-up, target acceptance 0.8, ForwardDiff, Tsit5 (3 segments of the MTK problem via `remake`). Runtime 250 s.
Diagnostics: R-hat 1.001-1.007, bulk ESS 700-1900 (all parameters), 1 divergent transition out of 4000, mean tree depth 5.6. Note the priors are centred on the published values, which favours them; a second run with generic priors not taken from the paper (medians alpha_T 0.1, alpha_V 0.1, K 1, Nmax 3, delta_V 0.05, sd 1.5 on the log scale; `turing_*_generic.csv`) also finished (R-hat <= 1.006, ESS >= 595, no divergences). Its medians [90% CI]: alpha_T 0.101 [0.067, 0.208], alpha_V 0.157 [0.084, 0.470], K 1.06 [0.85, 1.40], Nmax 8.7 [1.9, 72], delta_V bev 0.116 [0.053, 0.331], delta_V van 0.169 [0.080, 0.472], sigma 0.109. Same picture with wider intervals; all published values are inside.

## Comparison (Table 1 vs this work)

| parameter | published | PEtab MLE | PEtab profile 95% (pointwise, 2 dNLL < 3.84) | Turing median [90% CI] |
|---|---|---|---|---|
| alpha_T (1/day) | 0.109 | 0.090 | about 0.055 to 0.37 | 0.097 [0.070, 0.150] |
| alpha_V (1/day) | 0.119 | 0.196 | above 0.067, no upper limit | 0.167 [0.103, 0.376] |
| K | 1.14 | 1.12 | about 0.1 to 2.0 | 1.08 [0.87, 1.41] |
| Nmax | 6.7 | at bound (1e4) | above about 1.5, no upper limit | 10.9 [2.8, 44] |
| delta_V bev (1/day) | 0.113 | 0.149 | above 0.039, no upper limit | 0.123 [0.071, 0.271] |
| delta_V van (1/day) | 0.115 | 0.214 | above 0.063, no upper limit | 0.179 [0.106, 0.386] |
| sigma | 0.185 (individual-level b) | 0.094 | | 0.109 [0.087, 0.142] |

Profile limits are read from the grid by linear interpolation (`petab_profiles.csv`); K lower limit is approximate. Numbers also in `estimates.toml`.

Plots: `fig_fits.png` (data, published simulation, PEtab MLE, Turing 90% band per group), `fig_profiles.png`, `fig_posterior.png` (posterior vs prior), `fig_pairs.png` (posterior correlations).

## Conclusion

- Each published value lies inside both the Turing 90% interval and the PEtab pointwise 95% profile interval (alpha_V and delta_V van are near the lower edge of the Turing interval). In that sense no single published value is contradicted by the means.
- Taken together the published vector fits the group means worse than the MLE: 2 (nllh_published - nllh_MLE) = 33.9 (published vector with sigma re-optimised, nllh 120.3 vs 103.4). This is exploratory, not a formal test: the data are digitised group means fitted with a typical-individual model, repeated measures on the same animals are treated as independent, some parameters are unidentified or at a bound, and sigma is re-estimated under the published vector, so a chi-square reference distribution with 7 degrees of freedom does not strictly apply. A simulation-calibrated test, population-mean predictions and digitising uncertainty are needed before drawing a conclusion. The simulation with Table 1 values undershoots the Bevacizumab mean after day 59 (about 110 vs 140-180) and the early Vanucizumab means (see `fig_fits.png`). The published delta_V values are almost equal (0.113, 0.115) while the Vanucizumab mean is clearly lower than Bevacizumab. This does not show that Table 1 is wrong (see limitations).
- Not identifiable from three mean curves: Nmax (only a lower bound near 1.5; Nmax = 1, i.e. no window, is rejected, 2 dNLL = 16; the posterior is mostly the prior); alpha_V, delta_V bev and delta_V van form one ridge (posterior correlations 0.98-0.99 on the log scale, alpha_V vs alpha_T -0.90): the data constrain the net vessel balance, not the individual rates, and upper limits are open. alpha_T and K are constrained within a factor of about 2-3. sigma is well determined.

## Limitations

- Group means are not typical individuals. Table 1 values are population typical values with large between-mouse variability (IIV 0.42 on alpha_T, 0.71 on K); the mean of a nonlinear population is not the typical individual, so a difference between the fit and Table 1 is expected even if both are right.
- Digitisation error about 6 mm^3 per point, plus three hand-handled or hidden markers (above). The 10 mice per group are summarised by one mean per day, no error bars, so sigma is a model-plus-digitising error, not mouse variability.
- Window edges and T0 fixed; Nmax is confounded with the window edges (not tested). A T0-estimated sensitivity run (`T0MODE=shared|group` in `petab.jl`) was **not run in phase 1**.
- Treatment reading (full effect from day 38, dosing not modelled) is the paper's text, not verified against the original data file.
- Turing priors centred on the published values bias toward them; the generic-prior run (above) gives the same conclusions with wider intervals.
- Profiles hold for the fixed window; Nmax profile is cut at the bound 1e4.

## Files

`common.jl`, `petab.jl`, `turing.jl`, `plots.jl`, `digitise_fig1d.jl`; results in this folder; planning notes in `../notes.md`.
