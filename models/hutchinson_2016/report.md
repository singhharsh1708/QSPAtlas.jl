# hutchinson_2016: grade C (Figure 2: B, Figure 6: C)

Hutchinson LG, Mueller HJ, Gaffney EA, Maini PK, Wagg J, Phipps A, Boetsch C, Byrne HM, Ribba B. Modeling longitudinal preclinical tumor size data to identify transient dynamics in tumor response to antiangiogenic drugs. CPT Pharmacometrics Syst Pharmacol 2016;5:636-645. DOI 10.1002/psp4.12142, PMC5192995.

Grade C overall. Figure 2 (22 curves) reproduces after one documented typo fix (Figure 2 alone would be B, with no input read from a figure). Figure 6 (8 curves) does not reproduce even with two further fixes (T0 read from the figure, resistance clock from the start of chemotherapy): the printed alpha~_C = 1 wipes out the tumour in panels b and d, and panel a (no drug) misses its peak. A value of alpha~_C = 0.01 brings panels b to d close, but it is a model parameter that can only be inferred from Figure 6 itself, so it is reported as a sensitivity run and not as a fix (see Fixes, Issues 3-6, notes.md). Both figures were frozen as targets together before any comparison, so both are graded; the per-figure split is in `score.toml` (`figure_grades`). Note: `score.toml` shows Figure 2 as B† because the library attaches a run's fix list, which includes the Figure 6 T0 fix, to every figure; no Figure 2 curve uses it (Figure 2 numbers are identical to the earlier run).

## Model

Two ODEs, tumour volume T and vessel-dependent carrying capacity V (both mm^3), time in days after tumour cell injection, start and first dose at day 38.

- Control (Eqs 1-2): dT/dt = alpha_T T (1 - T/V); dV/dt = alpha_V T^(2/3) V^(1/3).
- Treatment (Eqs 3-5): dT/dt = alpha_T T (1 - T/(N V)); dV/dt = alpha_V T^(2/3) V^(1/3) - delta_V V; N(t) = Nmax for tnorm1 < t < tnorm2, otherwise 1. Antiangiogenic drug is not modelled, it is assumed fully effective from day 38.
- Initial conditions: T(38) = T0, V(38) = K T0 (the paper prints K = T0/V0, see Fixes).
- Cytotoxic drug (Eqs 7-11): dT/dt gets an extra - delta_T C T exp(-lambda t); C = alphaC(t) N V exp(-k mod(t,1)); alphaC(t) = alpha~_C between tCon and tCoff, otherwise 0.

`model.jl` has one system with all of these terms (parameters switch them off). Fig 2 parameters come from the caption (alpha_T 0.1, alpha_V 0.09, K 2, T0 70; delta_V 0.04-0.12 with Nmax = 1; tnorm 52/62, delta_V 0.1, Nmax 2/3/6/11/21, read from the legends). The ODE right-hand side is discontinuous (window edges, drug sawtooth), so the solver is forced to stop on those times (`breakpoint_times`, implemented with a zero-amount `bolus_callback`).

## Targets

- Figure 2 panels a-f, PSP4-5-636-g002.jpg (667 px wide): 22 simulation lines (a: T, b: V control; c, d: five delta_V values; e, f: five Nmax values). Traced by colour (MATLAB default colours) with `reference/digitise_hutchinson.jl` (the repository digitiser plus rectangle exclusions for the legends), settings in `reference/digitise.toml`. Curves hidden under the top (green) curve are only digitised where visible. Calibration: two ticks per axis, checked on 2-3 other ticks per panel (max error 3 mm^3 on the 100-1000 V axis, otherwise below 0.5 % of an axis); overlay in `figures/digitise_check.png` (local only: it contains the paper's figure).
- Tolerance: atlas default (RMS 0.015, P95 0.03 of the axis box), no per-curve changes.
- Frozen: 2026-11-04 about 09:48 IST, before any simulation was compared.
- Figure 6: T (black) and V (magenta) in the four panels, 8 curves, digitised the same way (`reference/digitise_fig6.toml`). The dashed blue drug curve C is not digitised (daily sawtooth, clipped at the axis top). Graded. The paper gives no T0 for this figure; the as-published run uses T0 = 70 mm^3 (the Methods mean tumour size at day 38) and the with-fixes run uses T0 = 80 read from the figure (Fixes item 3). One cleanup happened after the first comparison run: stray points on the blue drug sawtooth in 6b T were removed with an exclusion rectangle. It does not change any conclusion.

## Result

Source: `score.toml`, written by `reproduce()` (graded: Figures 2 and 6, 30 curves).

| curve | as published rms | p95 | pass | with fix rms | p95 | pass |
|---|---|---|---|---|---|---|
| fig2_a_T | 0.2787 | 0.5256 | no | 0.001811 | 0.003272 | yes |
| fig2_b_V | 0.2639 | 0.4613 | no | 0.001709 | 0.002951 | yes |
| fig2_c_T_d04 | 0.4219 | 0.6108 | no | 0.001472 | 0.002366 | yes |
| fig2_c_T_d06 | 0.3615 | 0.4702 | no | 0.001847 | 0.003072 | yes |
| fig2_c_T_d08 | 0.303 | 0.3539 | no | 0.002185 | 0.00514 | yes |
| fig2_c_T_d10 | 0.2634 | 0.2753 | no | 0.005971 | 0.0042 | yes |
| fig2_c_T_d12 | 0.1988 | 0.2375 | no | 0.002432 | 0.004262 | yes |
| fig2_d_V_d04 | 0.4258 | 0.5967 | no | 0.001799 | 0.003519 | yes |
| fig2_d_V_d06 | 0.3343 | 0.3906 | no | 0.002308 | 0.003818 | yes |
| fig2_d_V_d08 | 0.2627 | 0.3151 | no | 0.002185 | 0.003324 | yes |
| fig2_d_V_d10 | 0.2066 | 0.2795 | no | 0.002287 | 0.003807 | yes |
| fig2_d_V_d12 | 0.1814 | 0.3057 | no | 0.002036 | 0.003328 | yes |
| fig2_e_T_N11 | 0.3538 | 0.4785 | no | 0.002554 | 0.004483 | yes |
| fig2_e_T_N2 | 0.2695 | 0.3017 | no | 0.0072 | 0.01188 | yes |
| fig2_e_T_N21 | 0.2738 | 0.4769 | no | 0.006069 | 0.01161 | yes |
| fig2_e_T_N3 | 0.3006 | 0.3633 | no | 0.004661 | 0.009077 | yes |
| fig2_e_T_N6 | 0.3143 | 0.4139 | no | 0.006418 | 0.009809 | yes |
| fig2_f_V_N11 | 0.7034 | 0.7527 | no | 0.005366 | 0.008739 | yes |
| fig2_f_V_N2 | 0.5835 | 0.5966 | no | 0.006759 | 0.01232 | yes |
| fig2_f_V_N21 | 0.6991 | 0.9282 | no | 0.008787 | 0.01917 | yes |
| fig2_f_V_N3 | 0.6136 | 0.6582 | no | 0.004831 | 0.008304 | yes |
| fig2_f_V_N6 | 0.6407 | 0.7062 | no | 0.004077 | 0.007156 | yes |

As published all 22 fail because V0 = T0/K = 35 mm^3 instead of 140 mm^3. With the fix all 22 pass, worst rms 0.0088 (fig2_f_V_N21), worst p95 0.019. The lines for Nmax = 2 and 21 in panels e and f deviate slightly in a systematic way at the peak and tail (about 2-4 mm^3), well inside tolerance; I did not chase it.

Figure 6 (8 curves), 2026-11-04 update. As published = printed inputs, T0 = 70, K as printed. With fixes = K fix, resistance clock from tCon, T0 = 80, alpha~_C = 1 as printed.

| curve | as published rms | p95 | pass | with fixes rms | p95 | pass |
|---|---|---|---|---|---|---|
| fig6_a_T | 0.1678 | 0.2621 | no | 0.0487 | 0.0878 | no |
| fig6_a_V | 0.1814 | 0.2447 | no | 0.0374 | 0.0746 | no |
| fig6_b_T | 0.1898 | 0.2631 | no | 0.3864 | 0.5975 | no |
| fig6_b_V | 0.1939 | 0.2302 | no | 0.3404 | 0.5460 | no |
| fig6_c_T | 0.0827 | 0.1312 | no | 0.0986 | 0.1855 | no |
| fig6_c_V | 0.1402 | 0.1917 | no | 0.0798 | 0.1627 | no |
| fig6_d_T | 0.2314 | 0.2750 | no | 0.3195 | 0.3852 | no |
| fig6_d_V | 0.2215 | 0.2595 | no | 0.2951 | 0.3804 | no |

Figure 6 does not match with the allowed fixes. Panel a (no drug, Table 1 parameters only) improves with T0 = 80 but still fails (peak about 147.5 mm^3 against about 165 in the figure). Panels b and d get worse than as published: with alpha~_C = 1 the tumour is eradicated, and the resistance clock from tCon removes the factor that had been switching the drug off.

Sensitivity runs for Figure 6 (not fixes, not in `score.toml`; `sensitivity()` in `reproduce.jl`), all with the same tolerances:

| run | a T | a V | b T | b V | c T | c V | d T | d V |
|---|---|---|---|---|---|---|---|---|
| alpha~_C = 0.01, T0 = 80, clock from tCon (rms) | 0.049 | 0.037 | 0.011 | 0.014 | 0.019 | 0.022 | 0.019 | 0.015 |
| pass | no | no | yes | yes | no | no | no | no |
| alpha~_C = 0.01, T0 = 80, absolute clock (rms) | 0.049 | 0.037 | 0.100 | 0.073 | 0.312 | 0.189 | 0.226 | 0.145 |
| alpha~_C = 0.01, T0 = 70, clock from tCon (rms) | 0.109 | 0.106 | 0.042 | 0.062 | 0.055 | 0.073 | 0.022 | 0.054 |

With alpha~_C = 0.01 and both other choices, panels b to d are within a few mm^3 of the figure (an independent rebuild found the same), but only 6b T and 6b V pass the atlas tolerance (p95 0.019 and 0.026); 6c T/V (rms 0.019/0.022) and 6d T/V (rms 0.019/0.015, p95 0.030/0.027) just fail, and panel a fails regardless. So even if alpha~_C = 0.01 were admitted as a fix, Figure 6 would stay C. The three choices only work together: the absolute clock with alpha~_C = 0.01 fails all of b to d, as does T0 = 70. Diagnostic only: reading N = Nmax + 1 in the chemotherapy panels on top of alpha~_C = 0.01 would pass 6b T, 6c T, 6d T and 6d V, but this is not a fix (see Issues 2) and Figure 2 rules out N = Nmax + 1.

## Fixes

1. `K_is_V_over_T` (typo): V0 = K T0 instead of the printed K = T0/V0. Evidence, independent of Fig 2: supplement MLXTRAN code `PSP4-5-636-s004.txt` has `T_0=X`, `V_0 = K*X`; Table 1 describes K as "Initial ratio of V/T". Consistency check: Fig 2b starts at V = 140 = 2 x 70.

Figure 6 only (the "with fixes" run applies these to the Figure 6 scenarios; Figure 2 is simulated exactly as before):

2. `resistance_clock_from_tCon` (:ambiguity): the resistance factor of Eq. 7 is exp(-lambda (t - tCon)), timed from the start of chemotherapy, not exp(-lambda t) with t in days since inoculation. Evidence: Eq. 7 writes exp(-lambda t) and never defines where t starts; the text says resistance is likely to influence Fig 6(d) compared with Fig 6(c), which fits a clock that runs from the start of the course (6d has a three-week course, 6c one week). Both readings are textually possible; no number is changed. Caveat: the choice was checked against the figure, and with the printed alpha~_C it does not improve the match (see Result), so it earns no grade by itself.
3. `T0_fig6_from_figure` (:from_target): T0 = 80 mm^3 for all four Figure 6 panels. The paper gives no T0 for this figure (population parameters exclude it, Methods give a mean of 70). The value is read from the starting point of the T curves in Figure 6 itself, so it is a protocol input read from the target and carries the dagger. The independent rebuild read the same 80. A model parameter was not read from the figure.

Not a fix, sensitivity run only: alpha~_C = 0.01 in place of the printed 1 (see notes.md for the reasoning). No other fix. In particular N = Nmax (Eq. 5 as printed) is what Fig 2 shows (see Issues item 2), so the supplement code reading is not applied.

## Issues found in the paper

1. Methods, before Eq. 6: K = T0/V0 is inverted (should be V0/T0), see Fix 1.
2. Eq. 5 gives N = Nmax in the window, but the supplement code (`s004.txt`) has `N = Nmax*(Bev+Van)+1`. With a unit dose amount that is N = Nmax + 1, so Table 1 Nmax = 6.7 might mean a factor 7.7. Diagnostic (not a fix): Fig 2e/f with N = Nmax + 1 fails for Nmax = 2, 3 (rms 0.058, 0.025) and 6 (0.014-0.016, borderline), so Fig 2 was made with Eq. 5 as printed. Also in the code Bev and Van have ddt = 0 (accumulating dose compartments), so with repeated weekly doses both N = Nmax*dose + 1 and the delta_V*dose term would step up with each dose, unless the data file holds a single dose record of amount 1; and N = Nmax*dose + 1 gives N = 1 for controls and applies the window to treated animals only, whereas Eq. 5 is written for everyone. Any implementation must take a side. The dose-amount coding is not in the supplement, so the Table 1 delta_V and Nmax cannot be checked against the data fit.
3. Figure 6: T0 is not given. The population parameters do not include it. Methods give a mean size of 70 mm^3, but the curves start at T about 80 and V about 91-96 (K = 1.14 would imply T0 about 83). With T0 = 70 the start is 10 mm^3 off. The with-fixes run uses 80 (Fix 3, :from_target).
4. Figure 6: the drug scale does not match Eq. 9 and alpha~_C = 1 mg/mL/mm^3. Eq. 9 gives C of about 78-93 for the daily peak in 6b (C = V at the dose, N = 1); the figure's C axis runs 0-5 and the pre-window peak is about 0.9, a factor of about 100. Scaling alpha~_C by 0.01 matches that C scale but then the drug has almost no effect on T (6b: T(49) = 89 vs 65 in the figure); the printed value kills far too much (T(49) = 33). Neither reading reproduces both C and T with T0 = 70 and the absolute resistance clock (diagnostic only, not fitted). Update 2026-11-04: with T0 = 80 and the resistance clock from tCon, alpha~_C = 0.01 does give both the C scale of the figure and T and V in panels b to d within a few mm^3 (sensitivity table in Result). It stays a sensitivity run: alpha~_C is a model parameter, the printed value is explicit (text, 1 mg/mL/mm^3), and the only evidence for 0.01 is the C axis and the tumour curves of Figure 6, the figure being reproduced. The journal units of alpha~_C and delta_T are also garbled (the preprint text s001 has the consistent forms), so a unit slip is possible, but no source gives 0.01.
5. Eq. 7: exp(-lambda t) is written with time t of undefined origin. Earlier note: with the printed alpha~_C = 1 and T0 = 70, starting the clock at day 38 kills the tumour completely. With alpha~_C = 0.01 (sensitivity run) the opposite holds: the absolute clock switches the drug off (exp(-0.08 x 54) = 0.013) and fails panels b to d, while the clock from tCon matches. The with-fixes run uses the clock from tCon (Fix 2).
6. Figure 6a (no drug) with T0 inferred as 83 still gives T(60) = 148 vs about 165 in the figure and T(71) = 128 vs 137, so even the drug-free panel is not exactly reproduced from Table 1 (the figure may use other parameter values than the population estimates). With T0 = 80 (Fix 3) the peak is 147.5 at day 59.2 against about 165 in the figure; the rebuild found the same. Panel b has the same window and matches, and the figure reaches the 165 peak if the window ends near 60.5 or if Nmax is very large, so panel a was probably made with other window or Nmax settings. Not fixable from the paper.
7. BIC of the monotonic model is 2,759 in the text before Table 1 and 2797 later (BIC_M).
8. The LRT statistic is printed as -Delta 2LL = -32764.61 (sign and form unclear).
9. Fig 2 caption gives delta_V as a range 0.04-0.12 and Nmax as 2-21; the five values are only in the plot legends.

## Self-check (2026-11-04)

An independent rebuild from the paper alone (plain ODE solver, no code shared with this model) agreed with Figure 2 to within 2e-5 of each curve's range, and found the Figure 6 reading above: alpha~_C = 0.01 (the figure's own C axis implies it; the printed 1 wipes out the tumour), the resistance factor timed from the start of chemotherapy, and T0 = 80 read from the figure. With these panels b to d come within about 5 mm^3, panel a still peaks near 147.5 instead of about 165. This atlas adopts the clock and T0 (a documented ambiguity and a flagged protocol input) but not alpha~_C, which is a model parameter inferred from the target figure; hence Figure 6 stays C.

## Review

Independent review, 2026-11-04. The review found the equations, the K fix and the Figure 2 digitising sound. It rejected the original grade B, which graded Figure 2 only: Figures 2 and 6 were frozen as targets together at about 09:48 and Figure 6 was relabelled "secondary" only after the first comparison showed it failing. Excluding a target after seeing it fail is not allowed, so Figure 6 is graded and the model is C, with the per-figure split (Figure 2: B, Figure 6: C) reported. Other review changes: the Figure 2b V axis was declared 100-1000 but is drawn 100-900 (frame rows map to 103 and 897), corrected; an unsupported claim that panel 6a "fails with any T0" was removed from the notes (only T0 = 70 and 83 were tried).

## What would raise the grade

Figure 2 is at B: the printed K = T0/V0 fails every Figure 2 curve, and the fix is well evidenced. Figure 6 (and with it the model grade) would need a source for alpha~_C other than the figure itself (the authors' value, an erratum, or the code used for Fig 6), the authors' T0, the origin of the resistance clock, and the parameter set used for panel 6a. The code and values used for Fig 6 are not in the supplement.

## Literature check (2026-11-04)

The first author's DPhil thesis (L. G. Hutchinson, University of Oxford, Mathematical Institute; ORA uuid:0c1034ee-ff8e-4ed4-b168-f559f93f9142) contains Figure 6 as vector graphics. In it the plotted drug concentration is exactly 0.01 N V, that is Eq. 9 with alpha~_C = 0.01 rather than the printed 1, and T0 is 80 mm^3 with V0 = 96 (K = 1.2). This supports the sensitivity run above, but the value still comes from the figure being reproduced, so it is not used as a fix; Figure 6 remains C (panel a does not reproduce with either N = Nmax or Nmax + 1). Both texts state N = Nmax.
