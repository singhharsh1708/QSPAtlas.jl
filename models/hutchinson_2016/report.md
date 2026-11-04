# hutchinson_2016: grade C (Figure 2: B, Figure 6: C)

Hutchinson LG, Mueller HJ, Gaffney EA, Maini PK, Wagg J, Phipps A, Boetsch C, Byrne HM, Ribba B. Modeling longitudinal preclinical tumor size data to identify transient dynamics in tumor response to antiangiogenic drugs. CPT Pharmacometrics Syst Pharmacol 2016;5:636-645. DOI 10.1002/psp4.12142, PMC5192995.

Grade C overall. Figure 2 (22 curves) reproduces after one documented typo fix (Figure 2 alone would be B). Figure 6 (8 curves) does not reproduce: its initial tumour size is not stated, and the drug concentration scale is about 100 times off from Eq. 9 with the printed alpha~_C, so no reading found reproduces both drug and tumour (see Issues 3-6). Both figures were frozen as targets together before any comparison, so both are graded; the per-figure split is in `score.toml` (`figure_grades`).

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
- Figure 6: T (black) and V (magenta) in the four panels, 8 curves, digitised the same way (`reference/digitise_fig6.toml`). The dashed blue drug curve C is not digitised (daily sawtooth, clipped at the axis top). Graded. T0 = 70 mm^3 is the Methods mean tumour size at day 38; the paper gives no T0 for this figure. One cleanup happened after the first comparison run: stray points on the blue drug sawtooth in 6b T were removed with an exclusion rectangle. It does not change any conclusion.

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

Figure 6, as a check only (inputs as printed, T0 = 70, K fix applied): rms 0.05-0.24 for all eight curves (6c T 0.053 best, 6d T 0.241 worst). Even panel 6a, which has no drug and needs only Table 1 parameters, fails (rms 0.109 with T0 = 70).

## Fixes

1. `K_is_V_over_T` (typo): V0 = K T0 instead of the printed K = T0/V0. Evidence, independent of Fig 2: supplement MLXTRAN code `PSP4-5-636-s004.txt` has `T_0=X`, `V_0 = K*X`; Table 1 describes K as "Initial ratio of V/T". Consistency check: Fig 2b starts at V = 140 = 2 x 70.

No other fix. In particular N = Nmax (Eq. 5 as printed) is what Fig 2 shows (see Issues item 2), so the supplement code reading is not applied.

## Issues found in the paper

1. Methods, before Eq. 6: K = T0/V0 is inverted (should be V0/T0), see Fix 1.
2. Eq. 5 gives N = Nmax in the window, but the supplement code (`s004.txt`) has `N = Nmax*(Bev+Van)+1`. With a unit dose amount that is N = Nmax + 1, so Table 1 Nmax = 6.7 might mean a factor 7.7. Diagnostic (not a fix): Fig 2e/f with N = Nmax + 1 fails for Nmax = 2, 3 (rms 0.058, 0.025) and 6 (0.014-0.016, borderline), so Fig 2 was made with Eq. 5 as printed. Also in the code Bev and Van have ddt = 0, so repeated weekly doses would accumulate; the dose-amount coding is not in the supplement, so the Table 1 delta_V and Nmax cannot be checked against the data fit.
3. Figure 6: T0 is not given. The population parameters do not include it. Methods give a mean size of 70 mm^3, but the curves start at T about 80 and V about 94-96 (K = 1.14 would imply T0 about 83). With T0 = 70 the start is 10 mm^3 off.
4. Figure 6: the drug scale does not match Eq. 9 and alpha~_C = 1 mg/mL/mm^3. Eq. 9 gives C of about 78-93 for the daily peak in 6b (C = V at the dose, N = 1); the figure's C axis runs 0-5 and the pre-window peak is about 0.9, a factor of about 100. Scaling alpha~_C by 0.01 matches that C scale but then the drug has almost no effect on T (6b: T(49) = 89 vs 65 in the figure); the printed value kills far too much (T(49) = 33). Neither reading reproduces both C and T (diagnostic only, not fitted).
5. Eq. 7: exp(-lambda t) is written with the absolute time t (day 38 and later). Starting the resistance clock at day 38 instead kills the tumour completely, so the absolute reading is the plausible one, but it is never stated.
6. Figure 6a (no drug) with T0 inferred as 83 still gives T(60) = 148 vs about 165 in the figure and T(71) = 128 vs 137, so even the drug-free panel is not exactly reproduced from Table 1 (the figure may use other parameter values than the population estimates).
7. BIC of the monotonic model is 2,759 in the text before Table 1 and 2797 later (BIC_M).
8. The LRT statistic is printed as -Delta 2LL = -32764.61 (sign and form unclear).
9. Fig 2 caption gives delta_V as a range 0.04-0.12 and Nmax as 2-21; the five values are only in the plot legends.

## Review

Independent review, 2026-11-04. The review found the equations, the K fix and the Figure 2 digitising sound. It rejected the original grade B, which graded Figure 2 only: Figures 2 and 6 were frozen as targets together at about 09:48 and Figure 6 was relabelled "secondary" only after the first comparison showed it failing. Excluding a target after seeing it fail is not allowed, so Figure 6 is graded and the model is C, with the per-figure split (Figure 2: B, Figure 6: C) reported. Other review changes: the Figure 2b V axis was declared 100-1000 but is drawn 100-900 (frame rows map to 103 and 897), corrected; an unsupported claim that panel 6a "fails with any T0" was removed from the notes (only T0 = 70 and 83 were tried).

## What would raise the grade

Figure 2 is at B: the printed K = T0/V0 fails every Figure 2 curve, and the fix is well evidenced. Figure 6 (and with it the model grade) would need the authors' T0 for that figure, the actual C scale (or the value of alpha~_C and delta_T that were used), and the parameter set used for panel 6a. The code and values used for Fig 6 are not in the supplement.
