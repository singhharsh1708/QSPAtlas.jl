# hutchinson_2016 working notes (2026-11-04)

Paper: Hutchinson et al., CPT:PSP 2016;5:636-645, PMC5192995. Source material: article text, figures g001-g006 and supplements s001-s006 (not redistributed). s004.txt is MLXTRAN code (the only code in the supplement).

## Equations (main text)

- Eq 1-2 control: T' = aT T (1 - T/V); V' = aV T^(2/3) V^(1/3) (beta 2/3, gamma 1/3 fixed on geometric grounds).
- Eq 3-5 treatment: T' = aT T (1 - T/(N V)); V' = aV T^(2/3) V^(1/3) - dV V; N = Nmax in (tnorm1, tnorm2), else 1.
- Eq 7-11 chemo: T' = ... - dT C T exp(-lambda t); C = aC(t) N V exp(-k mod(t,1)); aC = alpha~_C in (tCon, tCoff).
- Initial: T(38) = T0 (observed), V(38) = K T0 (the text prints K = T0/V0: typo). Time axis in days since injection, t0 = 38.

## Parameters and sources

- Fig 2 (caption): aT 0.1, aV 0.09, K 2, T0 70. Panels c, d: Nmax = 1, dV 0.04-0.12 (legend 0.04, 0.06, 0.08, 0.10, 0.12). Panels e, f: tnorm 52/62, dV 0.1, Nmax 2, 3, 6, 11, 21 (legend).
- Fig 6: Table 1 final model, vanucizumab: aT 0.109, aV 0.119, K 1.14, tnorm1 53.4, tnorm2 59.2, Nmax 6.7, dV 0.115. Chemo section: alpha~_C 1, k 0.9, dT 0.12, lambda 0.08; windows (42, 49), (54, 61), (38, 59). T0 is not stated (as published: 70 from the Methods; with fixes: 80 read from the figure, see Update 2026-11-04).
- s004.txt: V_0 = K*X; N = Nmax*(Bev+Van)+1 in the window; dV term (dV_bev*Bev + dV_van*Van)*V; Bev/Van ddt = 0.

## Timeline

- Figure 2 and Figure 6 digitised (`reference/digitise_hutchinson.jl`, calibrated on ticks, checked on other ticks). Targets frozen 2026-11-04 about 09:48 IST, default tolerances, before any comparison.
- First comparison run: Figure 2 fails as printed (V0 = T0/K) and passes with the K fix. Figure 6 fails.
- After the first run: stray digitised points on the drug sawtooth in 6b T were removed (exclusion rectangle).

## Decisions

- Figure 6 was first labelled "secondary" (not graded) after the first comparison showed it failing, for these reasons: T0 not given; the C scale is off by about 100 against Eq 9 with the printed alpha~_C; 6a fails with T0 = 70 and with T0 = 83 (other T0 values were not tried). Review (2026-11-04): excluding a target after seeing it fail is not allowed, since Figures 2 and 6 were frozen together at about 09:48. Figure 6 is therefore graded and the model grade is C, with the per-figure split (Figure 2: B, Figure 6: C) reported in `report.md` and `score.toml`.
- (Superseded 2026-11-04, see Update below.) Earlier: no from_target fix was used because T0 alone would not make 6b to 6d pass.
- N = Nmax is kept (the diagnostic N = Nmax + 1 fails Figure 2e/f for small Nmax).

## Open questions

- Which C scale / alpha~_C and T0 did the authors use for Figure 6?
- Dose amounts in the Monolix data file (Bev/Van columns).

## Update 2026-11-04: classification of the Figure 6 candidate changes

Source of the candidates: an independent rebuild from the paper alone. It matched Figure 2 to within 2e-5 of each curve's range and reproduced Figure 6 b to d within about 5 mm^3 with three choices. Rules applied: docs/src/method.md (fixes need evidence; `:from_target` only for protocol inputs read from the figure being reproduced; model parameters may never be chosen by matching the target).

1. T0 = 80 mm^3 for Figure 6: `:from_target` (B† input). A protocol input the paper does not give, read from the starting point of the curves of the figure being reproduced. Allowed, flagged. The Methods mean of 70 is a population mean of the data, not the starting point of these simulations, and is kept for the as-published run.
2. Resistance clock exp(-lambda (t - tCon)): `:ambiguity`. Eq. 7 prints exp(-lambda t) without defining the origin of t; both readings (days since inoculation, days since start of chemotherapy) are textually possible, and the text says resistance is likely to influence Fig 6(d) compared with Fig 6(c). No parameter value changes; the reading is stated. Honest caveat: the choice between the two readings was made by looking at which reproduces the figure (the rebuild did so too), so it is not independent of the target; it is acceptable as an ambiguity because it is a model-structure reading with textual support and no free number, but it earns nothing alone (with alpha~_C = 1 it makes 6b and 6d worse).
3. alpha~_C = 0.01 instead of the printed 1: NOT a fix, kept as a sensitivity run. Reasoning:
   - It is a model parameter, not a protocol input. The rules say parameters are never read off the target figure.
   - The paper does state a value (1 mg/mL/mm^3, chemotherapy section); this is not a missing value that another section fills. Replacing a stated value needs evidence of a typo or unit error from outside the figure. None found: the journal units are garbled (the preprint text s001 has the consistent forms) but give no factor of 100; the supplement has no Figure 6 code.
   - The evidence offered is the C axis of Figure 6 (peaks near 0.9 in 6b, which implies C about 0.01 N V). C is not a scored target, but it is drawn on the same figure and the same panels as the scored T and V curves, and the same value also makes the scored curves match. Reading a parameter from the figure being reproduced is exactly the circular step the dagger cannot cover (the dagger is for protocol inputs only). Also, alpha~_C and delta_T enter only as a product in the drug kill term, so the axis scale fixes the product only by assuming the plotted C is the C of Eq. 9; a different plotted scaling (for example C divided by 100) would equally fit.
   - Hence it is reported in report.md as a sensitivity run (`sensitivity(; overrides = (; alphaC_on = 0.01))` in `reproduce.jl`), not applied in the "with fixes" run, and it does not change the grade.
   - Result of the sensitivity run: with T0 = 80 and the clock from tCon, 6b T and 6b V pass (rms 0.011, 0.014), 6c and 6d are close but fail (rms 0.015 to 0.022), 6a fails (rms 0.049, 0.037). So Figure 6 would be C even if the parameter were admitted.
4. Panel 6a: not reproducible from Table 1. Peak 147.5 mm^3 at day 59.2 against about 165 in the figure (T0 = 80 does not change it). Not a fix candidate: the only ways to reach 165 (window end near 60.5, or Nmax very large) are parameter changes chosen by matching the target.
5. N = Nmax*dose + 1 with accumulating dose compartments (s004 MLXTRAN): not applied. A code-versus-text discrepancy that cannot be settled (dose records are not in the supplement); Figure 2 e/f rule out N = Nmax + 1 for small Nmax. Run as a diagnostic only (report.md, Result).

Outcome: with fixes (K, clock, T0) Figure 6 fails all 8 curves with the printed alpha~_C = 1 (6b and 6d worse than as published, 6a and 6c V better). Model grade C; figure grades from `score.toml`: Figure 6 C. Figure 2 numbers are unchanged. The library computes figure grades from the run's whole fix list, so `score.toml` labels Figure 2 B† although no Figure 2 curve uses the T0 fix; the Figure 2 grade is B on its own evidence.

## Provenance ledger (2026-11-05)

- The as-published run now lists its Figure 6 fill `T0_fig6_methods_mean` (T0 = 70, Methods mean); without a T0 fix Figure 6 throws NotRunnable.
- `resistance_clock_from_tCon` moved from the fixed run to a sensitivity run: it was adopted after a rebuild showed it matches Figure 6, and the text allows both readings, so no independent evidence justifies it on its own. With fixes Figure 6 now uses the printed exp(-lambda t); scores in report.md. Grade unchanged: C (Figure 2 B, Figure 6 C).
