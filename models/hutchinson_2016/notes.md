# hutchinson_2016 working notes (2026-11-04)

Paper: Hutchinson et al., CPT:PSP 2016;5:636-645, PMC5192995. Source material: article text, figures g001-g006 and supplements s001-s006 (not redistributed). s004.txt is MLXTRAN code (the only code in the supplement).

## Equations (main text)

- Eq 1-2 control: T' = aT T (1 - T/V); V' = aV T^(2/3) V^(1/3) (beta 2/3, gamma 1/3 fixed on geometric grounds).
- Eq 3-5 treatment: T' = aT T (1 - T/(N V)); V' = aV T^(2/3) V^(1/3) - dV V; N = Nmax in (tnorm1, tnorm2), else 1.
- Eq 7-11 chemo: T' = ... - dT C T exp(-lambda t); C = aC(t) N V exp(-k mod(t,1)); aC = alpha~_C in (tCon, tCoff).
- Initial: T(38) = T0 (observed), V(38) = K T0 (the text prints K = T0/V0: typo). Time axis in days since injection, t0 = 38.

## Parameters and sources

- Fig 2 (caption): aT 0.1, aV 0.09, K 2, T0 70. Panels c, d: Nmax = 1, dV 0.04-0.12 (legend 0.04, 0.06, 0.08, 0.10, 0.12). Panels e, f: tnorm 52/62, dV 0.1, Nmax 2, 3, 6, 11, 21 (legend).
- Fig 6: Table 1 final model, vanucizumab: aT 0.109, aV 0.119, K 1.14, tnorm1 53.4, tnorm2 59.2, Nmax 6.7, dV 0.115. Chemo section: alpha~_C 1, k 0.9, dT 0.12, lambda 0.08; windows (42, 49), (54, 61), (38, 59). T0 is not stated (70 used, from the Methods).
- s004.txt: V_0 = K*X; N = Nmax*(Bev+Van)+1 in the window; dV term (dV_bev*Bev + dV_van*Van)*V; Bev/Van ddt = 0.

## Timeline

- Figure 2 and Figure 6 digitised (`reference/digitise_hutchinson.jl`, calibrated on ticks, checked on other ticks). Targets frozen 2026-11-04 about 09:48 IST, default tolerances, before any comparison.
- First comparison run: Figure 2 fails as printed (V0 = T0/K) and passes with the K fix. Figure 6 fails.
- After the first run: stray digitised points on the drug sawtooth in 6b T were removed (exclusion rectangle).

## Decisions

- Figure 6 was first labelled "secondary" (not graded) after the first comparison showed it failing, for these reasons: T0 not given; the C scale is off by about 100 against Eq 9 with the printed alpha~_C; 6a fails with T0 = 70 and with T0 = 83 (other T0 values were not tried). Review (2026-11-04): excluding a target after seeing it fail is not allowed, since Figures 2 and 6 were frozen together at about 09:48. Figure 6 is therefore graded and the model grade is C, with the per-figure split (Figure 2: B, Figure 6: C) reported in `report.md` and `score.toml`.
- No from_target fix is used: the only Figure 6 input that could be taken from the figure (T0) would be one, and it still would not make 6b to 6d pass.
- N = Nmax is kept (the diagnostic N = Nmax + 1 fails Figure 2e/f for small Nmax).

## Open questions

- Which C scale / alpha~_C and T0 did the authors use for Figure 6?
- Dose amounts in the Monolix data file (Bev/Van columns).
