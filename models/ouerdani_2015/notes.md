# ouerdani_2015 working notes (2026-11-04)

Paper: Ouerdani et al. CPT:PSP 2015;4:660-668, PMC4716582. Target: Figure 4 (typical P solid, K dashed; left mice, right patients).

## Equations (Eq. 1-3, text)

dP/dt = lambda P (1 - P/K) - alpha exp(-delta t) P ; dK/dt = b P^n - gamma K ; alpha = alpha0 AUC^beta_alpha ; gamma = gamma0 AUC^beta_gamma. Eq. 4: P = K = (gamma/b)^(1/(n-1)).

NONMEM code in the supplement: s004.txt = preclinical (DADT(2) = B*A(1), n = 1; THETA(8) beta_alpha FIX 0), s001.txt = clinical (DADT(2) = B*A(1)**0.5). Time is in days. The `$THETA` values in the code are initial estimates (not Table 1). P0 is a data column (VB / BASE).

## Parameters (Table 1, typical)

| | lambda | K0 | b | gamma0 | beta_gamma | alpha0 | beta_alpha | delta | n |
|---|---|---|---|---|---|---|---|---|---|
| Mice | 0.166 | 543 | 0.0183 | 0.007 | 0.332 | 0.251 | 0 (fixed) | 0.196 | 1 |
| Patients | 0.0021 | 329 | 0.0392 | 0.0023 | 0.142 | 0.0032 | 0.125 | 0.0153 | 0.5 |

Mice: AUC 220.2 / 656.8 / 1140.8 for 10 / 30 / 100 mg/kg. Patients: mean AUC 771.6 (range 629.4-802.4, mean dose 727 mg); Supplementary Table 2 uses a fixed 800.

## Open questions and decisions

- P0 (typical): not stated. Read from Figure 4 at t = 0: mice 150 (coincides with a tick), patients 66.8 (top of the curve start on the y axis, pixel row about 29). Declared as a fix (:from_target; circular, see report.md).
- Mouse dose: not stated. All three tried; 100 mg/kg matches (rms 0.001), 30 mg/kg gives 0.074, 10 mg/kg gives 0.27. Declared as a fix (:from_target: chosen by matching the figure).
- Patient AUC: 771.6 (stated) passes (rms 0.009 / 0.011 for P / K); 800 and 802.4 fit K about six times tighter (0.002). 771.6 was kept and the alternatives are reported.
- Days per month on the clinical axis: 365.25/12 used; 30 also passes, 28 fails.
- Treatment is continuous in the model (no stop at 24 d for mice).

## Timeline

- 09:42 Both panels digitised (rule-based, black curves); overlay checked (the overlay image contains the paper's figure and is kept locally only).
- 09:44 Targets frozen: `reference/*.csv`, `reference/digitise.toml`, default tolerances (RMS 0.015 / P95 0.03) in `protocol.jl`. No simulation had been compared yet.
- After the freeze: the unstated inputs (P0, mouse dose, patient exposure) were tried, the fixes above chosen, the sensitivity runs summarised in `report.md` made.

Self-check (2026-11-04): an independent rebuild from the paper alone was compared with these results; see the Self-check section in `report.md`.

## Provenance ledger (2026-11-05)

The as-published run now lists its five fills as Fixes (`p0_text_guess`, `mouse_dose_text_guess`, `patient_auc_mean`, `days_per_month`, `continuous_treatment`); the last three are kept in the fixed run. Without a P0, dose or month fix the simulation throws NotRunnable. No number changed; grade B†, strict D.
