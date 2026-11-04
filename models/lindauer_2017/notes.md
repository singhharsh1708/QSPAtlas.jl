# lindauer_2017 research notes

Lindauer A et al. CPT Pharmacometrics Syst Pharmacol 2017;6:11-20, DOI 10.1002/psp4.12130, PMC5270293.

Licence of the article and supplement: CC BY-NC-ND 4.0 (XML `<permissions>`). The mouse dataset (s007.csv) is NOT copied into the repository.

## Target decision (written 2026-11-04 09:5x IST, before any simulation was run)

Fig 3a (reduction of tumour diameter vs dose, six scenarios, bootstrap bands):

- The plot is a dose-response at (presumably) 6 months, not a time course as the caption says ("over 6 months").
- Methods: "Human PK parameters were taken directly from the bootstrap results ... Uncertainty from the mouse parameters ... propagated". The caption says the bands come from the bootstrap analysis. Nothing says the central line is the nominal (typical-value) simulation; it is most likely the median of the bootstrap replicates. Fig S1 states that the *typical response (no uncertainty)* is used only for the sensitivity analyses (Fig 4).
- Cross-check from the paper's own numbers: the Fig 4 base case at 5 mg/kg is a log10 ratio of -1.29, i.e. a diameter reduction of 1 - 10^(-1.29/3) = 63 %; at 0.5 mg/kg -0.39, i.e. 26 %. The Fig 3a lines read roughly 55-80 % at 5 mg/kg and 40-55 % at 0.5 mg/kg, i.e. different from a nominal simulation of any one scenario.
- Therefore the Fig 3a central lines are NOT used as grade-bearing targets (they are not nominal simulations). Not digitised.

Fig 2b: the black line is the median of OBSERVED RO in bins; the model line is a grey 90 % band from simulation with uncertainty. No deterministic typical line. Not a target.

Fig 2a (VPC): observed median and band. Not a target.

Fig 4 (tornado, 5x and 0.2x of 16 parameters, 5 and 0.5 mg/kg): deterministic one-at-a-time perturbations around the base case ("typical response (no uncertainty)", Fig S1). Grade-bearing targets: the 64 bar tips and the two caption base values (-1.29, -0.39). Encoded as kind: points when the targets were written (the reference CSVs use kind: bars, see `report.md`), one CSV per panel x multiplier (16 bars, in plotted order) plus `fig4_base_case.csv` (2 values). The parameter order is in the CSV headers and `reference/digitise.toml`.

Tolerance: default (rms 0.015, p95 0.03 of the axis box). Left panel box = 3.5 log units (0.05 log units rms), right panel box = 2.5 (0.04). Digitisation error is about 1 px = 0.015 log units (left) / 0.010 (right).

Scorer note: the scorer measures point-to-polyline distance, so the simulated bars are returned as a step function (plateau per slot of width 0.5 in x, box 0..8), so that the distance equals the vertical error up to half a slot (0.031 of the box). Exact vertical errors are also reported in `report.md` and below.

Freeze: `reference/*.csv` and `reference/digitise.toml` written 2026-11-04 09:47 IST. Overlay check (local only, it contains the paper's figure): calibration check ticks -1 and 0 agree to < 0.3 px; the caption base values map to 158.75 / 485.56 px against axis lines at 159 / 486.

## Protocol choices not stated in the paper (declared before running anything)

The paper does not say, for the Fig 4 base case: which of the six scenarios is used (Fig S1: "Take parameter set for one scenario (melanoma, 1)"), the dosing schedule, the body weight, the human V_blood, and whether Table 1 (rounded) or Table S2/S3 values are used.

As-published run (declared a priori):

- scenario 1 of Table S3 = fast growth / allometric scaling (L0 0.0088/day, SLtg 2.575E-06, W0 41.5 mL = 41500 mm3, exponential growth only: L1 not used);
- Q3W bolus doses at 0, 21, ..., 168 days (9 doses), response = V(182 d)/V(0), 5 and 0.5 mg/kg x 70 kg (70 kg is the human body weight stated for the allometric scaling);
- human PK = Table 1 typical values ("Actual human PK used"), human binding / FcRn constants from Table 1, species-invariant parameters from the Table S2 estimates (full precision: gamma 2.28, Tmulti 4.32, ...);
- V_blood = 1400 uL (the only value printed anywhere: s006 control stream; Table 1 has no V_blood row);
- KdegPD1 human = 0.00246 (Table 1 as printed; the stated allometric formula gives 0.00252).

Candidate fixes (menu fixed in advance; selected only on the caption base values -1.29 / -0.39, so that the 64 bars act as an independent check):

- F-Vblood: human V_blood = 5 L (physiological; 1400 uL is a 20 g mouse blood volume);
- F-scenario: other rows of Table S3;
- F-schedule: Q2W.

## Results (2026-11-04)

Identification on the two caption base values only (V_blood = 5 L; log10 volume ratio at 26 wk, 5 mg/kg / 0.5 mg/kg; targets -1.29 / -0.39):

| scenario | Q3W | Q2W |
|---|---|---|
| slow_allom | -1.707 / -0.731 | -1.873 / -0.882 |
| fast_growth | -1.441 / -0.456 | -1.613 / -0.621 |
| fast_allom | -1.541 / -0.554 | -1.712 / -0.719 |
| med_growth | -1.302 / -0.401 | -1.424 / -0.558 |
| slow_growth | -0.692 / -0.318 | -0.693 / -0.460 |
| med_allom | -1.653 / -0.672 | -1.822 / -0.829 |

As published (scenario 1, V_blood 1400 uL): -6.0 / -4.7. With V_blood 5 L, med_growth, Q3W: -1.302 / -0.401 (caption -1.29 / -0.39); with KdegPD1 by the allometric formula: -1.294 / -0.391.

Sensitivity of that base case to rounding and code choices (5 / 0.5 mg/kg): base -1.302 / -0.401; Table 1 rounding (gamma 2.3, Tmulti 4.3) -1.340 / -0.418; KdegPD1 allometric formula -1.294 / -0.391; SLtg from the L0 ratio (0.0036/0.113 * 1.98e-5) -1.306 / -0.402; the control-stream RO <= 0.01 guard changes nothing (-1.302 / -0.401).

Solver check on the same base case: Rodas5P (abstol 1e-14, reltol 1e-9), Rodas5P (1e-16 / 1e-11) and FBDF give identical log ratios to 6 digits (-1.301768).

Mouse growth check: the mouse vehicle run (170 mm3 start) equals a plain ODE of the s006 growth structure dTV/dt = L0 TV / (1 + (L0/L1 TV)^20)^(1/20) to 4 decimals (day 5: 299.1061; day 10: 526.2616; day 15: 925.9298; day 60, transition active: 9101.857).

Fixes: human_V_blood_5L (:missing), KdegPD1_allometric (:typo), scenario_med_growth (:from_target). Grade C: 58 of 64 bars within 0.049 log10; six bars off (gamma x5 two bars, the L row four bars, which look like swapped colours).

Bars that are off (log10 ratio, paper vs simulation with fixes):

| panel | bar | paper | simulated | difference |
|---|---|---|---|---|
| 5 mg/kg x5 | gamma | -2.768 | -7.027 | -4.260 |
| 5 mg/kg x5 | L | -0.822 | -1.505 | -0.684 |
| 5 mg/kg x0.2 | L | -1.540 | -0.833 | +0.707 |
| 0.5 mg/kg x5 | gamma | -1.802 | -2.673 | -0.871 |
| 0.5 mg/kg x5 | L | 0.000 | -0.978 | -0.978 |
| 0.5 mg/kg x0.2 | L | -0.969 | -0.005 | +0.964 |

The other 58 bars: rms 0.0049 of the axis box (0.017 log10), maximum 0.0139 of the box (0.049 log10). Base case: paper -1.290 / -0.390, simulated -1.294 / -0.391. Per-curve scores are in the table of `report.md`.

gamma x5 solver diagnostics: both runs finish successfully (retcode Success at t = 182 d). 5 mg/kg: TV at the end 0.0039 mm3 (log ratio -7.03); 0.5 mg/kg: 87.8 mm3 (-2.67). The RO <= 0.01 % linear branch of s006 makes no difference to these values (-7.0273 / -2.6732 with and without the guard).

L row with the colours swapped (diagnostic only, not used for any fix): simulated L at 5 mg/kg is x5 = -1.505, x0.2 = -0.844; at 0.5 mg/kg x5 = -0.987, x0.2 = -0.011. These are close to the paper's bars of the opposite multiplier (5 mg/kg: x0.2 -1.540, x5 -0.822; 0.5 mg/kg: x0.2 -0.969, x5 0.000).

Other findings:

- The equation cross-check findings are in `report.md` (supplement typos for Ceb / FcRn, truncated central-compartment equation, s006 is the initial-estimate stream, RO > 0.01 branch).
- s007.csv has no observations (DV all "."), so it cannot feed a Turing recalibration.
- The Fig 3a nominal dose-response was not simulated (not a target).

Self-check (2026-11-04): an independent rebuild from the paper alone was compared with these results; see the Self-check section in `report.md`.
