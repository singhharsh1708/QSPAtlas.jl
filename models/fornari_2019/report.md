# fornari_2019: grade C

Fornari C, Oplustil O'Connor L, Pin C, Smith A, Yates JWT, Cheung SYA, Jodrell DI, Mettetal JT, Collins TA. Quantifying drug-induced bone marrow toxicity using a novel haematopoiesis systems pharmacology model. CPT Pharmacometrics Syst Pharmacol 2019;8:858-868. DOI 10.1002/psp4.12459, PMC6875710.

**Grade C.** The model cannot be run as printed (equations not closed). After documented typo/missing-equation fixes and first-principles gap fills it runs, but with the rat dose converted from the stated body weight (0.320 kg) and carboplatin MW it fails all five Figure 3 lines (rms 0.05-0.13, the model is far too mild). Of the human Figure 4 lines (two 21-day cycles, per the figure caption) only 4b neutrophils pass; 4a neutrophils, 4a platelets and 4b platelets fail. Grade follows mechanically from these scores. Working notes: `notes.md`.

## Model

Haematopoiesis as ODE chains (hours, 1e9 cells/L): free carboplatin platinum PK (2 compartments, S5) drives a Damage state (Eq 3: dDamage/dt = kdam Cp - krep Damage). Damage kills the proliferating flux of MPP, CMP, MEP (rate Slope x Damage x kprol x cells, Eqs 1, 2, 4) and of the first two (proliferative) reticulocyte and platelet transit compartments (delta x Slope_MEP, Eq 7). States: HSC source into MPP (Eq 1, zero-order kstem with feedback), MPP -> CMP and MEP (Eqs 1, 4); CMP -> 3 non-proliferative transit compartments -> neutrophils, and -> monocytes (Eqs 5, 6, 9); MEP -> 2 proliferative + 1 non-proliferative transit compartments -> reticulocytes (-> RBC) and platelets (Eqs 7-9). Feedbacks (Eqs 11-13) are powers of baseline/current ratios of mature cells: stem source from CMP and MEP, MPP maturation towards CMP from neutrophils and monocytes and towards MEP from RBC, transit proliferation from reticulocytes/RBC and platelets. All rate constants (ktr, kprol, kstem, a, kcirc_Ret) are derived from baselines, MTTs, kcirc and lambda1..5 by the steady-state algebra of S6/S7; the derived steady state is exact (relative drift below 1e-15 with no drug). 23 states.

Files: `model.jl` (parameters with sources, `derive_rates`, `haematopoiesis_system`), `protocol.jl` (fixes, dosing, targets, `benchmark_setups()`), `reproduce.jl`.

## Targets

Frozen 2026-11-04 (about 09:49 IST) before any simulation of this model was run; default tolerance (RMS 0.015, P95 0.03 of the axis box), no overrides. All kind `curve`: the red "mean model prediction" lines of the authors, digitised from the 709 px JPEGs by a red-pixel mask (median row per column), `reference/digitise.toml`, `reference/digitise_run.jl`.

| id | figure panel | points | note |
|---|---|---|---|
| fig3_MPP, fig3_CMP, fig3_MEP | Fig 3 bone marrow, rat 40 mg/kg x8 q14d | 131, 102, 128 | log y |
| fig3_Neut, fig3_Plt | Fig 3 blood | 54, 129 | neutrophil line is faint under data symbols, first ~20 d and some stretches not traced |
| fig4a_neut, fig4a_plt | Fig 4a, human, rat slope values | 60 each | red line hidden under the blue Schmitt line for the first ~3-4 d (and where overlapped) |
| fig4b_neut, fig4b_plt | Fig 4b, human, IC50-scaled slopes (Table 2) | 60 each | |

Not targeted: Fig 3 monocytes, reticulocytes, RBC (simulated and plotted, not graded); observation symbols and percentile bands; Fig 4 blue Schmitt/Friberg lines, bands and bar chart (c); Supplement Figs S1, S2, S4-S7.

Calibration: x from tick-label/tick centroids, log y from tick-label glyph centroids (-0.5 px, template-checked on the minor ticks) for Fig 3 and from the frame lines (axis limits are decades) for Fig 4. Verification (`figures/digitise_check.png`, local only because it contains the paper's figure; points on the lines): Fig 3 x tick "60" within 0.45 px of the interpolated position, MEP 10^1 label within 0.23 px of the two-point interpolation; Fig 4 day-20 tick within 0.2 px; Fig 4 grade 3/4 threshold lines (1, 0.5; 50, 25 x10^9/L) at 129.5/145.0 px, 326.0/342.0 px against 129.5/145.2, 326.5/342.2 predicted from the calibration. Residual calibration uncertainty is about 0.5 px, 0.4 % of the box.

## Result

As published: **error, model not runnable** ("Eq 1 contains fdbk_mat_MPP ... Eq 9 omits the RBC equation ..."). With fixes (all from `score.toml`):

| curve | rms | p95 | pass |
|---|---|---|---|
| fig3_MPP | 0.0676 | 0.1595 | no |
| fig3_CMP | 0.0917 | 0.1229 | no |
| fig3_MEP | 0.1291 | 0.2846 | no |
| fig3_Neut | 0.1336 | 0.1597 | no |
| fig3_Plt | 0.0465 | 0.0880 | no |
| fig4a_neut | 0.0693 | 0.1036 | no |
| fig4a_plt | 0.0323 | 0.0540 | no |
| fig4b_neut | 0.0099 | 0.0179 | yes |
| fig4b_plt | 0.0275 | 0.0526 | no |

Figures: `figures/comparison.png`, `figures/fig3_regenerated.png`, `figures/fig4_regenerated.png`. At 0.862 umol per mg/kg the rat model depletes the progenitors far less than the paper's lines (MEP nadir ~20 vs ~5 x10^9/L at 30 mg/kg single dose).

### Sensitivity (diagnostic only, `sensitivity.jl`; not fitted, does not change the grade)

- Rat dose conversion. Main (0.862 umol per mg/kg) fails all Fig 3 lines. The Fig S3 PK curve implies 3.90 umol per mg/kg (117 umol for 30 mg/kg, a ~1.45 kg rat at MW 371.25), i.e. the paper's dose/PK units are internally inconsistent or unstated. With 3.90 the rms would be 0.010/0.032/0.010/0.036/0.008 (MPP/CMP/MEP/Neut/Plt: three of five pass); 4.41 (from the Fig S3 Damage curve) gives 0.013/0.041/0.012/0.048/0.008; 3.0 would pass all five. These values are fitted to figure curves (parameter estimation) and are not used for the grade.
- First dose at t = 0 instead of 24 h: same result.
- Fig 4 (two cycles): bolus vs infusion, a third dose, Neut0/Mono0/Ret0/RBC0 anywhere in the printed ranges and CL 9 L/h do not change the pass/fail pattern much (4b neutrophils 0.009-0.014, 4a neutrophils 0.064-0.071). Eq 10 slopes in (b) give 4b neutrophils 0.014, 4b platelets 0.011 (pass). A single dose (S11/S12 wording) fails all four (0.040-0.069).

## Fixes (used in the "with fixes" run)

| id | kind | what | evidence |
|---|---|---|---|
| eq1_maturation | typo | Eq 1 maturation terms -ktr_CMP fdbk_mat_MPP - ktr_MEP fdbk_mat_MPP replaced by -ktr_CMP fdbk_mat_CMP MPP - ktr_MEP fdbk_mat_MEP MPP | fdbk_mat_MPP defined nowhere; Eq 4 gain terms are the same fluxes; S4 kstem needs outflow proportional to MPP; text says Eq 12 feedbacks enter Eqs 1 and 4 |
| rbc_equation | missing | add dRBC/dt = kcirc_Ret Ret - kcirc_RBC RBC | Eq 9 prints four equations; Table 1 kcirc_RBC; S4 kcirc_Ret = kcirc_RBC RBC0/Ret0; Eqs 12-13 use RBC |
| plt_transit_rate | typo | platelet T1 loss uses a_Plt, not a_Ret | Eq 8, T2_Plt equation, S3/S4 steady state |
| dose_time | ambiguity | rat doses at t = 24 h + 14 d k (k = 0..7); human first dose at t = 0 | Fig S3 hour axis starts at 24 h with the dose at the left edge; Table S1 "day 1"; Fig 3 shows eight cycles |
| dose_scale | missing | rat dose enters Cen as d x 0.320 kg x 1000 / 371.25 umol (0.862 umol per mg/kg) | S1: rats "a minimum 320g of body weight at the start of dosing" (only a minimum stated, 0.320 kg used); carboplatin MW 371.25 g/mol (PubChem CID 426756). First principles, not fitted. The Fig S3 PK curve implies 3.9 umol per mg/kg (~1.45 kg), diagnostic only |
| human_schedule | ambiguity | two 750 mg doses, day 0 and day 21 | Fig 4 caption: "Two cycles of 21 days of carboplatin were simulated. Doses were calculated with the Calvert formula"; S11/S12 say "a single dose" |
| human_pk | missing | Zandvliet 2008 two-compartment PK: CL 0.76 x 7.5 + 1.5 = 7.2 L/h, V1 15.5, Q 3.46, V2 9.86; 750 mg (S11 Calvert) per dose as 30 min infusion | S11 cites the model but gives no values; values from Zandvliet Table 2 (PMC2561111); GFR = CrCl and 30 min infusion assumed |
| human_baselines | ambiguity | midpoints of printed ranges: Neut0 4.5, Mono0 5.1, Ret0 77.5, RBC0 5000 | Table 1 prints ranges; rule declared before any simulation (not read off Fig 4) |

No parameter was fitted to a target or to a figure curve in the graded run. Table 2 human slopes are used as printed even though they disagree with Eq 10 (below).

## Issues found in the paper

1. Eq 1: `fdbk_mat_MPP` undefined, missing factor MPP on both maturation terms (main text, Eq 1).
2. Eq 4: second line is labelled dCMP/dt, is dMEP/dt.
3. Eq 7: platelet T1 equation uses a_Ret and the label "T1_Ptt".
4. Eq 9: no RBC equation. Table 1 and S4 imply it. Table 1 does not list the reticulocyte kcirc (derived).
5. Supplementary model code (s013): the PDF stops after the parameter mapping (page 4 ends at `k12_Carbo = QC/Vcen_Carbo;`); ODE right-hand sides, dosing and solver call are absent.
6. Dose-to-amount conversion (body weight, molecular weight) is not given anywhere. The physical value (0.32 kg, 371.25 g/mol: 0.862 umol per mg/kg) is 4.5x smaller than what Fig S3 needs (3.9).
7. Fig S3: the Damage curve (peak 4.77 for 30 mg/kg) is about 13 % higher than the printed Damage ODE gives for the PK curve in the same figure (peak 4.13) and decays faster (0.0154/h vs 0.0127/h) although kdam = krep = 0.017.
8. Table S3 gives the units of QC and CL as 1/h; the code says L/h (the values 0.06 and 0.42 only work as L/h).
9. Table S1 (PDF text) lists doses on days 1-85 for study E (seven doses) while Fig 3 and the text show eight cycles.
10. Table 1 human baselines for Neut, Mono, Ret, RBC are literature ranges, not values (Mono 0.2-10 spans a factor 50). Human CMP0 20.9 and MEP0 15 differ slightly from 16 x MPP0 and 11 x MPP0 (S9).
11. Table 2 human slopes are not Eq 10 applied to Table S2: Eq 10 gives CMP 1.03 and MEP 1.21 (Slope_R x IC50_R/IC50_H with the printed IC50 values), Table 2 prints 0.57 and 0.66, i.e. a consistent extra factor of about 0.55 on all three slopes (0.79 for MPP as well, although Table S2 has no MPP IC50 and the CMP ratio seems to have been used). The extra factor is unexplained.
12. Fig 4 contradicts the slope values: Table 2 human slopes are 2.6-3.3x smaller than the rat slopes, so the (b) predictions should be less toxic than (a), but in the paper's Fig 4 and bar chart (c) (b) is clearly more toxic (neutrophil nadir 2.6 vs flat in (a)). Our model gives the expected ordering. Either panels/labels differ from the description or the scaling is not what Eq 10/Table 2 say.
13. S11/S12 say a single carboplatin dose was simulated, but the Fig 4 caption says two 21-day cycles and all Fig 4 lines, including the Schmitt model, show a second trough near day 35 (period about 21 d). The caption is used (fix human_schedule).
14. Figure 3 and 4 red lines are "average model predictions" of 1000 virtual subjects with log-additive residual error (Table S4 sigma 0.17-0.19 for neutrophils/CMP), but no inter-individual variability model is given, so only a typical-value simulation can be reproduced. Residual error of that size shifts the mean by only about 1.5-1.8 % and cannot explain a large offset.

## What would raise the grade

- The complete supplementary ODE file (right-hand sides, dose conversion and body weight, dosing and time grid) and the IIV/residual model behind the "average" lines. The rat dose conversion is decisive: from stated weight and MW the model is far too mild, while the Fig S3 PK curve implies a 4.5x larger amount; the paper does not say which was used for the PD fits.
- For Fig 4: the human baselines actually used (single values), the slopes used in panel (a), and an explanation of the extra factor 0.55 between Eq 10 and Table 2 (item 11) and of why (b) is more toxic than (a) (item 12).

## Review

Independent review, 2026-11-04. The review argued for D because the original dose conversion (3.90 umol per mg/kg) was fitted to the Fig S3 PK curve, which is parameter estimation, and noted that the Fig 4 caption specifies two 21-day cycles and that the report mis-called that schedule circular and over-claimed (line 81, line 86). Final decision: the fitted value is removed and replaced by a first-principles conversion under the documented gap rule (stated minimum body weight 0.320 kg, MW 371.25), and the grade follows from the resulting scores.

Changes this round: (1) `dose_scale` now 0.320 kg x 1000 / 371.25 (S1 weight, PubChem MW); the S3-implied 3.9 is kept only as an ungraded diagnostic (`dose_scale_from_S3.jl`, `sensitivity.jl`). (2) New fix `human_schedule` (:ambiguity): two doses, day 0 and 21, from the Fig 4 caption; the claim that this was circular is removed. (3) Grade recomputed from the new scores (C; Fig 3 now fails all five lines, Fig 4 passes only 4b neutrophils). (4) "No parameter was re-estimated" reworded, residual-error statement corrected (1.5-1.8 % mean shift), unsupported "likely pass" claim removed. Targets and tolerances were not changed.
