# stein_2019: grade C

Stein AM, Grupp SA, Levine JE, et al. Tisagenlecleucel Model-Based Cellular Kinetic Analysis of Chimeric Antigen Receptor-T Cells. CPT Pharmacometrics Syst Pharmacol 2019;8:285-295. DOI 10.1002/psp4.12388, PMC6539725.

Grade C: the model runs and, simulated under the most literal reading of the paper's design (90 patients observed at the scheduled visit days, Table 1 values, residual error), reproduces the right-panel 10th and 50th percentile bands (days 0-30) but not the 90th percentile at the peak nor the persistence phase of the left panel (days 30-185). Reading the design more densely (dense grid, band bin nodes; the "with fixes" run) does not change this. The paper does not state the VPC design (who is sampled when, bins, BLQ and quantile handling, base vs full covariate model); those choices are listed in FIXES with their evidence and are not published facts, so the result is "fails under this reading", not a definitive failure of the authors' simulation. The misses are mostly several times the digitising uncertainty and far above Monte Carlo error. Grade C means "runs, does not match", not "cannot run".

## Model

Cellular kinetic model of CAR-T transgene copies/ug DNA (Eq. 1): exponential expansion C = R0 exp(rho t) up to Tmax (rho = ln(foldx)/Tmax, R0 = Cmax/foldx), then C = Cmax(1-FB) exp(-alpha (t-Tmax)) + Cmax FB exp(-beta (t-Tmax)). Implemented as the paper's own compartment model (Figure 2, supplement lines 6-31): effector cells E (grow at rho, then decline at alpha) and memory cells M (fed at k = FB(alpha-beta), decline at beta), C = E + M, MTK equations with `ifelse` on t against Tmax. The tocilizumab/corticosteroid effect (Eq. 2, Figure 1b) is in the model (T1, T2, F1, F2) and used in a sensitivity run only. The ODE matches the closed form of Eq. 1/2 and the supplement's Monolix code to 1.4e-9 relative error (Vern9, tolerance 1e-10; `check_model`).

Parameters (Table 1, sources in `model.jl`): foldx 3900, Tmax 9.3 d, Cmax 24000 copies/ug, alpha 0.16 /d, FB 0.0079, beta 0.0032 /d (Ftoci 1.2, Fster 1.0 not used in the graded run). All six log-normal with omega foldx 2.4, Tmax 0.38, Cmax 0.65, alpha 0.91, FB 0.8, beta 0.86; residual: constant error in log space, a = 0.56. Reading of Table 1: omega is a standard deviation of eta (Monolix convention; supplement line 77 "eta with variance omega^2"; the table gives no CV% or variance label), diagonal Omega, a is the sd of the log-scale residual (`ExponentialError(0.56)`). Typical values are consistent with the text (doubling time 0.78 d, half-lives 4.3 and 217-220 d).

## Targets

Figure 3a (PSP4-8-285-g003.jpg), a Monolix VPC in two panels (0-380 d, and a zoom on 0-30 d), log y axis, LLOQ 50 copies/ug. Caption: blue lines are the 10/50/90th percentiles of the DATA (empirical, not graded); the shaded areas are the intervals of the percentiles from the model (blue 10th and 90th, pink 50th) and are the authors' simulation. Graded: six `kind: band` references (`reference/fig3a_{left,right}_{p10,p50,p90}.csv`, `x,lo,hi`), scored with the `:band` comparison (error 0 when the simulated percentile lies inside the band at that x, else the vertical distance to the nearest edge in axis-box units). Default tolerance (RMS 0.015, P95 0.03 of the box), which here is the allowance for digitising the band edges. Scored ranges: right panel p90 4.5-29.5 d, p50 3.4-4.4 and 6.4-29.5 d (red outlier marks cover 4.4-6.4 d), p10 6.4-29.5 d; left panel 30-185 d (before 30 d the band is under the data cloud and covered by the right panel; after ~187 d the band is a flat extension of the last bin). Figure 3b (individual fits) is not a simulation summary and is not a target.

Digitising: colour-run tracing of the band pixels (`reference/trace_fig3a.jl`, calibration and checks in `digitise.toml`, overlay `figures/digitise_check.png` is local only). y calibration uses the printed tick labels; the "100,000" label is displaced in the figure, and the grid lines do not coincide with the labels, so the y calibration has about 3 px (0.05-0.06 decade) of uncertainty; the LLOQ position (50 copies/ug) predicted from the labels agrees with the lowest quantifiable dots and the first simulated-BLQ dots within one dot radius. Digitising correction made after review (2026-11-04, after the first comparison): two points of the right p90 band (x = 8.8772 and 9.9825 d) had been traced from the lower (10th percentile) blue band where stacked dots split the upper band (found in review). They were dropped by a rule stated in `trace_fig3a.jl` and `digitise.toml` (a point whose whole interval lies away from its neighbours' interval), applied to all six targets; no other point was removed. Tolerance unchanged. Right p90 score before / after: rms 0.0515 / 0.0235 (with fixes), a fail either way.

Timeline: targets, tolerance and a first design (dense grid for all six bands, left-panel nodes read from the target) were frozen 2026-11-04 14:50 IST, before any simulation was compared with the bands (`notes.md`). The first comparison (graded design of that time; 400 replicates) gave rms left p10 0.107, p50 0.058, p90 0.033, right p10 0.000, p50 0.005, p90 0.052 (C). After review the graded runs were restructured (below) and the two misassigned reference points corrected; targets and tolerance were not otherwise changed.

## Simulation design

Two runs share one simulation (base model, 90 patients, residual error included, BLQ values kept, quantile type 7, every patient observed at every time):
- "as published" = the most literal reading of Methods "Sample analysis": patients observed at the SCHEDULED visits (days 1, 2, 7, 11, 14, 21, 28, months 3 and 6 = 91.3 and 182.6 d); the simulated percentile at those times is joined by straight segments in (day, log10 copies) and compared with the band at each scored x (for both panels).
- "with fixes" = a dense grid (0.25 d steps on (0, 30]) for the right panel, and for the left panel the percentile at three bin nodes (28, 68, 190 d) joined by chords (how Monolix draws bins).
- Choices not stated by the paper that BOTH runs need are recorded as Fixes (see below); the as-published run cannot be defined without them, so they are in the fixes list for transparency but do not differ between the two runs. Only `dense_grid` and `left_bin_nodes` differ.
- The CAR-T dose is not an input (no dose-exposure relationship; the infusion enters through R0). Time t = days after infusion.
- 800 replicates in blocks of 50 (seed 20260411 + block), 4 threads. `reproduce()` takes about 3.1 min in a fresh session (186 s including compilation); the earlier 400 replicates had a maximum MC SE of 0.0034 box, above the 0.003 rule set beforehand, so 800 were used. The chunked percentile code is identical to QSPAtlas `vpc` (same cohort and seed: maximum difference 0.0, `check_vpc_consistency`).
- Monte Carlo SE of the compared median (bootstrap over replicates, as in `vpc`): at most 0.0023 of the axis box at any scored point in either run (0.010 decades; mean 0.001-0.0017), below the 0.003 rule, 6 times below the 0.015 RMS tolerance and well below the discrepancies below. The 0.9 interval in `figures/fig3.png` is only a display (the paper does not state the level).

## Result

From `score.toml`. Grade: as published fails, with fixes fails, so C.

| band | scored x | as published rms / p95 | pass | with fixes rms / p95 | pass |
|---|---|---|---|---|---|
| right p10 (blue) | 6.4-29.5 d | 0.000 / 0.000 | yes | 0.000 / 0.000 | yes |
| right p50 (pink) | 3.4-29.5 d | 0.010 / 0.019 | yes | 0.006 / 0.013 | yes |
| right p90 (blue) | 4.5-29.5 d | 0.026 / 0.041 | no | 0.024 / 0.038 | no |
| left p10 (blue) | 30-185 d | 0.108 / 0.157 | no | 0.108 / 0.155 | no |
| left p50 (pink) | 30-185 d | 0.059 / 0.078 | no | 0.060 / 0.076 | no |
| left p90 (blue) | 30-185 d | 0.047 / 0.064 | no | 0.035 / 0.046 | no |

Where it misses (simulated median vs authors' band, with-fixes run): right p90 at the peak (day 12): 40,000 vs 63,000-150,000 copies/ug; left panel at day 182: p10 17 vs 73-120, p50 104 vs 223-379, p90 566 vs 914-3,230; at the 68-d node p50 241 vs 464-729. The authors' simulated bands are wider in the upper tail and much higher in the persistence phase than the Table 1 model gives. The authors' p10 never falls below the LLOQ in the left panel, while the Table 1 model has roughly a quarter of the late values below 50.

Sensitivity runs (not graded, none selected by evidence; 800 replicates each; `reproduce(sensitivity = true)`). Pass count out of six bands, as published / with fixes:
- omega read as variances (sd = sqrt): 0 / 1 of 6 (left p90 passes in the with-fixes read; the right panel fails).
- BLQ censored to 50: 2 / 2 (right p10, p50; left p10 improves 0.108 to 0.052, still fails).
- comedication model, assumed timing (36 % tocilizumab, 26 % steroids, Table 2 medians/ranges, Ftoci 1.2, Fster 1.0): 4 / 4 (left p90 and all right bands pass; left p10 0.093 and p50 0.036 fail).
- residual sd halved: 2 / 2.
- left panel on the dense percentile curve instead of nodes: all three left bands fail in every variant (rms 0.02-0.13).
None gets all six bands; the comedication model is the closest, but its design is an assumption (the paper gives only marginal proportions, medians and ranges), so it is a lead for a future fix, not a fix.

Local diagnostic (skipped when `PSP4-8-285-s002.csv` is absent): the paper's Data set S1 is described in the paper as a "simulated dataset based on the model". It has the real design (90 patients, 836 rows; 35.6 % tocilizumab at median 5.7 d, 25.6 % steroids at median 7.5 d, the Table 2 ranges) but its values are far above Figure 3a (median ~1.2e5 and p90 ~8e5 at days 9-12, maximum 2.4e8 copies/ug), so it is not patient data, cannot validate the model, and is not a target. `figures/digitise_check_s1_diagnostic.png` (local only) shows its percentiles against the simulated bands.

## Fixes

Seven design choices that the paper does not state (details and evidence in `protocol.jl`): `all_patients_all_times` (:missing), `residual_included` (:ambiguity), `blq_kept` (:ambiguity; the caption's red asterisks show simulated values are plotted, not how they enter the percentiles), `quantile_type7` (:missing), `base_model_no_comedication` (:ambiguity; Results links the base model to Figure 3 and says the final model is the base model, Methods selects the full model, and nothing shows that the fitted Ftoci = 1.2 was replaced by 1, so this is unresolved), `dense_grid` (:missing), `left_bin_nodes` (:from_target: nodes read from the band in the target figure). The first five apply to both runs (the as-published run cannot be defined without them); the last two only to "with fixes". Model parameters were never changed. Because `left_bin_nodes` is :from_target, a match could only have earned B-dagger; there was no match.

## Issues found in the paper

- Table 1 random effects have no label saying SD, CV% or variance; the footnote formula for eta shrinkage ("(1 - var(eta))/omega^2") is also garbled. Read as SD here.
- Supplement line 10 prints dM/dt = alpha E - beta M; consistency with C2 = Cmax k/(alpha-beta) = Cmax FB requires dM/dt = k E - beta M.
- Supplement Eq. 6 gives dichotomous covariate effects as z = theta^x (a multiplier), while Table 1 labels them "log Cmax covariate effect" (additive in log). The Methods say the full covariate model is the final model, the Results that the final model is the base model (Table 1 has both the covariate effects and the typical values for the reference category). Which model produced Figure 3a is not stated.
- Fig 1b / Eq. 2: F < 1 is described as a decrease of the expansion rate, but the estimate is Ftoci = 1.2 (> 1). Under the supplement's CmaxAdjust, Ftoci = 1.2 for a patient given tocilizumab at the median 5.7 d (Tmax 9.3 d) raises the peak by exp(0.2 rho (Tmax - T1)), about 1.9 times, which makes "Cmax" a typical value only for patients without comedication.
- The VPC design is not described: bins, prediction-interval level, which patients and times were simulated, how BLQ values enter the VPC (only "red asterisks denote simulated data points for BLQ"). The comedication regressors of the VPC simulation are not stated. Data set S1 is a simulated data set whose values are inconsistent with Figure 3a; the real data are not available.
- Figure 3a: the y axis tick labels and grid lines do not coincide (the 100,000 label is displaced); the bands left of the first bin node and right of the last bin node are flat extensions.
- Caption says "confidence intervals" for the shaded areas, the legend says "Prediction interval".
- Time origin: Methods call day 1 the infusion day and list sampling at "days 1, 2, 7..."; the figure shows the first samples at about 0.7-1 d after infusion, so day N is read as t = N d.

## Self-check (2026-11-04)

An independent rebuild from the paper and S1, with separate code (plain OrdinaryDiffEq, no shared code). The data set S1 (s002) was read at run time for its design only (ids, times, tocilizumab and steroid regressors); the values were ignored.

Structure: the closed form and the ODE agree to 3e-13 in log. The derived values match the text (doubling time 0.78 d, ln2/alpha 4.33 d, ln2/beta 217 d against "220").

Same unstated choices, made independently: omega read as SD (the footnote uses omega^2 for shrinkage, Monolix reports omega as an SD, and the variance reading puts the simulated p10 at days 10 to 13 at 350 to 600 against a band of about 700 to 3000), residual included, BLQ not censored, quantile bands from simulated replicates.

VPC (S1 design, 90 patients, 500 replicates, 90 % interval of each percentile; Table 1 values for the reference patient):
- p10 and p50 fall inside the bands to about day 90. At days 10 to 13: p10 1.3k to 1.8k (band 0.7k to 3k), p50 10.6k to 12.6k (band 9k to 20k).
- p90 at the peak is about 2x low: 44k to 56k against 75k to 170k. The covariate variant (Table 1 log-Cmax effects, covariates drawn from Table 2 frequencies) gives 68k to 90k, inside the band. This agrees with the comedication sensitivity here, where the right-panel bands pass.
- After day 150 all percentiles fall below the bands in every variant. In the pooled 185 to 380 d bin: p10 21 to 7.6 against 75 to 110, p50 122 to 55 against 200 to 420, p90 673 to 282 against 800 to 2800. The covariate model (p50 157 at day 183), omega as variance and LOQ handling do not close the gap. The published late bands would need FB x Cmax about 3x higher.

These are the same failures as here (p90 at the peak and the persistence phase). Together with the same structure and reading, they point to the paper side, not a transcription error. Grade unchanged.

## Review

Independent review, 2026-11-04: the review agreed with grade C (runnable, not matching; the result is failure under an assumed design, not definitive failure of the authors' simulation). Equations, parameters, the SD reading of omega and the s002.csv finding (S1 is simulated) were found sound. One fix round on the review findings: (1) the unstated design choices recorded as Fixes with evidence, as published redefined as the scheduled visits; (2) two misassigned points of the right p90 reference corrected by a stated rule, tolerance unchanged, old and new scores in Targets; (3) replicates 400 to 800 so the MC SE is below the 0.003 rule (max 0.0023); (4) report wording fixed: the assumed design is not published, the 'base model' wording is in Results (not in the caption), MC SE stated correctly. The grade comes from the new scores: still C.
- The mlxtran code uses `tmax` (lower case) in one branch while the parameter is `Tmax`.
- The mlxtran comment says non-recipients of a comedication have regressor value 99999; the data set S1 uses 9999. No effect, because the time is clipped at Tmax.
- Table 1 gives the Log Cmax covariate effects as negative numbers (-0.63, -0.36), which only make sense as exp(theta x), not as theta^x in Eq. 6.
- The text says Cmax is "twofold" higher with tocilizumab, while the Table 1 effect of 0.44 gives exp(0.44) = 1.55.
- The text links Figure 3a to the base model, but Table 1 holds the parameters of the full covariate model. No base-model estimates are published, so the base model of Fig 3a cannot be run as described.


## What would raise the grade

An A or B needs the authors' VPC specification (comedication regressors per patient, bins, interval level, BLQ handling) or the final model's covariate set, so that the Table 1 values can be simulated under the design that produced the bands. The most likely missing piece is the covariate/comedication structure (the persistence-phase and upper-tail gaps point at the Cmax and FB/beta side of the model, which the full covariate model and Ftoci = 1.2 change). Individual regressors and the real data are not published.

## Files

`model.jl`, `protocol.jl`, `reproduce.jl`, `reference/` (6 band CSVs, `digitise.toml`, `trace_fig3a.jl`, `check_overlay.jl`), `figures/comparison.png` (as published), `figures/comparison_fixed.png`, `figures/fig3.png` (regenerated simulated bands with one simulated cohort at nominal visit days; the paper's data and empirical lines are not shown), `figures/digitise_check*.png` (local only), `score.toml`. Notes: `notes.md`.

## Literature check (2026-11-04)

The paper's Results tie the Figure 3 VPC to the base model. The Monolix project files in the supplement enter comedication per patient as time-of-first-dose regressors, and values below the limit of quantification as censored data that are simulated in the VPC. The caption's 10th, 50th and 90th percentiles with 90 % intervals match Monolix defaults (500 simulations). The supplement writes "variance omega^2", consistent with Table 1 reporting omega as a standard deviation. The FDA clinical pharmacology review of tisagenlecleucel is not publicly available.
