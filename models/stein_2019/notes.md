# stein_2019 working notes (2026-11-04)

Paper: Stein et al. CPT:PSP 2019;8:285-295, PMC6539725. Target: Figure 3a (VPC, two panels: 0-380 d and 0-30 d).

## Findings before simulating

- Caption: blue lines = empirical 10/50/90 percentiles of the DATA (not graded); blue shaded = prediction interval of the simulated 10th and 90th percentiles; pink shaded = prediction interval of the simulated 50th. Monolix VPC. Six band targets (three per panel), kind band.
- `PSP4-8-285-s002.csv` is the paper's "Data set S1: simulated dataset based on the model" (not patient data). It has the same design as the real study (90 patients, 836 rows, tocilizumab 35.6 % / median 5.75 d, steroids 25.6 % / median 7.45 d, identical to Table 2), but its values (median about 1.2e5 at days 9-12, maximum 2.4e8) are far above the data in Figure 3a (maximum about 3e5). Local-only diagnostic.
- The supplement equation dM/dt = alpha E - beta M is a typo (should be k E - beta M, with k = FB (alpha - beta), from C2 = Cmax k/(alpha - beta) = Cmax FB).
- Eq. 6 z = theta^x is inconsistent with the Table 1 label "log Cmax covariate effect" (additive in log). Full model vs base model: Methods says the full covariate model was selected as final, Results says the final model is the base model.
- Table 1 typical values are self-consistent with the text: doubling time 0.78 d (ln 2 / (ln 3900/9.3)), half-lives 4.3 and 217-220 d.
- Axis calibration: y tick labels (10,000, 1,000, 100) used; the "100,000" label is displaced; about 3 px uncertainty (0.06 decade). See `reference/digitise.toml`.

## Frozen at 2026-11-04 14:50 IST (before any simulation was compared with the bands)

Targets `reference/fig3a_{right,left}_{p10,p50,p90}.csv` (kind band, default tolerance RMS 0.015 / P95 0.03), `reference/digitise.toml`.

Scored parts: right p90 4.5-29.5 d, p50 3.4-4.4 and 6.4-29.5, p10 6.4-29.5; left all 30-185 d.

Design: base model, no covariates / comedication; Table 1 omega as SD, diagonal; residual log-normal 0.56; 90 patients, all observed at every grid time (0.25 d on (0, 30], 1 d on [31, 190]); quantile type 7; residual included; BLQ not censored (simulated values kept); 200 replicates in blocks of 50 (seed 20260411 + block); statistic = median over replicates of the replicate percentile. Left panel compared at nodes 28 / 68 / 190 d joined by chords (how Monolix draws bands), right panel dense.

Decision rule on Monte Carlo error (set before looking at any comparison): if the maximum MC standard error exceeds 0.003 box units, increase the replicates to 400.

Sensitivities pre-registered (not fixes): omega as variance; comedication model with assumed timings; BLQ censored to LLOQ; left panel dense instead of chord.

14:57 MC pilot with 200 replicates: maximum MC SE 0.0047 box units, above the 0.003 rule, so replicates = 400 (before any band comparison).

## First comparison (15:02-15:13, dense / chord design, 400 replicates, before review)

rms: left p10 0.107, p50 0.058, p90 0.033; right p10 0, p50 0.005, p90 0.052 (two spike reference points). Grade C. All sensitivities fail. Review (2026-11-04) agrees with C.

## Review round (15:15-15:45)

- The FIXES now record the unstated design choices (7, each with evidence). As published = scheduled visits (1, 2, 7, 11, 14, 21, 28, 91.3, 182.6 d) joined by chords; with fixes = dense grid on the right panel and left nodes 28 / 68 / 190 (:from_target).
- Reference: two misassigned points of right p90 (x = 8.8772 and 9.9825) dropped by an outlier rule (documented in `reference/digitise.toml` and `reference/trace_fig3a.jl`). Tolerance unchanged.
- Replicates 800 (400 gave a maximum MC SE of 0.0034, above 0.003; 800 gives 0.0023). `reproduce()` takes about 186 s in a fresh session.
- New scores: as published, rms left p10 0.108, p50 0.059, p90 0.047, right p10 0, p50 0.010, p90 0.026; with fixes, left 0.108 / 0.060 / 0.035, right 0 / 0.006 / 0.024. Grade C. The comedication sensitivity passes 4 of 6 bands (left p10 and p50 fail) and is a lead only.

Self-check (2026-11-04): an independent rebuild from the paper alone was compared with these results; see the Self-check section in `report.md`.

## Provenance ledger (2026-11-05)

The seven shared design fills are listed in both runs (two new entries: `omega_as_sd`, `omega_diagonal`, applied before without an entry); the as-published run also lists `scheduled_visits` and `days_per_month`, the fixed run `dense_grid` and `left_bin_nodes`. Without them the band reading throws NotRunnable. No number changed; grade C, strict D.

## VPC design reconstruction (2026-11-05)

Aim: an assumed reconstruction of the VPC, informed by the Monolix documentation (not a verified copy of
what the authors ran), from sources that are not Figure 3a, and grade it. It is not blind to the target: it
comes after the first comparison (2026-11-04) and after the previous design's scores were known, and the
comedication regressors (the source of its gain, 4 of 6 instead of 2 of 6) had been run on 2026-11-04 as the
comedication sensitivity (4 of 6 bands) and shown to improve agreement before they were adopted. The other
choices were fixed from the sources before the reconstruction was compared with the bands. The 4 of 6 is
conditional on the assumptions below.

Documented vs assumed (version matters; the paper does not state its Monolix version, probably 2016R1 to
2019R1, whose pages could not be retrieved):
- Monolix 2024R1 docs: 500 simulations, 90 % PI, all simulated values incl. BLQ in the percentiles
  (:reference); R5 percentiles (:assumption for the paper's era); binning equal width = default,
  least-squares = option.
- Monolix 4 (Lavielle & Bleakley 2011): least-squares automatic binning, 90 % intervals. Default of the
  paper's version not documented, so least-squares binning is an :assumption.
- No document: bin node at the bin mean time; comedication regressor spreads (log-normal, independent,
  sd of log 0.7 and 1.2, clamps) from Table 2 marginals, no independent evidence. :assumption.

Sources read:
- Paper: Monolix used, version not stated; BLQ fixed-point censored (CENSORING column); Results tie Fig 3
  to the base model and say the final model equals the base model; Results "longest follow-up ... 1 year";
  Methods visit schedule (days 1, 2, 7, 11, 14, 21, 28, 3-monthly months 3-12).
- Supplement S1 (s001.pdf lines 100-152): the mlxtran has only the model block (no VPC or task settings);
  parameter = {foldx, Tmax, Cmax, Ftoci, Fster, alpha, FB, beta}, regressor = {TOCI1T, STERSTT}. So the
  model that was fitted (base = final) carries the comedication effect through per-patient regressors.
- Data set S1 (s002.csv): 90 ids, 836 rows, times 0.88-379.8 d, TOCI1T for 32 patients, STERSTT for 23
  (9999 = none). Same design as the study (Table 2). Not redistributed (DATA_LICENSES.md), so CI cannot read it.
- Monolix 2024R1 VPC page: 10/50/90 percentiles, 90 % PI, 500 simulations, percentiles "R5"; binning
  equal width (default), equal size or least-squares, automatic bin count (rule not documented).
- Monolix 2024R1 "Correcting the VPC for missing observations": with censored data used and simulated, all
  simulated values (also below LOQ) enter the simulated percentiles; observed BLQ are drawn from the
  conditional distribution.
- Lavielle & Bleakley 2011 (JPKPD 38:861, HAL hal-00654842): Monolix 4 automatic binning, J = sum n_k
  var_k^beta, beta 0.2, K by log J + lambda beta K, lambda 0.3; 90 % intervals.
- Older version-specific Monolix pages (2016R1-2019R1) could not be retrieved (the old lixoft.com pages
  redirect to 2024R1).

Graded "with fixes" run (CI): Monolix settings above on the published design: 90 patients at the nominal
visits to 12 months, regressors drawn from Table 2 (the 2026-11-04 comedication sampler), base model with
Ftoci 1.2 / Fster 1.0, residual error, least-squares bins (one per visit on nominal times), type 5
percentiles, all simulated values kept, 90 % interval. 500 replicates unless the MC SE rule needs more.
Local check (not graded): the same on the S1 design (own times and regressors; least-squares gives 8 bins,
nodes 2.9, 6.5, 12.0, 20.3, 28.6, 83.5, 175, 311 d).

Results (2026-11-05, logs/2026-11-05-stein-vpc-*.log): graded reconstruction 4 of 6 (rms left p10 0.096,
p50 0.038 fail; left p90 0.011, right p10 0.000, p50 0.005, p90 0.007 pass), max MC SE 0.0027 box with 500
replicates, so 500 kept. As published unchanged (2 of 6). Previous fixed design 2 of 6 (reported only).
Grade C, strict D. S1 design (local, 2000 replicates, max MC SE 0.0025): 2 of 6 (right p10, p90); right p50
0.0155 and left p90 0.0153 just above 0.015. S1 without regressors (500 replicates): 1 of 6.
Fresh-module check: grade and every curve equal score.toml; 310 s with compilation (over the 4 min budget;
cost is per simulated subject, solver tolerance 1e-8 vs 1e-10 makes little difference). src/ not edited;
type 5 percentiles and least-squares binning are local to protocol.jl.

## Wording and metadata fix (2026-11-05, independent review finding 7)

The reconstruction was described as what Monolix builds by default. Reworded as an assumed reconstruction
informed by the Monolix documentation: the binning method (least-squares, a Monolix 4 / 2024R1 option, equal
width is the 2024R1 default) and the R5 convention (2024R1 only) are not documented for the paper's version,
and the regressor spreads have no independent evidence, so those fixes now carry basis :assumption (stage
after_comparison). Only text and Fix metadata changed; no number, design, target or tolerance. The settings
were chosen after the 2026-11-04 comedication sensitivity had shown that the regressors improve agreement,
so the 4 of 6 bands are conditional on the assumptions.
