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
