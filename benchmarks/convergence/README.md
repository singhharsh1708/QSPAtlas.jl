# Convergence of the atlas scores

Question: does any score (rms, p95, pass/fail of a curve) depend on solver settings? The solver benchmarks in
`benchmarks/solvers` measure state error and speed, not scores. This study reruns the scoring itself.

Run it from the repository root (the five deterministic models take about 20 minutes on a loaded 4-thread machine; Stein with 800 replicates about 45 minutes, `CONV_STEIN_REPS=100` about 5):

```
julia --startup-file=no --threads=4 --project=. benchmarks/convergence/run.jl            # everything
julia --startup-file=no --threads=4 --project=. benchmarks/convergence/run.jl stroh_2019   # one model
julia --startup-file=no --threads=4 --project=. benchmarks/convergence/run.jl report       # rebuild results.csv and tables.md
julia --startup-file=no --threads=4 --project=. benchmarks/convergence/run.jl nostein      # the five deterministic models, report, check (the CI job)
julia --startup-file=no --threads=4 --project=. benchmarks/convergence/run.jl check        # only the check, on the committed results
```

After a run the script checks its results and exits non-zero if a curve flips pass/fail in any variant, a solve ends without a successful retcode (except the three documented Lindauer published bars, listed in `ALLOWED_FAILED_SOLVES` in `run.jl`), or an rms moves by more than 10 % of the curve's rms tolerance (Stein's `reps x2` and `seed+1000` variants, which size the Monte Carlo noise, are exempt from the last test).

Status note: the study was first run on the code before the changes listed under "Recommended changes"; the tables, `results.csv` and `solves.csv` below are from that run (Stroh's published run with 1500 saved points). The first run of the `convergence` job in `.github/workflows/Benchmarks.yml` regenerates them. With the Stroh grid raised to 6000 points the dense x4 variant moves rms by about 1.5e-7 instead of 0.0025, so the check passes.

The study itself edits nothing in `src/`, `models/`, `docs/` or `test/`. The model modules are loaded as they are and their
problem builders are reused. Every solve call that hard-codes the solver is copied into `run.jl` with the solver
as a parameter (each copy names its source function).

## Variants

| variant | what changes |
|---|---|
| own | the settings in the model's `protocol.jl` (checked against `score.toml`, see "Fidelity") |
| tol/100 | abstol and reltol 100 times tighter |
| altA, altB | two other solvers: Tsit5 models get Rodas5P (stiff) and Vern9; Rodas5P models get FBDF and Rosenbrock23 (Lindauer, Stroh, stiff) or Vern9 and FBDF (Fornari); Stein (Vern9) gets Rodas5P and Tsit5. Same tolerances as own. |
| dense x4 | 4 times as many saved points (Ouerdani 1000, Hutchinson 0.05 d grid, Fornari 0.05 d, Stroh 1500, Lindauer 366) or observation times (Stein) |
| Stein only: reps x2, seed+1000 | twice the replicates, and another seed, to size the Monte Carlo noise the solver changes are compared with |

Both runs that produce a `score.toml` entry are covered ("published" and "fixed"). Fornari's published run is
`NotRunnable` by design (the printed equations are not closed), so it has no curves.

Per curve and variant `results.csv` holds rms, p95, max, pass/fail (the same `compare_curve` and tolerances as
`run_curves`), the own pass/fail, the own margin (max of rms/tol and p95/tol), and the largest shift of the
simulated curve against the own curve in axis-box units (the unit of the score). `solves.csv` has one row per
solve: retcode, finiteness, number of saved points, seconds.

## Result in one paragraph

No pass/fail result changes under any variant in any model, and every curve that sits close to a threshold
(own margin between 0.67 and 1.5) keeps its verdict. Solver error is far below the score tolerance (rms 0.015)
everywhere: tolerance, solver and dense-output changes move rms and p95 by at most 2e-6 for Ouerdani, Hutchinson,
Fornari, Lindauer and Stein (solver and grid changes), and for Stroh's solver changes by at most 5e-7 (rms) and 1.3e-6 (p95). Two findings are not about
tolerances: (1) Stroh's published run depends on the output grid for two curves (rms moves by 0.0025, both curves
fail by a wide margin), (2) Lindauer's published run has three bars whose solve fails in every variant
(`Unstable` / `MaxIters`), scored as the maximum error. Retcode checks were missing in Ouerdani, Hutchinson and
Fornari (added since, see the last section); nothing failed there. The Fornari human Fig 4b reference solve succeeds today (details below).

## Stein (Monte Carlo)

Same seed in every variant: the eta draws depend only on (seed, replicate, patient), and the residual draws are
re-applied in `run.jl` from the production streams at the 285 baseline times, so that a variant with 4 times as many
observation times sees identical residuals. The own variant equals the production `simulate_design` summary exactly
(max abs difference of the median 0). **Because of the time budget the Stein variants were run with 100 replicates
(`CONV_STEIN_REPS=100`), not the protocol's 800.** A first 800-replicate run (about 45 minutes under load, all variants)
was stopped before it finished. The comparison of variants is unaffected (same seed, same draws, only the solver or
grid changes), but the absolute rms/p95 here are 100-replicate values and are not those of `score.toml`
(the fidelity table below shows that difference for Stein). Rerun with `run.jl stein_2019` for the full 800.

- Solver and grid changes: the compared medians move by at most 4e-10 of the axis box (tol/100, Rodas5P, Tsit5,
  4 times as many observation times), rms and p95 by at most 6e-11 and 1.3e-10, and the shifts are below 2e-7 Monte Carlo
  standard errors (`results/stein_2019_mc_shifts.csv`). No pass/fail changes. The solver is irrelevant for Stein.
- Deterministic part: ipred of 300 sampled individuals (and the typical individual) against the closed form of Eq. 1
  (`results/stein_2019_deterministic.csv`): maximum relative error 6.6e-10 (own, Vern9 at 1e-10), 7.5e-12 (tol/100),
  1.0e-9 (Rodas5P), 8.3e-9 (Tsit5), 1.0e-9 (dense). The right-hand side has kinks at Tmax, T1, T2 (`ifelse` on time) that
  the solver is not told about; with these tolerances the error control resolves them.
- Monte Carlo noise, for scale: another seed (seed + 1000) moves rms by up to 0.0076 and p95 by 0.0105, and 200 instead
  of 100 replicates by up to 0.0035 and 0.0064 (the 100-replicate values; the noise shrinks with the square root of the
  replicates, so it is about 3 times smaller at 800). Every Stein pass/fail verdict is the same under all of these. No
  curve is within 1.5 times of a threshold (right_p50 is the closest to passing: rms 0.009, p95 0.015 against 0.015 and 0.03;
  right_p90 fails with rms 1.7 times the tolerance), so Monte Carlo noise does not decide any verdict either.
- 9000 trajectories per variant (18000 for reps x2) all end with `Success` and finite output.


## Tables (generated by `report`, see tables.md)

### Per model and variant: largest score change against the own settings

`d rms` and `d p95` are the largest absolute change of rms / p95 over the curves of the run (axis-box units; the pass rule is rms <= 0.015 and p95 <= 0.03 unless a curve has its own tolerance). `shift` is the largest change of any simulated point in axis-box units. `flips` counts curves whose pass/fail differs from own.

| model | run | variant | curves | max d rms | max d p95 | max shift (box) | flips | failed solves | non-finite |
|---|---|---|---|---|---|---|---|---|---|
| fornari_2019 | fixed | tol/100 | 9 | 4.9e-11 | 1.4e-10 | 2.17e-10 | 0 | 0 | 0 |
| fornari_2019 | fixed | altA Vern9 | 9 | 5e-11 | 1.4e-10 | 2.19e-10 | 0 | 0 | 0 |
| fornari_2019 | fixed | altB FBDF | 9 | 1.08e-09 | 1.92e-09 | 8.21e-09 | 0 | 0 | 0 |
| fornari_2019 | fixed | dense x4 | 9 | 4.96e-07 | 1e-06 | 0 | 0 | 0 | 0 |
| hutchinson_2016 | published | tol/100 | 30 | 2.53e-08 | 4.07e-08 | 5.04e-08 | 0 | 0 | 0 |
| hutchinson_2016 | published | altA Rodas5P | 30 | 2.78e-08 | 4.15e-08 | 5.14e-08 | 0 | 0 | 0 |
| hutchinson_2016 | published | altB Vern9 | 30 | 2.41e-08 | 3.83e-08 | 4.75e-08 | 0 | 0 | 0 |
| hutchinson_2016 | published | dense x4 | 30 | 2.32e-07 | 3.18e-07 | 0 | 0 | 0 | 0 |
| hutchinson_2016 | fixed | tol/100 | 30 | 4.25e-08 | 5.79e-08 | 6.56e-08 | 0 | 0 | 0 |
| hutchinson_2016 | fixed | altA Rodas5P | 30 | 4.25e-08 | 5.86e-08 | 7.48e-08 | 0 | 0 | 0 |
| hutchinson_2016 | fixed | altB Vern9 | 30 | 4.05e-08 | 5.52e-08 | 6.16e-08 | 0 | 0 | 0 |
| hutchinson_2016 | fixed | dense x4 | 30 | 6.68e-07 | 1.98e-06 | 0 | 0 | 0 | 0 |
| lindauer_2017 | published | tol/100 | 5 | 0 | 0 | 3.58e-09 | 0 | 3 | 0 |
| lindauer_2017 | published | altA FBDF | 5 | 0 | 0 | inf | 0 | 2 | 0 |
| lindauer_2017 | published | altB Rosenbrock23 | 5 | 1e-10 | 0 | 0.00306 | 0 | 3 | 0 |
| lindauer_2017 | published | dense x4 | 5 | 0 | 0 | 0 | 0 | 3 | 0 |
| lindauer_2017 | fixed | tol/100 | 5 | 1e-10 | 6e-10 | 1.42e-09 | 0 | 0 | 0 |
| lindauer_2017 | fixed | altA FBDF | 5 | 3.3e-09 | 1.51e-08 | 4.13e-08 | 0 | 0 | 0 |
| lindauer_2017 | fixed | altB Rosenbrock23 | 5 | 1.04e-08 | 4.71e-08 | 1.06e-07 | 0 | 0 | 0 |
| lindauer_2017 | fixed | dense x4 | 5 | 0 | 0 | 0 | 0 | 0 | 0 |
| ouerdani_2015 | published | tol/100 | 4 | 1e-11 | 1e-10 | 2.63e-10 | 0 | 0 | 0 |
| ouerdani_2015 | published | altA Rodas5P | 4 | 2e-11 | 0 | 2.42e-10 | 0 | 0 | 0 |
| ouerdani_2015 | published | altB Vern9 | 4 | 1e-11 | 1e-10 | 2.62e-10 | 0 | 0 | 0 |
| ouerdani_2015 | published | dense x4 | 4 | 1.57e-06 | 3.35e-07 | 0 | 0 | 0 | 0 |
| ouerdani_2015 | fixed | tol/100 | 4 | 4e-13 | 0 | 7.61e-11 | 0 | 0 | 0 |
| ouerdani_2015 | fixed | altA Rodas5P | 4 | 1e-11 | 1.4e-11 | 1.44e-10 | 0 | 0 | 0 |
| ouerdani_2015 | fixed | altB Vern9 | 4 | 4e-13 | 1e-12 | 7.53e-11 | 0 | 0 | 0 |
| ouerdani_2015 | fixed | dense x4 | 4 | 1.72e-07 | 1.12e-07 | 0 | 0 | 0 | 0 |
| stein_2019 | published | tol/100 | 6 | 1e-11 | 3e-11 | 2.89e-11 | 0 | 0 | 0 |
| stein_2019 | published | altA Rodas5P | 6 | 0 | 1e-11 | 2.19e-11 | 0 | 0 | 0 |
| stein_2019 | published | altB Tsit5 | 6 | 6e-11 | 1.3e-10 | 1.68e-10 | 0 | 0 | 0 |
| stein_2019 | published | dense x4 | 6 | 1e-11 | 2e-11 | 2.75e-11 | 0 | 0 | 0 |
| stein_2019 | published | reps x2 | 6 | 0.00227 | 0.00638 | 0.0103 | 0 | 0 | 0 |
| stein_2019 | published | seed+1000 | 6 | 0.00764 | 0.0105 | 0.0126 | 0 | 0 | 0 |
| stein_2019 | fixed | tol/100 | 6 | 1e-11 | 1e-11 | 5.34e-11 | 0 | 0 | 0 |
| stein_2019 | fixed | altA Rodas5P | 6 | 0 | 1e-11 | 4.14e-11 | 0 | 0 | 0 |
| stein_2019 | fixed | altB Tsit5 | 6 | 1e-11 | 1e-10 | 3.77e-10 | 0 | 0 | 0 |
| stein_2019 | fixed | dense x4 | 6 | 1e-11 | 2e-11 | 7.48e-11 | 0 | 0 | 0 |
| stein_2019 | fixed | reps x2 | 6 | 0.00345 | 0.00202 | 0.0116 | 0 | 0 | 0 |
| stein_2019 | fixed | seed+1000 | 6 | 0.00587 | 0.0042 | 0.0156 | 0 | 0 | 0 |
| stroh_2019 | published | tol/100 | 35 | 1e-10 | 5.3e-10 | 7.31e-05 | 0 | 0 | 0 |
| stroh_2019 | published | altA FBDF | 35 | 1.53e-08 | 1.59e-08 | 0.00757 | 0 | 0 | 0 |
| stroh_2019 | published | altB Rosenbrock23 | 35 | 5.06e-07 | 1.31e-06 | 0.00966 | 0 | 0 | 0 |
| stroh_2019 | published | dense x4 | 35 | 0.00254 | 0.00431 | 0 | 0 | 0 | 0 |
| stroh_2019 | fixed | tol/100 | 35 | 5.4e-11 | 5.3e-10 | 7.31e-05 | 0 | 0 | 0 |
| stroh_2019 | fixed | altA FBDF | 35 | 1.53e-08 | 1.59e-08 | 0.00757 | 0 | 0 | 0 |
| stroh_2019 | fixed | altB Rosenbrock23 | 35 | 5.06e-07 | 1.31e-06 | 0.00966 | 0 | 0 | 0 |
| stroh_2019 | fixed | dense x4 | 35 | 1.52e-07 | 2.56e-07 | 0 | 0 | 0 | 0 |

### Per curve: the three curves with the largest change of rms in each model and run

Largest change over all variants; the full per-curve, per-variant table is results.csv.

| model | run | curve | own rms | own p95 | max d rms | max d p95 | max shift (box) | variant of largest d rms |
|---|---|---|---|---|---|---|---|---|
| fornari_2019 | fixed | fig3_Plt | 0.00722 | 0.0155 | 4.96e-07 | 1e-06 | 8.21e-09 | dense x4 |
| fornari_2019 | fixed | fig4a_plt | 0.0323 | 0.054 | 3.83e-07 | 9.5e-07 | 8.53e-10 | dense x4 |
| fornari_2019 | fixed | fig4b_neut | 0.00993 | 0.0179 | 1.93e-07 | 3.1e-08 | 2.76e-09 | dense x4 |
| hutchinson_2016 | published | fig6_b_T | 0.19 | 0.263 | 2.32e-07 | 2.91e-08 | 2.93e-08 | dense x4 |
| hutchinson_2016 | published | fig2_f_V_N21 | 0.699 | 0.928 | 2.18e-07 | 3.18e-07 | 2.96e-08 | dense x4 |
| hutchinson_2016 | published | fig6_d_T | 0.231 | 0.275 | 2.06e-07 | 7.94e-08 | 2.92e-08 | dense x4 |
| hutchinson_2016 | fixed | fig2_e_T_N21 | 0.00607 | 0.0116 | 6.68e-07 | 5.94e-07 | 7.48e-08 | dense x4 |
| hutchinson_2016 | fixed | fig6_d_T | 0.244 | 0.29 | 5.87e-07 | 2.9e-09 | 3.43e-08 | dense x4 |
| hutchinson_2016 | fixed | fig2_f_V_N11 | 0.00537 | 0.00874 | 5.24e-07 | 4.09e-07 | 5.94e-08 | dense x4 |
| lindauer_2017 | published | fig4_5mgkg_x5 | 1 | 1 | 1e-10 | 0 | inf | altB Rosenbrock23 |
| lindauer_2017 | published | fig4_0p5mgkg_x0p2 | 0.939 | 1 | 0 | 0 | 4.25e-09 | tol/100 |
| lindauer_2017 | published | fig4_5mgkg_x0p2 | 0.937 | 1 | 0 | 0 | 2.03e-09 | tol/100 |
| lindauer_2017 | fixed | fig4_0p5mgkg_x5 | 0.131 | 0.359 | 1.04e-08 | 4.71e-08 | 6.28e-08 | altB Rosenbrock23 |
| lindauer_2017 | fixed | fig4_0p5mgkg_x0p2 | 0.0964 | 0.0992 | 8e-11 | 8e-11 | 4.47e-10 | altB Rosenbrock23 |
| lindauer_2017 | fixed | fig4_base_case | 0.000906 | 0.00117 | 1.74e-11 | 1e-11 | 5.51e-11 | altB Rosenbrock23 |
| ouerdani_2015 | published | fig4_right_P | 0.064 | 0.111 | 1.57e-06 | 7.5e-08 | 1.42e-10 | dense x4 |
| ouerdani_2015 | published | fig4_right_K | 0.0331 | 0.0424 | 2.4e-07 | 3.35e-07 | 4.07e-11 | dense x4 |
| ouerdani_2015 | published | fig4_left_P | 0.081 | 0.0985 | 2.11e-07 | 7.29e-08 | 2.63e-10 | dense x4 |
| ouerdani_2015 | fixed | fig4_right_K | 0.0112 | 0.016 | 1.72e-07 | 2.92e-08 | 3.96e-11 | dense x4 |
| ouerdani_2015 | fixed | fig4_right_P | 0.00882 | 0.0138 | 1.2e-07 | 8.06e-08 | 1.44e-10 | dense x4 |
| ouerdani_2015 | fixed | fig4_left_P | 0.000865 | 0.00156 | 6.71e-08 | 4.07e-08 | 7.22e-11 | dense x4 |
| stein_2019 | published | fig3a_left_p90 | 0.0441 | 0.0613 | 0.00764 | 0.0105 | 0.0119 | seed+1000 |
| stein_2019 | published | fig3a_left_p10 | 0.107 | 0.157 | 0.00298 | 0.00512 | 0.0121 | seed+1000 |
| stein_2019 | published | fig3a_right_p90 | 0.0261 | 0.0401 | 0.00192 | 0.00542 | 0.0113 | seed+1000 |
| stein_2019 | fixed | fig3a_left_p90 | 0.0322 | 0.0463 | 0.00587 | 0.0023 | 0.0129 | seed+1000 |
| stein_2019 | fixed | fig3a_left_p10 | 0.106 | 0.156 | 0.00369 | 0.00211 | 0.0061 | seed+1000 |
| stein_2019 | fixed | fig3a_left_p50 | 0.0582 | 0.0749 | 0.00204 | 0.00169 | 0.0126 | seed+1000 |
| stroh_2019 | published | fig4c_100x | 0.0573 | 0.0822 | 0.00254 | 0.00431 | 3.08e-05 | dense x4 |
| stroh_2019 | published | fig4c_10x | 0.157 | 0.178 | 0.00231 | 0.00355 | 0.00242 | dense x4 |
| stroh_2019 | published | fig4a_100x | 0.0262 | 0.0457 | 5.06e-07 | 1.31e-06 | 3.08e-05 | altB Rosenbrock23 |
| stroh_2019 | fixed | fig4a_100x | 0.0262 | 0.0457 | 5.06e-07 | 1.31e-06 | 3.08e-05 | altB Rosenbrock23 |
| stroh_2019 | fixed | fig4b_100x | 0.0307 | 0.0516 | 3.49e-07 | 1.06e-06 | 3.08e-05 | altB Rosenbrock23 |
| stroh_2019 | fixed | fig4a_10x | 0.00287 | 0.00712 | 2.01e-07 | 4.74e-08 | 0.00242 | altB Rosenbrock23 |

### Curves whose pass/fail flips in some variant, or that sit within 1.5x of a threshold

| model | run | curve | own rms | own p95 | own pass | own margin (max of rms/tol, p95/tol) | variants that flip |
|---|---|---|---|---|---|---|---|
| fornari_2019 | fixed | fig3_Neut | 0.011 | 0.0167 | true | 0.736 | none |
| ouerdani_2015 | fixed | fig4_right_K | 0.0112 | 0.016 | true | 0.745 | none |
| stroh_2019 | published | fig3a_M1S1 | 0.0192 | 0.0313 | false | 1.28 | none |
| stroh_2019 | published | fig3b_M1S2 | 0.0177 | 0.0285 | false | 1.18 | none |
| stroh_2019 | published | fig3c_M1S2 | 0.0151 | 0.0259 | false | 1.01 | none |
| stroh_2019 | published | fig3c_mAb00 | 0.00756 | 0.0235 | true | 0.784 | none |
| stroh_2019 | published | fig4b_10x | 0.0164 | 0.0212 | false | 1.1 | none |
| stroh_2019 | published | fig4d_100x | 0.0143 | 0.018 | true | 0.954 | none |
| stroh_2019 | fixed | fig3a_M1S1 | 0.0192 | 0.0313 | false | 1.28 | none |
| stroh_2019 | fixed | fig3b_M1S2 | 0.0177 | 0.0285 | false | 1.18 | none |
| stroh_2019 | fixed | fig3c_M1S2 | 0.0151 | 0.0259 | false | 1.01 | none |
| stroh_2019 | fixed | fig3c_mAb00 | 0.00756 | 0.0235 | true | 0.784 | none |
| stroh_2019 | fixed | fig4d_100x | 0.0143 | 0.018 | true | 0.954 | none |

### Solves that did not finish with a successful retcode or gave non-finite output

| model | run | solve | variant | retcode | finite |
|---|---|---|---|---|---|
| lindauer_2017 | published | 0.5:gamma:5.0 | own | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | own | Unstable | true |
| lindauer_2017 | published | 5.0:SLtg:5.0 | own | MaxIters | true |
| lindauer_2017 | published | 0.5:gamma:5.0 | tol/100 | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | tol/100 | Unstable | true |
| lindauer_2017 | published | 5.0:SLtg:5.0 | tol/100 | MaxIters | true |
| lindauer_2017 | published | 0.5:gamma:5.0 | altA FBDF | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | altA FBDF | Unstable | true |
| lindauer_2017 | published | 0.5:gamma:5.0 | altB Rosenbrock23 | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | altB Rosenbrock23 | Unstable | true |
| lindauer_2017 | published | 5.0:SLtg:5.0 | altB Rosenbrock23 | MaxIters | true |
| lindauer_2017 | published | 0.5:gamma:5.0 | dense x4 | Unstable | true |
| lindauer_2017 | published | 5.0:gamma:5.0 | dense x4 | Unstable | true |
| lindauer_2017 | published | 5.0:SLtg:5.0 | dense x4 | MaxIters | true |

### Reproduction of score.toml by the own variant

| model | curves compared | max relative rms difference | pass/fail differences |
|---|---|---|---|
| fornari_2019 | 9 | 0.000425 | 0 |
| hutchinson_2016 | 60 | 0.000372 | 0 |
| lindauer_2017 | 10 | 0.000166 | 0 |
| ouerdani_2015 | 8 | 0.000357 | 0 |
| stein_2019 | 12 | 1 | 0 |
| stroh_2019 | 70 | 0.000394 | 0 |


## Fidelity of the re-implementation

The "own" variant reproduces `score.toml` for every curve of the five deterministic models (table above): the rms differs by
at most 4.3e-4 relative, which is the 4 significant digits `score.toml` stores, and no pass/fail differs. For Stein the own
summary computed here equals the production `simulate_design` result exactly at the same replicate count (max abs
difference of the median 0, first row of `results/stein_2019_mc_shifts.csv`); the table row for Stein compares 100 replicates
with `score.toml` (800 replicates), so its rms difference is Monte Carlo noise, with no pass/fail difference.

## Retcodes and finite output

Every solve is in `solves.csv`. Outside Lindauer's published run and the Fornari FBDF probes below, all solves in
all variants end with a successful retcode and finite output. Lindauer's published run:

| bar set | solve | retcode (own, Rodas5P) | cause |
|---|---|---|---|
| fig4_5mgkg_x5, fig4_0p5mgkg_x5 | gamma x5 | Unstable | tumour volume falls from 41500 to about 1e-11 mm^3 within the first 0.5 h (RO_tumour 73 to 100 %, drug effect SL RO^(5 gamma)); tumour-interstitial volume Vis is proportional to TV, so every concentration (1/Vis) diverges |
| fig4_5mgkg_x5 | SLtg x5 | MaxIters (t = 1912 of 4368 h) | same mechanism, reached later; Rosenbrock23 also `MaxIters`; FBDF at the own tolerances does finish and returns a value from a state with TV near 1e-12 (the shift against the NaN own curve shows as inf in the table), but at abstol 1e-10 it stops `Unstable` at t = 2329 h |

These are not tolerance problems: tol/100, Rosenbrock23 and dense output fail the same solves, FBDF the two gamma solves, and
the same model with the other parameters (base case, fixed run) is stable to 1e-8 in log10 units. The
protocol maps a failed solve to NaN (scored as the maximum error 1.0), which is the safe direction, but for these
bars "NaN" really means "tumour eradicated, log ratio below -11", and the paper's bar may extend that far. The
base case of the published run (rms 1.0) is not a solver failure: the solve succeeds and log10(TV26w/TV0) = -6.0
(V_blood = 1400 uL), as report.md says.

## Fornari human Fig 4b reference solve

`benchmarks/solvers/results/summary.md` lists "reference solve failed" for `fornari_2019 human, two 21-day cycles,
Table 2 slopes (Fig 4b)`. `run.jl fornari_diag` repeats what `benchmarks/solvers/run.jl` `try_reference` does on the
problem `benchmark_setups()` returns (Vern9, maxiters 1e5, then Rodas5P, maxiters 1e7, both at abstol = reltol =
1e-12, with the callback passed a second time) and 214 variations (`fornari_human_reference.csv`).

- Result: both reference attempts succeed on the current code (Vern9 793 steps, Rodas5P 3647 steps, same
  final state), as do Tsit5, Vern9, Rodas5P and FBDF at 1e-6 to 1e-10, and Tsit5, Vern9, Rodas5P at 1e-12 and 1e-13;
  the Fig 4a (rat slopes) and Eq 10 slope problems, the single-dose and bolus variants, and all 81 corner and
  midpoint combinations of the printed human baseline ranges (Neut0 2-7, Mono0 0.2-10, Ret0 40-115, RBC0 4100-5900) also
  succeed. Passing the callback twice (the problem already carries it, `benchmark_setups` passes it again) is
  harmless: the saved points differ by a handful of duplicates, the solution does not change. The maxiters cap of
  1e5 is not reached (Tsit5 at 1e-12 needs 3748 steps). The problem is not stiff in any way that matters: the
  explicit solvers need no more steps than Rodas5P at the same tolerance.
- The only failures found: FBDF at abstol = reltol <= 1e-12 returns `Unstable` at t = 504.0000001 h, the start of
  the second infusion (`rate_callback` sets `rin` from 0 to 4040 umol/h; 115 255 steps spent in front of it).
  Flagging the event with `derivative_discontinuity!` does not help. This is the discontinuity (a rate jump) that a
  multistep method cannot cross at that tolerance; one-step methods (Vern9, Rodas5P, Tsit5) cross it without
  trouble. The benchmark's reference never uses FBDF, and its tolerance sweep stops at 1e-11, so this does not
  explain the committed summary.
- Conclusion: the failure in `summary.md` cannot be reproduced with the current model code, solvers or
  benchmark script. `summary.md` (2026-11-04 19:14) is older than the last edits of `models/fornari_2019/protocol.jl`
  and `model.jl` (2026-11-05 11:03 and 11:19), so it may describe an earlier version of the problem; I could not
  recover that version. What is certain is that `try_reference` hides the reason: its `catch` returns `nothing`
  and a non-successful retcode is not recorded, so the summary cannot say whether it was an exception, `Unstable`
  or `MaxIters`. Rerun the benchmark (recommendation 7) to refresh the summary.

## Recommended changes and what was done with them

Applied (the score of every curve is unchanged except Stroh's two published curves, see 3):

1. Retcode checks: applied. Ouerdani `simulate_curve`, Hutchinson `simulate_scenario` and Fornari `solved`/`series` check
   `successful_retcode(sol)` before the solution is used and return NaN values on failure (scored as a gap, as in Lindauer).
   Stein `replicate_percentiles` throws an error with the count of subjects whose retcode is not `:Success` after each block
   of replicates, because a failed subject in a grade-bearing VPC must not pass silently. Stroh already errors and is unchanged.
   No solve fails today, so no score changes.
3. Stroh published run: applied. `make_problem` takes `npts = 6000` instead of 1500 (the copy in `run.jl` follows). The fixed
   run is unchanged to 1e-6; `fig4c_10x` rms 0.1573 to 0.155 and `fig4c_100x` rms 0.0573 to 0.0548 (about 0.0025 each), both
   still fail. The sampling-dependent diagnostics `rev_n`, `rev_rms`, `rev_p95` and `sim_coverage` in `score.toml` move with the grid
   (fifth digit, not graded).
6. Fornari `benchmark_setups`: applied. The callback is no longer passed in `BenchmarkSetup`, the problems already carry it.
7. `benchmarks/solvers/run.jl` `try_reference`: applied. Each failed reference attempt records the solver, `maxiters` and the
   exception message or the retcode, and `summary.md` prints them in the row instead of a bare "reference solve failed". The committed
   `summary.md` still has the old text until the Benchmarks workflow is rerun.
8. Workflow: applied. `.github/workflows/Benchmarks.yml` has a `convergence` job (`workflow_dispatch`, Julia 1.10, root project) that runs
   `run.jl nostein` and uploads the results as an artifact. `run.jl` exits non-zero on a flip, a non-successful retcode (except the three
   documented Lindauer published bars) or an rms move above 10 % of the tolerance.

Not applied, documented only:

2. Lindauer: no terminal callback or floor. The three bars (gamma x5 at both doses, SLtg x5 at 5 mg/kg) stay NaN, scored as complete
   misses; a floor would not change their pass/fail (the paper's bars are far above the simulated -11 log10 units). One sentence is
   added to `models/lindauer_2017/report.md` (Issues). The retcode check stays.
4. Solver and tolerance choice: no change; the current choices are validated by this study, loosening the tolerances is optional and
   was not requested.
5. Fornari and FBDF below abstol 1e-11: nothing to change, the benchmark sweep already stops at 1e-11 and no model code uses FBDF.

## Check rule (2026-11-05)

The check fails on any pass/fail flip, any unexpected failed solve, and any rms move larger than
10 % of the tolerance for curves whose rms is within 3 times the tolerance. Clearly failing curves are
exempt from the rms-move rule because they cannot flip and their rms can depend on the output grid: the
first run of the check on GitHub flagged Stroh published `fig4c_10x` (rms about 0.155, ten times the
tolerance) moving by 0.0145 between 6,000 and 24,000 output points, where the simulation dips to zero on
a log axis and the start of the gap moves with the sampling.

