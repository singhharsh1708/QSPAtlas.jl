<p align="center"><img src="docs/src/assets/logo.png" alt="QSPAtlas logo" width="120"></p>

# QSPAtlas.jl

**Can published pharmacometric and QSP models be rebuilt from their papers?** This
repository rebuilds published PK/PD and quantitative systems pharmacology models in
Julia ([ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl),
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)), compares them with
the curves in the papers' own figures, and grades each one.

## First results: 6 oncology models from CPT: Pharmacometrics & Systems Pharmacology

**1 of 6 reproduces its target figures, and only after two unstated inputs were read from
the figure itself. The other 5 run once their gaps are filled from documented sources, but
do not match.** No model needed a parameter that had no source at all, so none is graded D.
Five are graded on typical-value curves and one on a visual predictive check.

| model | what it is | grade | curves passing (as printed / with fixes) | main reason |
|---|---|---|---|---|
| [Ouerdani 2015](models/ouerdani_2015/report.md) | pazopanib, tumour growth + angiogenesis (mouse, human) | **B†** | 0 of 4 / 4 of 4 | initial tumour size and mouse dose group not stated; taken from the figure |
| [Hutchinson 2016](models/hutchinson_2016/report.md) | antiangiogenic transient vessel normalisation | **C** (Fig 2: B, Fig 6: C) | 0 of 30 / 22 of 30 | K = T0/V0 printed inverted; Fig 6 drug scale ~100x off Eq. 9 |
| [Lindauer 2017](models/lindauer_2017/report.md) | pembrolizumab translational PK/RO/tumour growth | **C** | 0 of 5 / 1 of 5 | human blood volume and base scenario not stated; 6 of 64 tornado bars off |
| [Fornari 2019](models/fornari_2019/report.md) | carboplatin bone-marrow toxicity QSP (rat, human) | **C** | did not run / 1 of 9 | equations not closed as printed; dose units inconsistent with its own PK figure |
| [Stroh 2019](models/stroh_2019/report.md) | masked tumour-activated antibody (Probody) QSP | **C** | 14 of 35 / 18 of 35 | human cleavage rate given only as an upper bound; "Total Drug" readout undefined |

**Population figures.** Papers whose figures show population spread (visual predictive
checks) are graded on the authors' *simulated* percentile intervals, using the population
simulation layer below; observed-data lines are a diagnostic only.

| model | what it is | grade | bands passing (as printed / with fixes) | main reason |
|---|---|---|---|---|
| [Stein 2019](models/stein_2019/report.md) | tisagenlecleucel CAR-T cellular kinetics, VPC | **C** | 2 of 6 / 2 of 6 | early expansion reproduces; upper tail and persistence phase fall below the authors' bands; VPC design (bins, BLQ, comedication) not stated |

Grades: **A** matches as published; **B** matches after documented fixes, each with
evidence; **B†** a B that needed an input read from the figure being checked; **C** runs
but does not match; **D** cannot be run from published information. How the figures were
digitised, how curves are scored, and what may fill a gap are on the method page.

## See it

**A model that reproduces.** Ouerdani 2015, Figure 4: our simulation (blue) against points
digitised from the paper's own curves (black). All four curves pass.

![Ouerdani 2015 reproduction](models/ouerdani_2015/figures/comparison.png)

```julia
using QSPAtlas
m = Module(); Base.include(m, joinpath(models_dir(), "ouerdani_2015", "reproduce.jl"))
rep = m.reproduce()        # simulates, scores every curve, writes score.toml and the figure above
grade_label(rep)           # "B†"
```

**A population figure that does not.** Stein 2019, Figure 3a: percentile bands from 800
simulated replicates of the 90-patient study. The early expansion matches the authors'
bands; the upper tail and the persistence phase fall below them.

![Stein 2019 VPC](models/stein_2019/figures/fig3.png)

**Population simulation.** 500 virtual mice from the Ouerdani model with the paper's
between-subject variability and residual error.

![500 virtual mice](examples/population_ouerdani.png)

```julia
mice = virtual_population(500; regimen = Regimen(Dose[]), obstimes = collect(0.0:1.0:75.0),
                          covariate_sampler = rng -> (AUC = 1140.8,))
res = simulate(model, mice; seed = 2015)                  # threaded, reproducible
p = percentiles(res, :P; q = (0.05, 0.5, 0.95))
```

Full script: [`examples/population_ouerdani.jl`](examples/population_ouerdani.jl).

**Performance.** Measured on standard GitHub Actions runners (4 threads, Julia 1.10) by the
Benchmarks workflow ([`benchmarks/population/run.jl`](benchmarks/population/run.jl)); median
of 3 runs after compilation. Runner hardware varies, so three runs are shown as a range.

| model | 100 subjects | 1,000 subjects | 10,000 subjects |
|---|---|---|---|
| Ouerdani 2015 tumour model (2 states) | 0.10 to 0.15 s | 0.9 to 1.4 s | 11 to 17 s |
| Two-compartment PK, 28 infusions, weight covariate | 0.11 to 0.17 s | 0.9 to 2.0 s | 11 to 18 s |

About 1 to 2 ms per simulated subject, linear in cohort size. A visual predictive check with
200 replicates of a 100-subject cohort (20,000 simulated subjects, 3 quantiles, 14 bins)
takes 25 to 39 s.

**Solver choice.** Work-precision benchmarks on the exact problems each model uses
([full results](benchmarks/solvers/results/summary.md)). Explicit Runge-Kutta methods are
fastest on the non-stiff models (Vern7, Tsit5 or BS3, 0.06 to 2.6 ms per solve to a relative
error of 1e-6); on the stiff antibody models (Lindauer 2017, Stroh 2019) they fail, and the
BDF methods QNDF and FBDF are fastest (2 to 4 ms), with Rodas5P best at tight tolerances.

![Work-precision diagram, Stroh 2019 human projection](benchmarks/solvers/results/stroh_2019__human_4_5mgkg_1000x.png)

## How it works

```
models/<author_year>/
├── model.jl       # equations and parameter values as printed, each value with its source
├── protocol.jl    # dosing, targets, tolerances, documented fixes, benchmark problems
├── reproduce.jl   # simulate, compare with digitised figure points, write score.toml + figures
├── reference/     # digitised points (one CSV per curve) and the settings that regenerate them
├── figures/       # comparison plots, regenerated paper figures, digitising overlays
├── score.toml     # grade and per-curve errors
├── report.md      # targets, results, fixes with evidence, problems found in the paper, review
└── notes.md       # research notes: transcribed equations, parameter sources, target-freeze record
```

- Targets are the authors' own simulation lines, chosen and frozen before any comparison.
- A curve passes when the RMS distance from the digitised points to the simulated curve is
  at most 1.5 % of the figure's axis box and the 95th percentile at most 3 %.
- Re-estimating parameters never counts as a fix; recalibration is reported separately.
- Each report has a review section recording what was checked, what was disputed and
  how it was settled.

```julia
using QSPAtlas
m = Module(); Base.include(m, joinpath(models_dir(), "ouerdani_2015", "reproduce.jl"))
rep = m.reproduce()   # simulates, compares, writes score.toml and figures
```

`test/runtests.jl` reruns every model and fails if any grade or curve result changes.

## Population simulation

An open simulation layer in the spirit of NONMEM/Pumas, built on ModelingToolkit:

- dosing regimens: bolus, repeated doses (`ii`, `addl`), infusions (overlapping),
  bioavailability, lag time, steady state (`ss`), resets;
- `Subject`s with baseline and time-varying covariates, `PopulationModel` with
  between-subject variability and additive, proportional, combined or exponential error;
- `simulate` for one subject or a cohort, reproducible per (replicate, subject) and
  independent of thread count;
- `percentiles`, `vpc` with Monte Carlo standard errors, and plots.

```julia
# `model` is a PopulationModel (see examples/population_ouerdani.jl for a full one)
regimen = Regimen(100.0; cmt = :central, ii = 24.0, addl = 6)      # 100 every 24 h, 7 doses
pop = virtual_population(500; regimen, obstimes = collect(0:1.0:168))
res = simulate(model, pop; seed = 1)
v = vpc(model, pop; endpoint = :conc, q = (0.05, 0.5, 0.95), bins = collect(0:12.0:168), replicates = 200)
```

See `examples/population_ouerdani.jl` for 500 virtual mice from the Ouerdani model.

## Also in the repository

- `benchmarks/solvers/`: work-precision benchmarks of OrdinaryDiffEq solvers on the exact
  problems used for each figure.
- `benchmarks/calibration/`: Turing.jl (Bayesian) and PEtab.jl (maximum likelihood)
  recalibration of the Hutchinson model against the paper's own mean tumour data.
- `tools/digitise/`: the pixel digitiser used for every reference curve.

## Licences and sources

Code is MIT. Papers are open access on PubMed Central; no figure images or supplementary
files are redistributed here. Reference CSVs hold coordinates digitised from the published
figures, each with its source figure and method in the header.
