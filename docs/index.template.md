# QSP Reproducibility Atlas

> **Work in progress.** These readings of the papers have not yet been checked with the
> papers' authors. Corrections are welcome: please open an issue on [GitHub](https://github.com/singhharsh1708/QSPAtlas.jl/issues).

Published PK/PD and quantitative systems pharmacology (QSP) models, rebuilt from their
papers in Julia ([ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl),
[OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl)), checked against the
papers' own figures and graded.

## Result

**{{HEADLINE}}**

{{TABLE}}

Each model page lists the targets, the per-curve errors, every fix with its evidence,
and the problems found in the paper. The [method](method.md) page explains the
digitising and scoring.

## First release

Six oncology models from *CPT: Pharmacometrics & Systems Pharmacology*, chosen to
span model types: tumour growth with angiogenesis (pazopanib), transient vascular
normalisation (antiangiogenic antibodies), translational antibody PK/PD with tumour
growth inhibition (pembrolizumab), a haematopoiesis systems pharmacology model of
bone-marrow toxicity (carboplatin), a QSP model of a masked tumour-activated antibody,
and CAR-T cellular kinetics (tisagenlecleucel).

Five are graded on typical-value curves. The CAR-T model is graded on a visual predictive
check: the authors' shaded intervals of simulated percentiles are compared with
percentiles from our own population simulation (see
[Population simulation](simulation.md)). Lines computed from observed data are never
grade-bearing.

## What it looks like

A model that reproduces: Ouerdani 2015, Figure 4, our simulation (blue) against points
digitised from the paper's own curves (black).

![Ouerdani 2015 reproduction](assets/models/ouerdani_2015/comparison.png)

A population figure that does not: Stein 2019, Figure 3a, percentile bands from 800
simulated replicates of the 90-patient study. The early expansion matches the authors'
bands; the upper tail and the persistence phase fall below them.

![Stein 2019 VPC](assets/models/stein_2019/fig3.png)

Population simulation: 500 virtual mice from the Ouerdani model with the paper's
between-subject variability and residual error. Simulating 10,000 subjects takes 1 to 4.3 s
on a standard 4-thread GitHub runner ([benchmarks](benchmarks.md)).

![500 virtual mice](assets/population_ouerdani.png)

## Using it

To rerun the atlas, clone the repository (each `reproduce()` writes its score file and
figures into its own model folder):

```bash
git clone https://github.com/singhharsh1708/QSPAtlas.jl
cd QSPAtlas.jl
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'   # reruns every model
```

```julia
using QSPAtlas
m = Module(); Base.include(m, joinpath(models_dir(), "ouerdani_2015", "reproduce.jl"))
rep = m.reproduce()        # simulates, scores every curve, writes score.toml and figures
grade_label(rep)           # "B†"
```

For the population simulation layer on your own models, add the package:

```julia
using Pkg
Pkg.add(url = "https://github.com/singhharsh1708/QSPAtlas.jl")
```
