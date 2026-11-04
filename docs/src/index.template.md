# QSP Reproducibility Atlas

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

## Using it

```julia
using Pkg
Pkg.add(url = "https://github.com/singhharsh1708/QSPAtlas.jl")
using QSPAtlas
m = Module(); Base.include(m, joinpath(models_dir(), "ouerdani_2015", "reproduce.jl"))
rep = m.reproduce()   # simulates, compares, writes score.toml and figures
```
