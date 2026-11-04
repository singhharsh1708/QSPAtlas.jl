# Contributing a model


Every model folder `models/<author_year>/` has the same files. The test suite and
the docs build rely on this layout.

```
models/<author_year>/
├── model.jl       # MTK system(s), parameter values and initial conditions, as printed in the paper
├── protocol.jl    # dosing, simulation settings, targets, tolerances, FIXES, benchmark_setups()
├── reproduce.jl   # defines reproduce(); regenerates the paper's figures and writes score.toml
├── reference/     # digitised points, one CSV per curve, plus digitise.toml
├── figures/       # comparison.png and regenerated paper figures (committed); digitise_check.png is local only
├── score.toml     # written by reproduce(); grade and per-curve metrics
└── report.md      # grade, targets, fixes, issues found, notes
```

## Rules

1. **Model as printed.** `model.jl` uses the equations and parameter values exactly as
   published (main text, tables, supplement). Every parameter has a comment giving
   its source (`# Table 1`, `# Eq. 4`, `# supplement S2 line 40`).
2. **Targets are frozen before comparing.** Write `reference/*.csv` and the
   tolerances in `protocol.jl` before looking at how the simulation compares.
   Loosening a tolerance afterwards must be justified in `report.md`.
3. **Digitise authors' simulation lines**, not observations, wherever the figure has
   them (`kind: curve`). Observed data may be added as `kind: points` with a looser,
   stated tolerance, but grade-bearing targets should be simulation lines.
4. **Fixes need evidence.** A `Fix` records what was changed and where the evidence
   comes from (a supplement file, an equation that is dimensionally inconsistent, a
   figure that only matches one reading). Re-estimating parameters is never a fix.
5. **Grades.** A: as-published run matches all targets. B: matches only with fixes.
   C: runs but does not match. D: cannot run from published information (set
   `blocker`). Unfinished work is not D.
6. **Self-contained reproduce.jl.** It must work when included into a fresh module:
   `m = Module(); Base.include(m, "models/<id>/reproduce.jl"); m.reproduce(plots = false)`.

## reproduce.jl skeleton

```julia
using QSPAtlas
using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

Base.include(@__MODULE__, joinpath(@__DIR__, "model.jl"))  # plain include is undefined in a bare Module()
Base.include(@__MODULE__, joinpath(@__DIR__, "protocol.jl"))

const ID = "ouerdani_2015"
const CITATION = "Ouerdani A, et al. CPT Pharmacometrics Syst Pharmacol 2015;4:660-668"
const DOI = "10.1002/psp4.12001"
const PMCID = "PMC4716582"

function reproduce(; plots = true, figdir = joinpath(@__DIR__, "figures"))
    refs = sort(collect(values(load_references(joinpath(@__DIR__, "reference")))); by = r -> r.id)
    published = run_curves(ref -> simulate_curve(ref; fixes = Fix[]), refs, "as published";
        tolerance = TOLERANCE, tolerances = TOLERANCES)
    fixed = isempty(FIXES) ? nothing :
        run_curves(ref -> simulate_curve(ref; fixes = FIXES), refs, "with fixes";
            fixes = FIXES, tolerance = TOLERANCE, tolerances = TOLERANCES)
    rep = Reproduction(; id = ID, citation = CITATION, doi = DOI, pmcid = PMCID, published, fixed)
    write_score(joinpath(@__DIR__, "score.toml"), rep)
    if plots
        plot_reproduction(rep, joinpath(figdir, "comparison.png"))
        # optional: regenerate the paper's figure in its own layout
    end
    return rep
end

(abspath(PROGRAM_FILE) == (@__FILE__)) && display(reproduce())
```

`simulate_curve(ref; fixes)` returns `(x, y)` on a dense grid (use 500+ points or
`saveat` that resolves dosing spikes) in the figure's units for reference `ref.id`.

## Reference CSV

```
# id: fig4_left_P
# figure: Figure 4 left, preclinical typical tumour volume P
# x: time [day] linear 0 75
# y: tumour volume P [mm^3] linear 0 250
# kind: curve
# source: PSP4-4-660-g004.jpg (PMC4716582)
# method: digitise.jl colour trace ...
x,y
0,150
...
```

Axis limits in `x:` and `y:` are the limits drawn in the figure; errors are measured
as fractions of that box. Default tolerance: RMS 0.015 and 95th percentile 0.03 of
the box.

## Digitising

`tools/digitise/digitise.jl` traces a curve by colour inside a pixel box using a
two-tick calibration per axis, and records the settings in `reference/digitise.toml`
so the points can be regenerated. The overlay `figures/digitise_check.png` draws the
points on the published figure image, so it is generated locally and never committed
or published (it contains the paper's figure). Manual point picking is allowed when tracing is
impossible (overlapping black curves, dashed lines), as long as `method:` says so
and the pixel coordinates are kept in `digitise.toml`.
