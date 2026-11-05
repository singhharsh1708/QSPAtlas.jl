# Method

## What counts as reproducing a model

A paper's model is reproduced when the equations and parameter values it publishes,
simulated under the protocol it describes, give the curves it shows. The atlas
checks this against the authors' own simulation lines, not against observed data:
a model can fit its data badly and still be perfectly reproducible, and the reverse.

Each model is graded:

| grade | meaning |
|---|---|
| **A** | the model as published runs and matches every target curve |
| **B** | it matches only after documented fixes (typos, unstated choices, text vs code) |
| **C** | it runs but does not match, with or without fixes |
| **D** | it cannot be run from the published information |

**B†** marks a B that needed an input the paper never states (an initial condition,
or which of several published scenarios a figure shows) and that could only be read
from the figure being reproduced. That makes the match circular for that one value,
so the dagger flags it; the model's equations and parameters are still tested over
the whole time course. Parameters read off the target figure are never allowed: that
model is D. Every fix must cite its evidence, or it cannot earn a B.

**What may fill a gap.** A value the paper leaves out may be supplied only from the
paper itself (another section, table or figure that is not a scored target), its
supplement, a reference the paper cites for that value, or a standard physical or
physiological constant (a molecular weight, IgG mass of about 150 kDa, human blood
volume of about 5 L). Each such fill is a documented `:missing` fix. A model is D
only when an essential equation or value has none of these sources. Values chosen
by matching the scored figure are `:from_target` (B†) for protocol inputs and never
allowed for parameters. This rule was written after reviews had argued for D on three
of the models (Fornari, Lindauer, Stroh), and it is what keeps them at C. To show its
effect, every model also has a **strict grade**: under the strict reading only the
article and its supplement count as published information, and a model that needs any
value from elsewhere (a cited reference, a standard value, an assumption or the target
figure) is D. Both grades are listed on each model page and in the
[sensitivity tables](sensitivity.md).

**Provenance of every choice.** Each fix records its basis (paper, supplement,
interpretation, cited reference, standard value, assumption, target figure) and its
stage: whether it was chosen before or after the simulation was first compared with
the target. A published run that needs gap fills lists them and can never be graded A.

Grades describe the paper, not the effort spent. Re-estimating parameters never
counts as a fix, and recalibration results are reported separately without changing
a grade. Unfinished work is not graded D.

## Targets

For each paper, the grade-bearing targets are the simulation lines in chosen figure
panels. They are chosen, digitised and given tolerances before the simulation is
compared with them, and the choice is recorded in the model's notes and report.

## Digitising

Figures are taken from the open-access article (PubMed Central). `tools/digitise`
traces a curve by colour inside a pixel box and maps pixels to data with a two-tick
calibration per axis, linear or logarithmic. The settings are saved in each model's
`reference/digitise.toml`, so every point can be regenerated from the source image.
Where colour tracing cannot separate curves (black solid and dashed lines that
overlap), points are picked from pixel inspection and the method line says so. Every
digitisation is checked by overlaying the points on the original image
(`figures/digitise_check.png`) and against an axis tick not used for calibration.

## Scoring

Both axes of the figure panel are scaled to `[0, 1]` (log axes in decades), using the
axis limits as drawn. For each digitised point the error is its shortest distance to
the simulated curve in that box. Distance to the curve, rather than vertical error,
keeps steep sections (dosing spikes, rapid expansion) from being penalised for tiny
horizontal digitising error.

A curve passes when the root mean square error is at most **0.015** and the 95th
percentile is at most **0.03** of the axis box. The plot panels digitised so far are 85
to 294 px wide or high, so one pixel is 0.3 to 1.2 % of an axis, and the tolerance
corresponds to roughly 1.3 to 3.5 px on the vertical axis. Careful tracing of the
cleanest curves scores about 0.001 to 0.002. The tolerance is a judgement, not yet a
measured digitising error distribution; a measured error study is planned, and the
[sensitivity of every grade to the tolerance](sensitivity.md) is published. On a log
axis spanning many decades the same tolerance allows a larger relative error (on a
six-decade axis, rms 0.015 corresponds to a uniform factor of about 1.23). The 95th
percentile instead of the maximum means one badly traced point cannot fail a curve on
its own; the maximum is still reported. A model may set a different tolerance for a
curve, with the reason stated before comparing.

The nearest-distance metric scores only the regions where the figure has digitised
points: a simulated feature with no reference points nearby costs nothing. Per-curve
coverage of the x range is part of the planned metric upgrades.

## Population figures

Figures that show population spread (visual predictive checks) are graded only on the
authors' simulated quantities: shaded intervals of simulated percentiles are `band`
references, scored as zero inside the interval and as the vertical distance to the nearer
edge outside it; simulated percentile lines use the `vertical` metric, compared at matched
times. Lines computed from observed data are never grade-bearing. The cohort design, bins,
quantile definition, residual error and number of replicates are fixed before comparing,
and replicates are added until the Monte Carlo error is small next to the tolerance. See
[Population simulation](simulation.md).

## Regression testing

`test/runtests.jl` runs every `models/*/reproduce.jl` in a fresh module and checks
that the grade and every curve's pass/fail equal the committed `score.toml`. A change
in a solver, ModelingToolkit or the model code that moves a result fails CI.
