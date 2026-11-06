# Method

## How models are built

The models are rebuilt by AI coding agents. The agents transcribe each paper's equations,
parameter values and dosing into code, digitise the target figures, run the simulations and
comparisons, and draft the reports. The author reviews every equation, parameter value, fix
and grade against the paper before it is committed. The second implementations on the
[validation page](validation.md) were written by separate agent sessions with no access to
the atlas code. The next phase, a pre-registered study of a random sample of CPT:PSP models,
adds blind rebuilds by outside pharmacometricians, to compare with the agent rebuilds.

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

**‡** marks a B (or B†) that rests on at least one fix made after the first comparison
with the target figure (`stage = :after_comparison`). Such a fix was chosen with the
mismatch in view, so it carries more risk of fitting the figure than a fix decided
before any comparison, even when its evidence is independent of the figure. Each fix's
stage is listed on the model page.

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
corresponds to roughly 1.3 to 3.5 px on the vertical axis.

**How the tolerance was measured.** A digitising error study
([`tools/digitise/error_study/`](https://github.com/singhharsh1708/QSPAtlas.jl/tree/main/tools/digitise/error_study))
drew about 7,650 known curves at the real panel sizes of the atlas figures, saved them as
JPEG at the quality of the PMC files (92), traced them with the atlas tools and rules, and
scored them against the known curves with the atlas metric. For traces that pass an
overlay check, the 99th percentile of the error ranges from 0.0013 to 0.0077 (RMS) and
from 0.0029 to 0.0157 (P95) across panel types; the worst type is a small log-scale panel
(180 x 149 px). The atlas references were also compared with the same curves taken from
the first authors' theses (vector or higher-resolution figures) for 34 curves: worst RMS
0.0080, worst P95 0.0107. The rule "twice the worst panel type's 99th percentile, at two
significant figures" gives RMS 0.015 and P95 0.031; the atlas keeps 0.015 and 0.03. The
factor of two covers the gap between real and synthetic errors (real errors were 0.9 to
1.4 times the synthetic medians) and error sources neither part measures (publisher
resampling, calibration slips). No committed result changes under the measured rule.

About 4 % of synthetic traces failed grossly (the tracer jumped to another curve or mixed
two lines). Only a visual overlay of the traced points on the figure catches this, so
every reference trace is checked that way before it is frozen. On a log axis spanning many
decades the same tolerance allows a larger relative error (on a six-decade axis, RMS 0.015
corresponds to a uniform factor of about 1.23); the [sensitivity of every grade to the
tolerance](sensitivity.md) is published. The 95th percentile instead of the maximum limits
the effect of a single badly traced point, but does not remove it: one point 0.12 of the
axis away fails the RMS criterion on its own for a 60-point curve. The maximum is
reported. A model may set a different tolerance for a curve, with the reason stated before
comparing.

The nearest-distance metric scores only the regions where the figure has digitised
points: a simulated feature with no reference points nearby costs nothing. Every curve
therefore reports diagnostics alongside its pass or fail (`score.toml` and the model
pages): `coverage`, the share of the axis x range spanned by the digitised reference;
`sim_coverage`, the share where the simulation can be placed on the axes; a reverse
distance (`rev_rms`, `rev_p95`) from the simulated curve back to the reference; and for
bands, the deviation from the band centre in band half-widths and, where a model supplies
its own simulated interval edges, the error of those edges. Some passes rest on little of
the axis (Lindauer's base case is two points; Stroh `fig4c_mAb` covers 11 % of its axis).

These diagnostics are **report-only for the first six models**: they were added after
those models were graded, and turning them into pass criteria now would change grades
after the results were known. Every curve also reports its `fraction_checked`: the arc
length of the simulated curve inside the axis box over the x range of the digitised
reference, divided by its whole arc length inside the box, measured in the normalised axis
box (log axes in decades). It shows how much of the simulated curve a pass rests on and is
report-only for every model; it is left out for bar charts. Phase 3 models also apply the
reverse-distance criterion: a curve reference passes only if, in addition to the tolerance
above, its whole-trace reverse-distance 95th percentile (`rev_p95`) is at most **0.31** of the
axis box, as fixed in Amendment 1 of the Phase 3 protocol. It is set per model with
`Tolerance(reverse_p95_max = 0.31)`, applies to curve references only (never bands, bars or
points), counts as failed when no simulated point lies inside the reference's x range, and is
recorded as `failed_reverse` in `score.toml`. Coverage stays report-only.

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
