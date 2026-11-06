# API

```@docs
QSPAtlas
Axis
ReferenceCurve
load_reference
load_references
write_reference
Tolerance
CurveComparison
compare_curve
CurveDiagnostics
curve_diagnostics
Fix
CurveResult
Run
run_curves
Reproduction
grade
write_score
models_dir
model_ids
bolus_callback
repeated_doses
plot_reproduction
BenchmarkSetup
benchmark_problem
```

## Population simulation

```@docs
Route
Dose
Regimen
DoseControl
DoseEvent
doses
resolve
dosing_callback
Reset
reset_times
steady_state_cycle
with_steady_state
SteadyStateOptions
SteadyStateError
NMData
read_nmdata
write_nmdata
ErrorModel
AdditiveError
ProportionalError
CombinedError
ExponentialError
Subject
occasion_at
PopulationModel
SimResult
simulate
zero_re_prediction
virtual_population
percentiles
vpc
VPCResult
StratifiedVPC
mc_check
failure_rate
spaghetti
bands
plot_vpc
grade_label
grade_strict
NotRunnable
from_target
post_comparison
figure_grades
```

## Internals

Helpers that are not exported but are documented for contributors.

```@autodocs
Modules = [QSPAtlas]
Public = false
```
