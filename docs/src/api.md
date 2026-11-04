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
ErrorModel
AdditiveError
ProportionalError
CombinedError
ExponentialError
Subject
PopulationModel
SimResult
simulate
zero_re_prediction
virtual_population
percentiles
vpc
VPCResult
spaghetti
bands
plot_vpc
grade_label
from_target
figure_grades
```

## Internals

Helpers that are not exported but are documented for contributors.

```@autodocs
Modules = [QSPAtlas]
Public = false
```
