# Population simulation

QSPAtlas includes a simulation layer in the style of NONMEM and other pharmacometric
tools, built on ModelingToolkit and OrdinaryDiffEq. It covers dosing regimens, virtual
subjects with covariates, between-subject variability, residual error, cohort simulation
and simulation-based visual predictive checks.

## Dosing

A model receives doses through `Route`s. A route names a compartment, the state that
receives bolus amounts and, for infusions, a rate parameter declared in the model:

```julia
@parameters CL = 5.0 V = 50.0 Rin = 0.0
@variables A(t) = 0.0
@mtkcompile sys = System([D(A) ~ -CL / V * A + Rin], t)
routes = [Route(:central, sys.A; rate = sys.Rin)]
```

`Regimen(amt; time, cmt, ii, addl, rate, duration, ss)` follows the NONMEM dose record:
`amt` at `time`, then `addl` more doses every `ii`. Regimens concatenate with
`Regimen(r1, r2)`. `DoseControl(bioav = Dict(:depot => 0.6), lag = Dict(:depot => 0.5))`
holds subject-specific bioavailability and lag time. `rate = :model` and `duration = :model`
(NONMEM RATE=-1 and -2) take the infusion rate or duration from the subject's
`DoseControl(rate = Dict(:central => 20.0), duration = Dict(:central => 5.0))` (NONMEM R1 and
D1); a modelled rate keeps the rate and F scales the duration, a modelled duration keeps the
duration and F scales the rate, as for entered values. A dose whose `DoseControl` lacks the
modelled value is an error. Values that carry random effects, covariate effects or
between-occasion variability are computed by the population model instead, per dose record
(see [Dose control from the model](#Dose-control-from-the-model)).

| situation | behaviour |
|---|---|
| bolus | `F * amt * route scale` is added to the route's state |
| fixed-rate infusion (`rate` given) | rate kept, F scales the duration |
| fixed-duration infusion (`duration` given) | duration kept, F scales the rate |
| lag time | shifts the start; an infusion keeps its length |
| overlapping infusions | the rate parameter is set to its baseline plus the sum of the infusions running at that moment, recomputed at every start and stop; the baseline is restored exactly after the solve |
| doses at the same time | grouped only when times are exactly equal; applied in the order infusion stops, infusion starts, boluses; all boluses at one time apply |
| events at the first time point | applied before the first solver step |
| events at the last time point | applied |
| infusion running at the first time point | switched on if it started before and stops after that time |
| two routes sharing one rate parameter | rejected |

### Steady state

`Regimen(...; ss = true, ii = ...)` starts the subject at the periodic state of that
regimen. The state at the start of a dosing interval is found by simulating repeated
intervals from the model's initial conditions until every state changes by at most
`atol + rtol * abs(x)` between intervals (`SteadyStateOptions`; defaults rtol 1e-8,
atol 1e-10, at most 1000 intervals). Doses whose lag pushes them past `ii` wrap into the
next interval, and infusions longer than `ii` are carried across intervals and into the
main simulation.

An SS=1 record replaces the state at its time, so everything given before it is superseded:
doses of earlier records that have not entered the model yet are dropped, running infusions
are stopped, and the `addl` doses of earlier records that are still to come are cancelled,
including those whose nominal time is after the SS=1 time. The time of the record decides,
not the time of each `addl` dose; the SS=1 record's own `addl` doses and later records are
kept. This follows the NONMEM help on SS ("infusions in progress or pending additional doses
are cancelled") and has not been compared with NONMEM output. mrgsolve 2.0.1 keeps those
`addl` doses (see [Validation](validation.md)).

`ss = 2` (NONMEM SS=2) adds the steady state of that dose to the current state instead of
replacing it, so regimens can be combined, for example 10 every morning and 15 every evening:
`Regimen(Regimen(10.0; time = 8, ii = 24, ss = 1, cmt), Regimen(15.0; time = 20, ii = 24, ss = 2, cmt))`.
The added amount is the periodic state of the SS=2 dose minus the periodic state without it
(both from the initial conditions). For linear models with zero initial conditions this is the
NONMEM definition ("the sum of the steady-state amounts resulting from the given dose plus
whatever those amounts would have been at the event time were the steady-state dose not
given"); subtracting the dose-free state keeps an endogenous baseline (a turnover response)
from being counted twice. Superposition is only exact for linear kinetics, as NONMEM also
notes. An SS=2 record cancels nothing; each SS=2 record has its own time and `ii`. A regimen
still has at most one SS=1 time.

A constant infusion at steady state (NONMEM AMT=0, RATE>0, SS=1 or 2) is
`Regimen(0.0; time, rate = R, ss = true, cmt)` with `ii = 0`: the state at `time` is the
equilibrium under the infusion (for a one-compartment model the concentration is
`R / CL`). As in NONMEM, the record marks the END of an infusion that started at minus
infinity: after `time` the infusion is off. To keep it running, add an ordinary infusion from
`time`, as in the NONMEM Users Guide example (`TIME 0, AMT 0, RATE 30, SS 1` followed by
`TIME 0, AMT 3000, RATE 30`). Bioavailability does not apply to it (NONMEM) and lag has no
effect on it. The equilibrium is found by the same iteration with cycle lengths 1, 2, 4, ...
time units (capped at 2^40), with the same convergence test.

Where the iteration starts is a deliberate choice: the model's initial conditions. mrgsolve
starts from the state at the steady-state record; how NONMEM initialises its steady-state
computation depends on the routine and has not been checked against NONMEM output here.
The starting point does not matter when the periodic state is unique (linear PK); it can
matter for models with endogenous baselines (cell counts, receptors) or several periodic
solutions. Steady state is refused for
models whose equations depend on time explicitly (for example a resistance term
`exp(-δ t)`), because a periodic state does not exist there; `SteadyStateOptions(allow_time_dependent = true)`
overrides this when the time dependence is known to be periodic with `ii`. If no steady
state is reached, a single-subject simulation throws `SteadyStateError`; in a cohort the
subject is kept with `retcode = :SteadyStateFailed` and empty observations, and the rest
of the cohort is simulated.

### Resets

`Reset(time)` returns every state to its initial conditions and cancels running infusions
and lagged doses that have not yet been given (in the spirit of NONMEM EVID 3); later
`addl` doses are kept. The treatment of lagged doses pending at a reset is our convention
and has not been checked against NONMEM output.
A reset followed by a dose at the same time (EVID 4) applies the reset first, then the dose
with its own lag. NONMEM resets the compartment amounts to zero; QSPAtlas resets to the
model's initial conditions (the same for models that start at zero).

### Observations at event times

By default an observation at the same time as a dose, an infusion change, a reset or a
covariate change sees the state after the event (mrgsolve `recsort = 3`). An observation can
be tagged `:pre` instead: `Subject(...; obstimes = [12.0, 12.0], obsorder = [:pre, :post])`
gives the trough just before the dose at 12 and the value just after it. A `:pre` observation
sees the state before every dose, reset and steady-state event at its time. Covariate values
at that time follow the covariate rule either way (below). [Datasets](#NONMEM-style-datasets)
set the tags from the record order.

### Validation

The dosing semantics above were checked case by case against mrgsolve; see
[Validation](validation.md).

## Subjects and population models

```julia
sub = Subject("1"; regimen, covariates = (WT = 70.0,), obstimes = collect(0:1.0:48),
              timevarying = [24.0 => (WT = 60.0,)])
```

Covariates are constant or piecewise constant. With the default `covariate_rule = :locf`
(last value carried forward) a change at time `t` re-evaluates the model parameters from `t`
on, and an observation at `t` sees the new value. With `covariate_rule = :nocb` (next
observation carried backward, the NONMEM convention for data records) the values listed at
`t` apply over the interval ending at `t`, `(previous change, t]`: an observation at `t` still
sees them and the next values start just after `t`; before the first entry its values apply,
after the last entry they remain. Infusion rate parameters belong to the dosing code and may
not be set by covariates. `endpoints = [...]` restricts each observation to one endpoint.

```julia
model = PopulationModel(sys; routes,
    θ = (CL = 5.0, V = 50.0),
    ηnames = [:CL, :V], Ω = [0.09 0.0; 0.0 0.04],
    params = (θ, η, c) -> [sys.CL => θ.CL * exp(η.CL) * (c.WT / 70)^0.75, sys.V => θ.V * exp(η.V)],
    observed = (conc = sys.A / sys.V,),
    error = (conc = CombinedError(0.1, 0.15),),
    tspan = (0.0, 48.0), solver = Tsit5())
```

`Ω` is the covariance of the random effects in the order of `ηnames`. It must be symmetric
positive semidefinite within a tolerance relative to its largest entry; singular matrices
are allowed and sampled exactly. Error models (σ are standard deviations; negative
observations are never clipped):

| error model | observation |
|---|---|
| `AdditiveError(σ)` | `ipred + σ ε` |
| `ProportionalError(σ)` | `ipred (1 + σ ε)` |
| `CombinedError(σadd, σprop)` | `ipred + sqrt(σadd² + (σprop ipred)²) ε` |
| `ExponentialError(σ)` | `ipred exp(σ ε)` (median-centred) |

## Dose control from the model

Bioavailability (NONMEM F), lag time (ALAG), modelled infusion rate (R, for `rate = :model`,
RATE=-1) and modelled duration (D, for `duration = :model`, RATE=-2) can be individual
parameters with between-subject variability and covariate effects. The model computes them
from `(θ, η, covariates)` like the other parameters:

```julia
model = PopulationModel(sys; routes, θ = (CL = 5.0, V = 50.0, F = 0.6, ALAG = 0.5),
    ηnames = [:CL, :F, :ALAG], Ω = Diagonal([0.09, 0.04, 0.16]),
    params = (θ, η, c) -> [sys.CL => θ.CL * exp(η.CL), sys.V => θ.V],
    dosecontrol = (θ, η, c) -> (bioav = Dict(:depot => θ.F * exp(η.F) * (c.WT / 70)^0.2),
                                lag = Dict(:depot => θ.ALAG * exp(η.ALAG))),
    observed = (conc = sys.A / sys.V,), error = (conc = ProportionalError(0.1),),
    tspan = (0.0, 48.0), solver = Tsit5())
```

`dosecontrol` returns a NamedTuple with any of `bioav`, `lag`, `rate`, `duration` (each a
`Dict(compartment => value)`) or a `DoseControl`. The values follow the rules of the
[Dosing](#Dosing) table (F scales a bolus, the duration of a fixed-rate infusion or the rate
of a fixed-duration infusion; lag shifts the start). A subject's own `DoseControl` is combined
with them; the same entry in both is an error.

**Which covariate value a dose uses.** `dosecontrol` is evaluated once per dose record, at the
time of the record, with the covariates in force there by the subject's covariate rule. For a
dataset dose record this is the value on that record (under NOCB a record's values apply up to
and including the record; under LOCF from the record on). The values then stay with the dose:

| situation | value used |
|---|---|
| covariate changes at the dose record | the record's value |
| lagged dose, covariate changes between the record and the moment the dose enters | the record's value |
| infusion running when a covariate changes | the record's value until the infusion ends |
| `addl` doses | the values of their record, also when covariates or the occasion change in between |
| steady-state dose (SS=1, SS=2) | the record's values for the whole periodic history |

This is the NONMEM convention for the dose record: the NONMEM help on dose event records
says "the value of ALAGn in effect at time T applies to the dose and to all subsequent
additional doses", and that later changes of F "have no affect on past doses"
([nmhelp](https://nmhelp.pages.dev/dose.html)). For F of an `addl` or lagged dose NONMEM uses
"the bioavailability parameter in effect at the time the dose is actually introduced", which
with the default (PK not called at additional or lagged dose times) depends on which event
records surround the dose; this case has not been checked against NONMEM output. mrgsolve
2.0.1 uses the dose record's F and ALAG for its `addl` and lagged doses, as QSPAtlas does
(checked on the [validation page](validation.md)). To give implied doses their own values,
write them as separate dose records.

## Between-occasion variability

Occasions split a subject's record into periods (dosing days, study periods) whose parameters
vary at random around the subject's values. A second random-effect vector `κ ~ N(0, Ω_IOV)`
is drawn independently for every occasion of every subject:

```julia
model = PopulationModel(sys; routes, θ, ηnames = [:CL], Ω = [0.09;;],
    κnames = [:CL], Ω_IOV = [0.0225;;],
    params = (θ, η, c, κ) -> [sys.CL => θ.CL * exp(η.CL + κ.CL), sys.V => θ.V],
    observed, error, tspan, solver)
sub = Subject("1"; regimen, obstimes, occasions = [0.0, 24.0, 48.0])        # occasions 1, 2, 3
```

With `κnames`, `params`, `init` and `dosecontrol` are called with a fourth argument, the κ of
the occasion in force (a function with only three arguments is called without it; at least one
of `params` and `dosecontrol` must take κ). `Ω_IOV` is validated like `Ω`.

  - **Occasions** are given on the subject as start times (labelled 1, 2, ...) or as
    `time => label` pairs, or by an occasion column in a dataset
    (`read_nmdata(...; occasion = :OCC)`). A label that comes back (1, 2, 1) is the same
    occasion with the same κ, as with NONMEM's occasion indicator variables. Without occasions
    the whole record is occasion 1.
  - **Timing** follows the subject's covariate rule, like a time-varying covariate. Under LOCF
    (the `Subject` default) an occasion starts at its time: an observation or a dose at that
    time is in the new occasion. Under NOCB (datasets) the occasion on a record applies over the
    interval ending at that record, so the new occasion's parameters start just after the last
    record of the previous occasion, as when NONMEM evaluates its PK block with each record's OCC.
    Records at one time must have the same occasion.
  - **Switching.** Parameters are re-evaluated at every occasion change (a callback, as for
    covariates); the state is continuous. Dose control is evaluated per dose record (above), so
    doses given in the new occasion use its κ, a dose given before keeps the old one, and `addl`
    doses keep the κ of their record.
  - **Random streams.** κ of occasion `j` comes from its own stream keyed by (seed, replicate,
    subject id, `j`), separate from the η and residual streams. Results therefore do not depend
    on threads or subject order, η is the same as without IOV, and adding an occasion does not
    change the κ of the others. `simulate(...; κ = :zero)`, a NamedTuple for every occasion or a
    `Dict(label => κ)` fixes them; `η = :zero` (and `zero_re_prediction`) also sets κ to zero.
    With IOV the subject rows of the result have a field `κ` (label => κ).
  - Steady-state doses use the occasion in force at the steady-state time for the whole
    periodic history.

## NONMEM-style datasets

```julia
data = read_nmdata("pk.csv"; cmt = Dict(1 => :depot, 2 => :central),
                   endpoints = Dict(2 => :conc))          # covariate_rule = :nocb by default
res = simulate(model, data.subjects; seed = 1)
write_nmdata("pk_sim.csv", data; result = res)          # DV replaced record by record
```

`read_nmdata` reads a CSV file (comma- or whitespace-separated, "." or empty for missing), a
vector of NamedTuples or a NamedTuple of columns, and returns an `NMData` with one `Subject` per
individual record. Every column that is not a data item below is a covariate (or listed in
`drop`). The NONMEM conventions used are those of the NONMEM Users Guide data item
descriptions (EVID, MDV, AMT, RATE, SS, II, ADDL, CMT), cited from the copy at
[nmusers.github.io](https://nmusers.github.io/docs/user-guide/nmtran/data/).

| item | read as |
|---|---|
| `ID` | a new individual where the ID changes; an ID that comes back after another ID is an error (NONMEM would make it a new individual) |
| `TIME` | must not decrease within an individual (NONMEM allows a decrease after a reset; not supported) |
| `EVID` | 0 observation, 1 dose, 2 other-type, 3 reset, 4 reset and dose; other values are errors. Without `EVID` (NONMEM help, EVID): a dose if `AMT`, `RATE` or `SS` is nonzero (so a steady-state infusion, `AMT` 0, `RATE` > 0, `SS` 1, is a dose), else other-type if `MDV` is 1, else an observation |
| `MDV` | must be 0 on observation records and 1 on all others (NONMEM); `DV` of non-observation records is not used |
| `AMT`, `RATE` | a dose needs `AMT` > 0, except a steady-state infusion (`AMT` 0, `RATE` > 0 or -1, `SS` 1 or 2, `II` 0); `RATE` 0 bolus, > 0 rate, -1 modelled rate, -2 modelled duration. Non-dose records must have 0 or missing |
| `SS`, `II`, `ADDL` | `SS` 0, 1 or 2 (3 is an error); `II` > 0 with `SS` or `ADDL`, and an error on a dose without either; `SS` on an EVID 4 record is an error |
| `CMT` | dose compartment through `cmt` (a text value without an entry is the route name); observation endpoint through `endpoints` (or `endpoints = :all`); 0 or missing on an observation means every endpoint; negative values (turning a compartment off) are errors |
| covariates | numbers, never missing (NONMEM would read a missing value as 0); records at the same time must agree |
| `PCMT`, `CALL`, `CONT`, `L1`, `L2`, `DATE`, `DAT1-3`, `XVID1-5`, `MRG_`, `RAW_`, `RPT_` | errors |

Comment records (starting with `#` or `@`) are errors rather than skipped. Nothing in a
dataset is ignored silently: anything that cannot be represented is an error that names the
record.

**Record order.** Records keep their order. An observation record placed before a dose,
reset or steady-state record at the same time is a `:pre` observation (it sees the state
before that record; NONMEM processes records in order), one placed after it is `:post`.
Doses implied by `ADDL` at the time of an observation record are given before the
observation, as in mrgsolve's default (`recsort = 1`, which "puts doses scheduled through addl
before observations at the same time"); the NONMEM documentation does not say how an implied
dose and an event record at the same time are ordered, so this has not been checked against
NONMEM. Orders the simulation cannot express are errors: an observation between two dose
records at one time, an observation before a dose record at a time where an `ADDL` dose also
falls, and a reset or SS=1 record after a dose record at the same time (by the NONMEM help on
EVID and SS that dose would be reset to zero; the simulation applies resets and steady states
first).

**Covariates.** `covariate_rule = :nocb` is the default for datasets and follows NONMEM: the
values on a record are used to advance the system from the previous record to this one, so
they apply over the interval ending at the record (and at the record itself, for its own
prediction). `covariate_rule = :locf` applies them from the record onward. Dose control
(bioavailability, lag, modelled rates and durations) that depends on covariates, random effects
or the occasion comes from the model's `dosecontrol`, evaluated on each dose record with that
record's covariates (see [Dose control from the model](#Dose-control-from-the-model)). The
`control` keyword still takes a constant `DoseControl` or a function
`(id, baseline covariates) -> DoseControl`.

**Occasions.** `read_nmdata(...; occasion = :OCC)` reads that column as the occasion of each
record for between-occasion variability (integer labels, never missing, equal on records at one
time), with the dataset's covariate rule (see [Between-occasion variability](#Between-occasion-variability)).
Without the keyword such a column is an ordinary covariate.

**EVID 2.** An other-type record changes nothing in the state. It carries covariate values
(under NOCB they apply over the interval ending at it), and it produces no output row:
NONMEM computes a prediction for every record, QSPAtlas only for observation records, so add
an observation record (with `DV` ".") where a prediction is wanted. `MDV` must be 1; `AMT`,
`RATE`, `SS`, `II`, `ADDL` must be 0 or missing. A positive `CMT` (turn a compartment on) has
no effect because every compartment is always on; a negative one is an error.

**Steady state.** SS=1, SS=2 and the steady-state infusion follow the rules of the
[steady state section](#Steady-state); `RATE = -1 / -2` the rules for modelled rates and
durations above.

`NMData.obs` lists the observation records with their endpoint, tag and `DV`.
`write_nmdata(path, data; result, replicate = 1, value = :dv)` writes the dataset back with the
original columns, `DV` replaced by the simulated values (matched record by record).

## Simulating

```julia
pop = virtual_population(500; regimen, obstimes = collect(0:1.0:48),
                         covariate_sampler = rng -> (WT = 50 + 40rand(rng),))
res = simulate(model, pop; replicates = 1, seed = 1)       # EnsembleThreads by default
typ = zero_re_prediction(model, sub)                         # η = 0, no residual error
```

Every (replicate, subject) pair has its own random stream derived from the seed, with
separate draws for random effects and residual error, so results do not depend on the
number of threads or on the order of subjects. A failed solve keeps its rows with `NaN`
values and its return code, and is reported once per cohort. Note that the zero-random-effect
prediction is generally not the population median of a nonlinear model.

## Summaries and visual predictive checks

`percentiles(res, :conc; q = (0.05, 0.5, 0.95), bins)` gives per-bin quantiles
(Hyndman-Fan type 7). `vpc(model, subjects; endpoint, q, bins, replicates, pi_level, seed)`
simulates `replicates` copies of the cohort and reports, per bin and quantile, the median of
the replicate quantiles, its `pi_level` interval and a bootstrap Monte Carlo standard error.
Values below a limit of quantification can be censored to the limit, dropped or kept
(`blq`); failed and censored values are counted in every row. `spaghetti`, `bands` and
`plot_vpc` draw the results.

## Population figures in the atlas

Figures that show population spread are graded only on the authors' *simulated*
quantities: the shaded intervals of simulated percentiles (`kind: band` references) or
simulated percentile lines (`kind: curve` with `metric: vertical`). Lines computed from
observed data are a diagnostic, never grade-bearing. The cohort design, bins, quantile
definition, residual error and number of replicates are fixed before the comparison, and
the number of replicates is raised until the Monte Carlo error is small next to the
tolerance; the tolerance itself is never widened for a noisy simulation.

<!-- VPC extensions: censored VPC, prediction correction, strata, weighting. Appended as its own section. -->

## VPC details: estimand, censoring, prediction correction, strata

### What the binned percentiles estimate

`percentiles` and `vpc` pool the values of a bin. With the default `weighting = :observation`,
a subject with more observations in a bin counts more times, in the simulated and in the
observed percentiles alike. This is the usual VPC and it answers "what is the distribution
of the observations falling in this bin?", not "what is the distribution over subjects?".
The two agree when everyone is sampled at the same times (`bins = :obstimes`) or has one
observation per bin; they differ with wide bins and unequal or informative sampling, where
the pooled percentile leans towards the subjects who were sampled most. `weighting = :subject`
gives each subject one value per bin and replicate, its `subject_value` (`:mean`, `:first` or
`:last` by time), so every subject weighs once. Simulated and observed data always use the
same weighting. The counts `n_failed` and `n_finite` stay counts of observations.

### Censored VPC

`vpc(...; blq = (lloq = x, handling = :censored))` follows the usual censored VPC.

  - The fraction of values below `lloq` is computed per bin (`VPCResult.blq_fraction`): observed,
    and simulated as the median and the `pi_level` interval of the per-replicate fractions. It is
    also computed with the other `handling` values. Observed values are censored when
    `dv < lloq`, or when the observed row has a `blq = true` field (for data that store
    censored records at the limit).
  - Quantiles keep the censored values as censored: they rank below every quantified value.
    The type 7 quantile interpolates two neighbouring order statistics; if the lower one is
    censored, the quantile is below the limit (or depends on a value known only to be below it)
    and is reported as `NaN` with `below_lloq = true`, never as a number. A replicate whose
    quantile is below the limit is the lowest value when the quantiles are summarised across
    replicates, so the median, `lo` and `hi` are `NaN` with their own flags
    (`below_lloq`, `lo_below`, `hi_below`) when their rank falls among those replicates;
    `n_cens` counts them. `mc_se` is computed with them set to `lloq`.
  - The older handlings (`:censor_to_lloq`, `:drop`, `:keep`) and their results are unchanged.

`plot_vpc(v; blq_fraction = true)` adds a panel with the fraction below the limit. M3 (Beal) is
a method for estimating a model from censored data, not a simulation option, and is not
implemented here.

### Prediction-corrected VPC

With `pc = true` every observed and simulated value is normalised (Bergstrand, Hooker,
Wallin and Karlsson, *Prediction-corrected visual predictive checks for diagnosing nonlinear
mixed-effects models*, AAPS J 13:143-151, 2011):

```math
\mathrm{pcY}_{ij} = Y_{ij}\,\frac{\mathrm{PRED}_{\mathrm{bin}}}{\mathrm{PRED}_{ij}}
```

`PRED_ij` is the population prediction of the observation (random effects zero, no residual
error, the covariates and doses of that subject; [`zero_re_prediction`](@ref)), and
`PRED_bin` the median of `PRED` over the design of the bin (each subject and time once,
taken from the simulated cohort and used for the observed data too, so both live on one scale).
`vpc(model, subjects; pc = true)` computes the predictions; for a `SimResult` pass
`pred = simulate(model, subjects; η = :zero, residual = false)` or a dictionary
`(id, time) => prediction`; observed rows may carry their own `pred` field. The correction
is proportional, so it needs positive predictions: values with a prediction that is not
positive and finite are left out and counted in `n_nopred`. The lower-bound term and the
log-scale form of the paper are not implemented. BLQ censoring is decided on the raw
value, before the correction.

`variability_correction = true` adds the paper's variability correction (pvcVPC):

```math
\mathrm{pvcY}_{ij} = \mathrm{PRED}_{\mathrm{bin}} + (\mathrm{pcY}_{ij} - \mathrm{PRED}_{\mathrm{bin}})\,\frac{\mathrm{sd}_{\mathrm{bin}}}{\mathrm{sd}_{ij}}
```

where `sd_ij` is the standard deviation, over the simulated replicates, of `pcY_ij` of that
observation, and `sd_bin` the median of the `sd_ij` in the bin. It needs at least two replicates
and every observed record must be in the simulated design.

For a model with dose-proportional kinetics the correction removes the dose: the quantiles of
`pcY / PRED_bin` of a mixed-dose cohort and of each dose group coincide, which the tests use
as a known answer.

### Strata

`strata = subject -> key` (or a dictionary `id => key`) gives one VPC per stratum, a
`StratifiedVPC` indexed by the key. The function receives the `Subject` (for a `SimResult` the
subject row, which has `id` and `covariates`). One simulation serves all strata; each stratum
is the VPC of its own subjects and its own observed rows (matched by `id`), with its own bins'
counts and, for a pcVPC, its own `PRED_bin`. `plot_vpc(strat)` draws one panel per stratum.

### Records, censoring limits under correction, subject means

These rules refine the subsections above; where they say "each subject and time", read "each
record".

  - **Records.** Predictions and the variability groups of a pcVPC are matched per record:
    (id, time, k), where k numbers the records of one subject at one time in their order
    (`:pre` before `:post`, as `Subject` stores them). A pre-dose and a post-dose record at a
    dose time are two records with two different `PRED`, and `PRED_bin` takes each record of the
    design once. Observed rows at one time are numbered the same way: `order = :pre` rows first
    when they have an `order` field (`read_nmdata` gives one), otherwise in the order given.
    `pred` may be a dictionary keyed by `(id, time, k)`; a key `(id, time)` applies to every
    record at that time and cannot tell pre-dose from post-dose. With `vpc(model, subjects;
    pc = true)` an observed record outside the cohort design is predicted from its subject
    observed at the observed records of that time, with their `:pre`/`:post` order.
  - **Typical predictions.** `PRED` is the typical population prediction of Bergstrand et al.
    (2011): η = 0 **and** κ = 0 (no between-occasion effects), no residual error.
    `vpc(model, subjects; pc = true)` uses η = 0 and κ = 0 for `PRED` whatever `η` or `κ` the
    simulation is given. For `vpc(r::SimResult; pred = ...)` a `SimResult` with nonzero η or κ,
    or simulated with residual error, is refused; pass
    `simulate(model, subjects; η = :zero, κ = :zero, residual = false)`.
  - **Censoring limits under correction.** The correction moves each record's limit with its
    value: a censored record is known only to be below
    `pcLLOQ_ij = LLOQ * PRED_bin / PRED_ij` (and the variability correction maps that limit like
    a value). The limits now differ between records, so a censored value need not lie below
    every quantified value. An order statistic is reported only when it is the same wherever
    each censored value lies below its own limit; a quantile that depends on one that is not is
    `NaN` with `below_lloq = true`, and the observed and per-replicate rows carry its `limit`,
    the largest value its lower order statistic can take. Replicate quantiles below a limit are
    summarised across replicates with those bounds. With a common limit and no correction this is
    exactly the rule above (every limit is `lloq`). `plot_vpc` draws no LLOQ line on a
    prediction-corrected axis, since no single line is the limit there.
  - **Subject means.** With `weighting = :subject, subject_value = :mean`, a subject whose values
    in a bin include a censored one has a censored mean: it is unknown, whatever placeholder the
    censored value has. With `:censored` it is ranked as censored below its own limit (the mean
    with the censored values at their limits), with `:drop` the subject is left out of the bin,
    with `:censor_to_lloq` and `:keep` its number is used; it counts in `n_blq` and in the
    subject-level fraction below the limit in every case.
