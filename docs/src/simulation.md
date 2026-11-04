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
holds subject-specific bioavailability and lag time.

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

Two differences from NONMEM are deliberate. The iteration starts from the model's initial
conditions rather than from zero (mrgsolve starts from the state at the steady-state
record), which matters for models with endogenous baselines (cell
counts, receptors) and gives the same answer for linear PK. Steady state is refused for
models whose equations depend on time explicitly (for example a resistance term
`exp(-δ t)`), because a periodic state does not exist there; `SteadyStateOptions(allow_time_dependent = true)`
overrides this when the time dependence is known to be periodic with `ii`. If no steady
state is reached, a single-subject simulation throws `SteadyStateError`; in a cohort the
subject is kept with `retcode = :SteadyStateFailed` and empty observations, and the rest
of the cohort is simulated.

### Resets

`Reset(time)` returns every state to its initial conditions and cancels running infusions
and lagged doses that have not yet been given (NONMEM EVID 3); later `addl` doses are kept.
A reset followed by a dose at the same time (EVID 4) applies the reset first, then the dose
with its own lag.

### Observations at event times

An observation at the same time as a dose, an infusion change, a reset or a covariate
change sees the state after the event (mrgsolve `recsort = 3`).

### Validation

The dosing semantics above were checked case by case against mrgsolve; see
[Validation](validation.md).

## Subjects and population models

```julia
sub = Subject("1"; regimen, covariates = (WT = 70.0,), obstimes = collect(0:1.0:48),
              timevarying = [24.0 => (WT = 60.0,)])
```

Covariates are constant or piecewise constant; a change at time `t` re-evaluates the
model parameters from then on. Infusion rate parameters belong to the dosing code and may
not be set by covariates.

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
