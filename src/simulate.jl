# Simulation: simulate, zero_re_prediction, SimResult. Semantics: docs/src/simulation.md.
import ModelingToolkit            # symbolic ODEProblem(sys, pairs, tspan) and remake with pairs
using DiffEqCallbacks: PresetTimeCallback
using Random: Xoshiro, rand
using SciMLBase: CallbackSet, EnsembleProblem, EnsembleThreads, EnsembleSerial, ODEProblem, ReturnCode, remake, solve, derivative_discontinuity!,
    init, reinit!, solve!, successful_retcode
using SymbolicIndexingInterface: getu, setp

"""
    SimResult

Output of [`simulate`](@ref), as tidy rows.

  - `obs`: one row per (replicate, subject, endpoint, observation time):
    `(replicate, id, time, endpoint, ipred, dv, retcode)`. `ipred` is the individual
    prediction, `dv` the observation with residual error (equal to `ipred` when
    `residual = false`). If the solver stopped early (`retcode` is not `:Success`), the
    observations it never reached have `ipred = dv = NaN`; failed subjects are never dropped.
  - `subjects`: one row per (replicate, subject): `(replicate, id, η, covariates, retcode)`.
    A cohort trajectory without a steady state has `retcode = :SteadyStateFailed` and an extra field `message`.
  - `settings`: seed, solver and tolerances, package versions.

An observation at the time of a dose or covariate change sees the state AFTER that event,
unless the subject tags it `:pre` (see [`Subject`](@ref)); rows follow the subject's observation
order (time, then `:pre` before `:post`), endpoint by endpoint, and an observation restricted to
one endpoint gives a row for that endpoint only.
"""
struct SimResult
    obs::Vector{NamedTuple}
    subjects::Vector{NamedTuple}
    settings::NamedTuple
end

# NaN-safe equality, so two identical runs that contain failed subjects compare equal.
Base.:(==)(a::SimResult, b::SimResult) = isequal(a.obs, b.obs) && isequal(a.subjects, b.subjects) && isequal(a.settings, b.settings)
Base.isequal(a::SimResult, b::SimResult) = a == b
Base.hash(a::SimResult, h::UInt) = hash(a.settings, hash(a.subjects, hash(a.obs, h)))
Base.length(r::SimResult) = length(r.obs)
function Base.show(io::IO, r::SimResult)
    nsub = length(r.subjects)
    nrep = length(unique(s.replicate for s in r.subjects))
    nbad = count(s -> s.retcode !== :Success, r.subjects)
    print(io, "SimResult(", length(r.obs), " observations, ", nsub, " subject runs in ", nrep, " replicate(s)",
        nbad == 0 ? "" : ", $nbad FAILED", ")")
end

# ---- Steady state -----------------------------------------------------------------------

"""
    SteadyStateOptions(; rtol = 1e-8, atol = 1e-10, maxcycles = 1000, allow_time_dependent = false)

Convergence control for steady-state doses (`ss = true`). The state at the start of a dosing
interval is found by simulating cycle after cycle (length `ii`, same model, parameters, F and
lag) from the initial conditions. It has converged when, for EVERY state, the change over a
cycle is at most `atol + rtol * |x|` (compared directly, so a state that stays exactly 0 is
fine even with `atol = 0`). After `maxcycles` cycles without convergence, or if the
state becomes non-finite or the solver fails in a cycle, [`SteadyStateError`](@ref) is thrown;
a non-converged state is never used. The change over one cycle underestimates the distance to
the steady state by about `1 / (1 - r)` where `r` is the per-cycle contraction factor (a state
that decays slowly compared to `ii` needs a tighter `rtol` or a higher `maxcycles`).

The cycle iteration restarts every cycle at the same absolute time `T`, so it is only valid
for AUTONOMOUS dynamics (plus the dosing). If the compiled system depends explicitly on time
(an `exp(-δ t)` or `sin(t)` in an equation or observed expression), a steady-state dose is
refused with an `ArgumentError`. `allow_time_dependent = true` overrides that: the user asserts
that the forcing is periodic with period `ii` (or a divisor of it) and phase-aligned with `T`.
"""
struct SteadyStateOptions
    rtol::Float64
    atol::Float64
    maxcycles::Int
    allow_time_dependent::Bool
    function SteadyStateOptions(; rtol::Real = 1e-8, atol::Real = 1e-10, maxcycles::Integer = 1000, allow_time_dependent::Bool = false)
        (isfinite(rtol) && rtol >= 0) || throw(ArgumentError("rtol must be finite and >= 0, got $rtol"))
        (isfinite(atol) && atol >= 0) || throw(ArgumentError("atol must be finite and >= 0, got $atol"))
        rtol > 0 || atol > 0 || throw(ArgumentError("rtol and atol cannot both be 0"))
        maxcycles >= 1 || throw(ArgumentError("maxcycles must be >= 1, got $maxcycles"))
        return new(rtol, atol, maxcycles, allow_time_dependent)
    end
end

"""
    SteadyStateError(msg)

Thrown when a steady-state dose has no steady state within the [`SteadyStateOptions`](@ref)
limits. The message starts with "no steady state reached after N cycles".
"""
struct SteadyStateError <: Exception
    msg::String
end
Base.showerror(io::IO, e::SteadyStateError) = print(io, "SteadyStateError: ", e.msg)

# Does `expr` mention the independent variable `iv` other than as the argument of a state
# (x(t)) or inside a derivative?
function mentions_time(expr, iv)
    SU = ModelingToolkit.SymbolicUtils
    ex = ModelingToolkit.unwrap(expr)
    isequal(ex, iv) && return true
    SU.iscall(ex) || return false
    op = SU.operation(ex)
    op isa ModelingToolkit.Differential && return false
    args = SU.arguments(ex)
    if !(op isa Function) && length(args) == 1 && isequal(args[1], iv)
        return false                       # a state x(t)
    end
    return any(a -> mentions_time(a, iv), args)
end

"""
    explicit_time_dependence(sys) -> Bool

True if an equation or observed expression of the compiled system `sys` depends on time
explicitly, beyond the states.
"""
function explicit_time_dependence(sys)
    iv = ModelingToolkit.get_iv(sys)
    eqs = vcat(ModelingToolkit.equations(sys), ModelingToolkit.observed(sys))
    return any(e -> mentions_time(e.lhs, iv) || mentions_time(e.rhs, iv), eqs)
end

# State (in the order of the system's unknowns) at the start of the cycle [T, T + ii) of the
# periodic regimen, as a fixed point of the cycle map found by plain iteration from the initial
# state of `prob`. `prob` carries the subject's parameters and initial conditions.
#
# A constant infusion (ii == 0) has no period: its steady state is the equilibrium under the
# infusion. Then the cycle length starts at 1 time unit and doubles every cycle (capped at 2^40),
# with the same convergence test, so slow and fast dynamics both converge in few cycles.
function steady_state_state(model::PopulationModel, prob, cycle, opts::SteadyStateOptions)
    T, ii, events = cycle
    constant = ii == 0
    len(n) = constant ? ldexp(1.0, min(n - 1, 40)) : ii
    sys = model.sys
    cb, _ = dosing_callback(sys, events, model.routes)
    p = remake(prob; tspan = (T, T + len(1)), callback = cb)
    kw = merge(model.solve_kwargs, (; save_on = false, save_start = false, save_end = false, save_everystep = false))
    integ = init(p, model.solver; kw...)
    getall = getu(sys, ModelingToolkit.unknowns(sys))
    x = Float64.(prob.u0)
    maxchange = Inf
    for n in 1:opts.maxcycles
        reinit!(integ, x; t0 = T, tf = T + len(n))
        solve!(integ)
        successful_retcode(integ.sol.retcode) || throw(SteadyStateError(
            "no steady state reached after $n cycles: the solver stopped with retcode $(integ.sol.retcode) ($(constant ? "constant infusion" : "ii = $ii"))"))
        xn = Float64.(integ.u)
        all(isfinite, xn) || throw(SteadyStateError(
            "no steady state reached after $n cycles: the state became non-finite (the model grows without bound over a cycle, $(constant ? "constant infusion" : "ii = $ii"))"))
        if all(i -> abs(xn[i] - x[i]) <= opts.atol + opts.rtol * abs(xn[i]), eachindex(xn))
            return Float64.(getall(integ))
        end
        maxchange = maximum(abs(xn[i] - x[i]) / max(opts.atol + opts.rtol * abs(xn[i]), floatmin()) for i in eachindex(xn))
        x = xn
    end
    throw(SteadyStateError("no steady state reached after $(opts.maxcycles) cycles ($(constant ? "constant infusion" : "ii = $ii"); the largest state change over the last cycle was $(round(maxchange; sigdigits = 3)) times the tolerance atol + rtol |x| = $(opts.atol) + $(opts.rtol) |x|)"))
end

# `f(k)` for every k, as `(true, value)` or `(false, exception)` in order; in threads for
# EnsembleThreads. Errors are captured per task. Any error except `SteadyStateError` (a
# malformed regimen or model, which affects everyone) is rethrown in the calling thread, the
# first one in k order, as itself and not as a wrapped task failure.
function ordered_map(f, ensemble, ks)
    safe(k) = try
        (true, f(k))
    catch e
        (false, e)
    end
    res = ensemble isa EnsembleThreads ? fetch.([Threads.@spawn(safe(k)) for k in ks]) : [safe(k) for k in ks]
    for r in res
        r[1] || r[2] isa SteadyStateError || throw(r[2])      # only a missing steady state is a trajectory outcome
    end
    return res
end

has_steady_state(events) = any(e -> e.kind === :steady_state || e.kind === :steady_state_add, events)

# ---- Base problem and per-subject problems ----------------------------------------------

# Build the base ODEProblem once per model (first call wins; guarded by a lock so threads
# starting together do not build it twice).
function base_problem(model::PopulationModel, subject::Subject)
    cache = model.cache
    lock(cache.lock) do
        cache.prob === nothing || return cache.prob
        t0 = model.tspan_default[1]
        c0 = covariates_at(subject, t0)
        η0, κ0 = zero_eta(model), zero_kappa(model)
        pp = collect(Pair, call_model(model.params, model, η0, c0, κ0))
        check_no_route_rates(model, [first(p) for p in pp])
        ii = collect(Pair, call_model(model.init, model, η0, c0, κ0))
        prob = ODEProblem(model.sys, vcat(ii, pp), model.tspan_default)
        cache.param_keys = Any[first(p) for p in pp]
        cache.init_keys = Any[first(p) for p in ii]
        cache.param_setter = isempty(pp) ? nothing : setp(model.sys, cache.param_keys)
        cache.getters = Any[getu(model.sys, o) for o in values(model.observed)]
        cache.time_dependent = explicit_time_dependence(model.sys)
        cache.prob = prob
        return prob
    end
end

function check_keys(what, got, want)
    ok = length(got) == length(want) && all(isequal(a, b) for (a, b) in zip(got, want))
    ok || throw(ArgumentError("`$what` must return the same $(what == "params" ? "parameters" : "states") in the same order for every subject"))
    return nothing
end

# Observation recorder. The `:post` recorder comes after the dosing callback in the callback
# set, so it sees the state after the events at its times; the `:pre` recorder comes before it.
function observer_callback(getters, ut::Vector{Float64}, buf::Matrix{Float64})
    function record!(integrator)
        j = searchsortedfirst(ut, integrator.t)
        if j <= length(ut) && ut[j] == integrator.t
            for (e, g) in enumerate(getters)
                buf[j, e] = g(integrator)
            end
        end
        return nothing
    end
    return PresetTimeCallback(ut, record!; save_positions = (false, false))
end

# Parameter changes at the subject's covariate times (and occasion changes, for a model with
# IOV), values computed up front. Under :locf a change at t applies from t on (it fires before the
# observers of t); under :nocb the values of the next interval start just after t (it fires after
# the observers of t). Returns nothing or (callback, times).
function covariate_callback(model::PopulationModel, subject::Subject, η, κs, tspan)
    cache = model.cache
    times, vals = covariate_switches(model, subject, η, κs, tspan)
    isempty(times) && return nothing
    setter = cache.param_setter
    function change!(integrator)
        i = findfirst(==(integrator.t), times)
        i === nothing && return nothing
        setter(integrator, vals[i])
        derivative_discontinuity!(integrator, true)     # the right-hand side changed
        return nothing
    end
    return PresetTimeCallback(times, change!; save_positions = (true, true)), times
end

# Times at which the parameters change, and their new values (same order as cache.param_keys).
function covariate_switches(model::PopulationModel, s::Subject, η, κs, tspan)
    nocb = s.covariate_rule === :nocb
    cand = has_iov(model) ? sort!(unique(vcat(s.covariate_times, s.occasion_times))) : s.covariate_times
    times = [t for t in cand if nocb ? tspan[1] <= t < tspan[2] : tspan[1] < t <= tspan[2]]
    vals = map(times) do t
        c, o = nocb ? (covariates_after(s, t), occasion_after(s, t)) : (covariates_at(s, t), occasion_at(s, t))
        pp = call_model(model.params, model, η, c, kappa_of(κs, o))
        check_keys("params", [first(p) for p in pp], model.cache.param_keys)
        Float64[last(p) for p in pp]
    end
    return times, vals
end

# Distinct observation times of a subject, separately for `:pre` and `:post` observations.
obs_slots(s::Subject) = (pre = unique(t for (t, o) in zip(s.obstimes, s.obsorder) if o === :pre),
    post = unique(t for (t, o) in zip(s.obstimes, s.obsorder) if o === :post))
obs_buffers(slots, nend) = (pre = fill(NaN, length(slots.pre), nend), post = fill(NaN, length(slots.post), nend))

function subject_tspan(model::PopulationModel, s::Subject, tspan)
    t0, t1 = tspan === nothing ? model.tspan_default : (Float64(tspan[1]), Float64(tspan[2]))
    if tspan === nothing && !isempty(s.obstimes)
        t1 = max(t1, s.obstimes[end])
    end
    t0 < t1 || throw(ArgumentError("tspan must have t0 < tend, got ($t0, $t1)"))
    if !isempty(s.obstimes) && (s.obstimes[1] < t0 || s.obstimes[end] > t1)
        throw(ArgumentError("subject $(s.id): observation times [$(s.obstimes[1]), $(s.obstimes[end])] are outside tspan ($t0, $t1)"))
    end
    return (t0, t1)
end

# κ of occasion `o` (an empty NamedTuple for a model without IOV).
kappa_of(κs, o) = isempty(κs) ? (;) : κs[o]

"""
    trajectory_control(model, subject, η, κs) -> DoseControl or function

The dose control of one trajectory: the subject's `control` when the model has no
`dosecontrol`, else a function `(dose, record_time) -> DoseControl` that evaluates the model's
`dosecontrol` with the covariates and the occasion's κ at the time of the dose's record
(`covariates_at`, `occasion_at`; `addl` doses use their record's time) and combines it with the
subject's `control`.
"""
function trajectory_control(model::PopulationModel, s::Subject, η::NamedTuple, κs)
    model.dosecontrol === nothing && return s.control
    return function (d::Dose, rt::Real)
        dc = as_dosecontrol(call_model(model.dosecontrol, model, η, covariates_at(s, rt), kappa_of(κs, occasion_at(s, rt))))
        return merge_control(s.control, dc)
    end
end

# Evaluate params/init for a subject once, in the calling thread, so that mistakes (a
# different set of parameters per subject, a dose into an unknown compartment) surface as
# ordinary errors before the possibly threaded solves start.
function check_subject(model::PopulationModel, s::Subject, η::NamedTuple, κs, tspan)
    cache = model.cache
    c0 = covariates_at(s, tspan[1])
    κ0 = kappa_of(κs, occasion_at(s, tspan[1]))
    check_keys("params", [first(p) for p in call_model(model.params, model, η, c0, κ0)], cache.param_keys)
    check_keys("init", [first(p) for p in call_model(model.init, model, η, c0, κ0)], cache.init_keys)
    covariate_switches(model, s, η, κs, tspan)
    for ep in s.obsendpoints
        ep === nothing || haskey(model.observed, ep) || throw(ArgumentError(
            "subject $(s.id): observation for endpoint :$ep, but the model's endpoints are $(collect(keys(model.observed)))"))
    end
    events = resolve(s.regimen, model.routes, trajectory_control(model, s, η, κs))
    for e in events
        (e.kind === :steady_state || e.kind === :steady_state_add) && e.time < tspan[1] && throw(ArgumentError(
            "subject $(s.id): steady-state dose at time $(e.time) is before the start of tspan ($(tspan[1])); start the simulation at or before the steady-state dose"))
    end
    return events
end

# Problem for one subject without callbacks (its parameters and initial conditions), and the
# dose events with the steady-state vector filled in when the regimen has a steady-state dose.
# The steady state uses the subject's parameters at the steady-state time (covariates and
# occasion frozen at their value there) and starts the cycle iteration from the subject's initial
# conditions.
function subject_start(model::PopulationModel, base, s::Subject, η::NamedTuple, κs, tspan, events, opts::SteadyStateOptions)
    cache = model.cache
    c0 = covariates_at(s, tspan[1])
    κ0 = kappa_of(κs, occasion_at(s, tspan[1]))
    pp = call_model(model.params, model, η, c0, κ0)
    ii = call_model(model.init, model, η, c0, κ0)
    check_keys("params", [first(p) for p in pp], cache.param_keys)
    check_keys("init", [first(p) for p in ii], cache.init_keys)
    prob0 = remake(base; u0 = collect(Pair, ii), p = collect(Pair, pp), tspan)
    if has_steady_state(events)
        blocks = [b for b in ss_blocks(s.regimen) if b.time <= tspan[2]]     # later ones are never reached
        events = filter(e -> !((e.kind === :steady_state || e.kind === :steady_state_add) && e.time > tspan[2]), events)
        !isempty(blocks) && cache.time_dependent && !opts.allow_time_dependent && throw(ArgumentError(
            "subject $(s.id): steady state (ss = 1 or 2) needs autonomous dynamics, but the model depends explicitly on time; the cycle iteration would repeat the forcing of time T in every cycle. " *
            "If the forcing is periodic with period ii and aligned with T, pass steady_state = SteadyStateOptions(allow_time_dependent = true)"))
        rmap = _route_map(model.routes)
        ctl = trajectory_control(model, s, η, κs)
        set = nothing
        adds = Vector{Float64}[]
        for b in blocks
            cT = covariates_at(s, b.time)
            κT = kappa_of(κs, occasion_at(s, b.time))
            probT = cT == c0 && κT == κ0 ? prob0 :
                remake(base; u0 = collect(Pair, ii), p = collect(Pair, call_model(model.params, model, η, cT, κT)), tspan)
            cyc = block_cycle(b, rmap, ctl)
            x = steady_state_state(model, probT, cyc, opts)
            if b.mode == 1
                set = x
            else
                # SS=2: the part of the periodic state that the dose adds
                push!(adds, x .- steady_state_state(model, probT, (cyc[1], cyc[2], DoseEvent[]), opts))
            end
        end
        events = fill_steady_states(events, set, adds)
    end
    return prob0, events
end

# Problem for one subject, plus the buffer its observer fills (rows = unique obstimes,
# columns = endpoints). `prob0` comes from `subject_start`.
function subject_problem(model::PopulationModel, prob0, s::Subject, η::NamedTuple, κs, tspan, events)
    cache = model.cache
    slots = obs_slots(s)
    buf = obs_buffers(slots, length(cache.getters))

    # Order matters for events at the same time: LOCF covariate changes, `:pre` observations,
    # doses / resets / steady states, `:post` observations, NOCB covariate changes.
    callbacks = Any[]
    tstops = Float64[]
    cov = covariate_callback(model, s, η, κs, tspan)
    if cov !== nothing && s.covariate_rule === :locf
        push!(callbacks, cov[1])
        append!(tstops, cov[2])
    end
    if !isempty(slots.pre)
        push!(callbacks, observer_callback(cache.getters, slots.pre, buf.pre))
        append!(tstops, slots.pre)
    end
    dcb, dstops = dosing_callback(model.sys, events, model.routes)
    push!(callbacks, dcb)
    append!(tstops, dstops)
    if !isempty(slots.post)
        push!(callbacks, observer_callback(cache.getters, slots.post, buf.post))
        append!(tstops, slots.post)
    end
    if cov !== nothing && s.covariate_rule === :nocb
        push!(callbacks, cov[1])
        append!(tstops, cov[2])
    end
    tstops = sort!(unique(filter(x -> tspan[1] <= x <= tspan[2], tstops)))

    # Only the final solution matters: observations are recorded by the observer callbacks.
    prob = remake(prob0; callback = CallbackSet(callbacks...), tstops, save_on = false, save_start = false, save_end = false)
    return prob, buf
end

# Turn a finished solve into obs rows. Residual draws (when `rng` is given) run endpoint by
# endpoint, in observation order, one draw per row.
function observation_rows(model::PopulationModel, s::Subject, replicate::Int, buf::NamedTuple, retcode::Symbol, rng)
    slots = obs_slots(s)
    rows = NamedTuple[]
    for (e, name) in enumerate(keys(model.observed))
        err = model.error[name]
        for (j, t) in enumerate(s.obstimes)
            ep = s.obsendpoints[j]
            (ep === nothing || ep === name) || continue
            pre = s.obsorder[j] === :pre
            ut, b = pre ? (slots.pre, buf.pre) : (slots.post, buf.post)
            ipred = b[searchsortedfirst(ut, t), e]
            dv = rng === nothing ? ipred : sample_dv(err, ipred, rng)
            push!(rows, (replicate = replicate, id = s.id, time = t, endpoint = name, ipred = ipred, dv = dv, retcode = retcode))
        end
    end
    return rows
end

subject_row(s::Subject, replicate::Int, η::NamedTuple, retcode::Symbol) =
    (replicate = replicate, id = s.id, η = η, covariates = s.covariates, retcode = retcode)
# With IOV the row also has `κ`, a Dict occasion label => κ, placed after η.
subject_row(s::Subject, replicate::Int, η::NamedTuple, κs, retcode::Symbol) =
    isempty(κs) ? subject_row(s, replicate, η, retcode) :
    (replicate = replicate, id = s.id, η = η, κ = κs, covariates = s.covariates, retcode = retcode)

retcode_symbol(sol) = Symbol(sol.retcode)

# ---- Public API -------------------------------------------------------------------------

function check_eta_mode(η)
    η === :sample || η === :zero || η isa NamedTuple || η isa AbstractVector{<:Real} ||
        throw(ArgumentError("η must be :sample, :zero, a NamedTuple or a vector, got $(repr(η))"))
    return η
end
eta_mode(η) = η isa Symbol ? η : :given

# κ mode: `nothing` (the default) is :zero when η is :zero, else :sample.
function kappa_mode(model::PopulationModel, κ, η)
    κ === nothing && return η === :zero ? :zero : :sample
    has_iov(model) || κ === :zero || κ === :sample ||
        throw(ArgumentError("κ values given, but the model has no between-occasion random effects (κnames)"))
    κ === :sample || κ === :zero || κ isa NamedTuple || κ isa AbstractVector{<:Real} || κ isa AbstractDict ||
        throw(ArgumentError("κ must be :sample, :zero, a NamedTuple or vector (every occasion), or a Dict occasion => value, got $(repr(κ))"))
    return κ
end

function sim_settings(model::PopulationModel; seed, replicates, residual, η, steady_state, κ = :sample)
    kw = model.solve_kwargs
    base = (; seed, replicates, residual, η = eta_mode(η), steady_state)
    has_iov(model) && (base = merge(base, (κ = κ isa Symbol ? κ : :given,)))
    return (; base...,
        solver = string(nameof(typeof(model.solver))),
        abstol = get(kw, :abstol, nothing), reltol = get(kw, :reltol, nothing),
        solve_kwargs = kw, tspan_default = model.tspan_default,
        versions = package_versions(model.solver))
end

function package_versions(solver)
    v(m) = (x = try pkgversion(m) catch; nothing end; x === nothing ? "unknown" : string(x))
    return (julia = string(VERSION), ModelingToolkit = v(ModelingToolkit), SciMLBase = v(parentmodule(ODEProblem)),
        solver_package = string(nameof(parentmodule(typeof(solver)))), solver_package_version = v(parentmodule(typeof(solver))))
end

"""
    simulate(model, subject; η = :sample, residual = true, rng = Xoshiro(0), tspan = nothing,
             steady_state = SteadyStateOptions(), κ = nothing) -> SimResult

Simulate one subject. `η = :sample` draws η ~ N(0, Ω) from `rng`, `:zero` uses η = 0 (the
zero-random-effect prediction), or pass a NamedTuple / vector of η values. `residual = false`
returns `dv = ipred`. The residual error uses its own stream, seeded from `rng` right after
the η draw (so η does not depend on `residual`). `tspan` defaults to the model's, extended
to the last observation time if needed.

Steady-state doses (`ss = true` or `ss = 2`, see [`Regimen`](@ref)) and resets are handled here: for a
steady-state dose at time `T` the state at `T` is found per subject with that subject's
parameters (covariates frozen at their value at `T`), by cycling from the subject's initial
conditions as set by [`SteadyStateOptions`](@ref); the solve then continues from that state with
the infusions and pending doses of the periodic history switched on. The steady state uses the
parameters, F and lag of the occasion in force at `T`. `T` must not be before
the start of `tspan` (an error); a steady-state dose after the end is ignored. If there is no
steady state, [`SteadyStateError`](@ref) is thrown by this single-subject method. (In a cohort
the same failure is a per-trajectory outcome, see below.) A reset sets every state to the subject's initial condition at the start of the solve.

Between-occasion random effects (a model with `κnames`): `κ = nothing` (the default) samples
them, or sets them to zero when `η = :zero`; `:sample`, `:zero`, one NamedTuple / vector for
every occasion, or a `Dict(occasion label => κ)` can be given. Sampled κ of occasion `j` come
from their own stream, seeded by one draw from `rng` taken after the residual seed (only for a
model with IOV, so `rng` is used as before otherwise), and keyed by the occasion label, so they
do not depend on the other occasions. The subject row then has a field `κ` (label => κ).
"""
function simulate(model::PopulationModel, subject::Subject; η = :sample, residual::Bool = true,
        rng = Xoshiro(0), tspan = nothing, steady_state::SteadyStateOptions = SteadyStateOptions(), κ = nothing)
    check_eta_mode(η)
    κm = kappa_mode(model, κ, η)
    ts = subject_tspan(model, subject, tspan)
    base = base_problem(model, subject)
    ηv = η === :sample ? sample_eta(model, rng) : η === :zero ? zero_eta(model) : eta_value(model, η)
    resid_rng = Xoshiro(rand(rng, UInt64))      # always taken, so rng ends in the same state
    κseed = has_iov(model) ? rand(rng, Int64) : 0
    κs = subject_kappas(model, subject, κm, κseed, 1)
    events = check_subject(model, subject, ηv, κs, ts)
    prob0, events = subject_start(model, base, subject, ηv, κs, ts, events, steady_state)
    prob, buf = subject_problem(model, prob0, subject, ηv, κs, ts, events)
    rc = retcode_symbol(solve(prob, model.solver; model.solve_kwargs...))
    rows = observation_rows(model, subject, 1, buf, rc, residual ? resid_rng : nothing)
    warn_failures(1, rc !== :Success ? 1 : 0)
    return SimResult(rows, [subject_row(subject, 1, ηv, κs, rc)],
        sim_settings(model; seed = nothing, replicates = 1, residual, η, steady_state, κ = κm))
end

"""
    simulate(model, subjects; replicates = 1, seed = 0, η = :sample, residual = true,
             ensemble = EnsembleThreads(), tspan = nothing, steady_state = SteadyStateOptions(),
             κ = nothing) -> SimResult

Simulate a cohort. `replicates` is the number of times the whole cohort is re-simulated
(new η and residual draws each time). Subject ids must be unique.

Reproducibility: every (replicate, subject id) has its own random streams, derived from
`seed`, one for η and a separate one for the residual error. The results therefore do not
depend on the number of threads, on `ensemble`, or on the order of `subjects` (the order
of the output rows follows the cohort order). Between-occasion random effects (`κ`, as for one
subject) are drawn from a third stream per (replicate, subject id, occasion label), so they do
not depend on the other occasions either.

A trajectory whose steady-state dose has no steady state ([`SteadyStateError`](@ref)) does not
abort the cohort: it is not solved, its subject row has `retcode = :SteadyStateFailed` and an
extra field `message`, its observation rows have `ipred = dv = NaN` and the same retcode, and
the usual single warning counts it. A malformed regimen or model (steady-state dose before
`tspan`, explicit time dependence, unknown compartment, ...) still throws for the whole call.
"""
function simulate(model::PopulationModel, subjects::AbstractVector{<:Subject}; replicates::Integer = 1, seed::Integer = 0,
        η = :sample, residual::Bool = true, ensemble = EnsembleThreads(), tspan = nothing,
        steady_state::SteadyStateOptions = SteadyStateOptions(), κ = nothing)
    replicates >= 1 || throw(ArgumentError("replicates must be >= 1"))
    check_eta_mode(η)
    κm = kappa_mode(model, κ, η)
    ids = [s.id for s in subjects]
    allunique(ids) || throw(ArgumentError("subject ids must be unique, duplicated: $(unique(filter(i -> count(==(i), ids) > 1, ids)))"))
    settings = sim_settings(model; seed, replicates, residual, η, steady_state, κ = κm)
    n = length(subjects)
    n == 0 && return SimResult(NamedTuple[], NamedTuple[], settings)

    spans = [subject_tspan(model, s, tspan) for s in subjects]
    base = base_problem(model, subjects[1])
    N = replicates * n
    rep_of(k) = (k - 1) ÷ n + 1
    sub_of(k) = (k - 1) % n + 1

    # η and κ for every trajectory, drawn up front from their own streams.
    ηs = Vector{NamedTuple}(undef, N)
    κss = Vector{Dict{Int, NamedTuple}}(undef, N)
    given = η isa Symbol ? nothing : eta_value(model, η)
    events = Vector{Any}(undef, N)                  # dose events (F, lag may depend on η and κ)
    for k in 1:N
        j = sub_of(k)
        s = subjects[j]
        ηs[k] = η === :sample ? sample_eta(model, stream_rng(seed, rep_of(k), s.id, :eta)) :
                η === :zero ? zero_eta(model) : given
        κss[k] = subject_kappas(model, s, κm, seed, rep_of(k))
        events[k] = check_subject(model, s, ηs[k], κss[k], spans[j])
    end

    # Steady-state doses: the state at the start of the interval depends on the trajectory's
    # parameters, so it is found per trajectory before the main solves (in threads for
    # EnsembleThreads). A trajectory without a steady state is not solved: it is reported with
    # retcode :SteadyStateFailed and NaN observations, the others go on. Any other error is
    # rethrown here, the first in trajectory order.
    starts = Vector{Any}(undef, N)
    ss_failed = Dict{Int, String}()
    ssks = [k for k in 1:N if has_steady_state(events[k])]
    if !isempty(ssks)
        out = ordered_map(ensemble, ssks) do k
            j = sub_of(k)
            subject_start(model, base, subjects[j], ηs[k], κss[k], spans[j], events[k], steady_state)
        end
        for (k, (ok, o)) in zip(ssks, out)
            ok ? (starts[k] = o) : (ss_failed[k] = o.msg)
        end
    end
    kruns = [k for k in 1:N if !haskey(ss_failed, k)]

    bufs = Vector{Any}(undef, N)                   # slot k is written by trajectory k only
    prob_func = function (prob, ctx)
        k = kruns[ctx.sim_id]
        j = sub_of(k)
        prob0, evs = isassigned(starts, k) ? starts[k] : subject_start(model, prob, subjects[j], ηs[k], κss[k], spans[j], events[k], steady_state)
        p, bufs[k] = subject_problem(model, prob0, subjects[j], ηs[k], κss[k], spans[j], evs)
        return p
    end
    output_func = function (sol, ctx)
        k = kruns[ctx.sim_id]
        s = subjects[sub_of(k)]
        rc = retcode_symbol(sol)
        rng = residual ? stream_rng(seed, rep_of(k), s.id, :residual) : nothing
        return (k = k, retcode = rc, rows = observation_rows(model, s, rep_of(k), bufs[k], rc, rng)), false
    end
    results = Vector{Any}(undef, N)
    if !isempty(kruns)
        ens = EnsembleProblem(base; prob_func, output_func, safetycopy = false)
        out = solve(ens, model.solver, ensemble; trajectories = length(kruns), model.solve_kwargs...)
        for o in out.u
            results[o.k] = o
        end
    end
    obs = NamedTuple[]
    subs = NamedTuple[]
    for k in 1:N
        s = subjects[sub_of(k)]
        if haskey(ss_failed, k)
            # no steady state: NaN observations (residual draws keep their stream positions)
            nan_buf = obs_buffers(obs_slots(s), length(model.cache.getters))
            rng = residual ? stream_rng(seed, rep_of(k), s.id, :residual) : nothing
            append!(obs, observation_rows(model, s, rep_of(k), nan_buf, :SteadyStateFailed, rng))
            push!(subs, merge(subject_row(s, rep_of(k), ηs[k], κss[k], :SteadyStateFailed), (; message = ss_failed[k])))
        else
            append!(obs, results[k].rows)
            push!(subs, subject_row(s, rep_of(k), ηs[k], κss[k], results[k].retcode))
        end
    end
    # Order rows by replicate, then cohort order (k already is that order).
    warn_failures(N, count(s -> s.retcode !== :Success, subs))
    return SimResult(obs, subs, settings)
end

function warn_failures(total, failed)
    failed > 0 && @warn "$failed of $total subject simulations did not finish with retcode Success; their unreached observations are NaN (see the retcode column)"
    return nothing
end

"""
    zero_re_prediction(model, subject) -> SimResult

Prediction with all random effects (η and, with IOV, κ) at zero and no residual error. This is the
typical-individual curve for the subject's covariates; it is generally not the
population mean.
"""
zero_re_prediction(model::PopulationModel, subject::Subject; kwargs...) = simulate(model, subject; η = :zero, κ = :zero, residual = false, kwargs...)
