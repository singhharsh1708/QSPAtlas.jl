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

An observation at the time of a dose or covariate change sees the state AFTER that event.
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
function steady_state_state(model::PopulationModel, prob, cycle, opts::SteadyStateOptions)
    T, ii, events = cycle
    sys = model.sys
    cb, _ = dosing_callback(sys, events, model.routes)
    p = remake(prob; tspan = (T, T + ii), callback = cb)
    kw = merge(model.solve_kwargs, (; save_on = false, save_start = false, save_end = false, save_everystep = false))
    integ = init(p, model.solver; kw...)
    getall = getu(sys, ModelingToolkit.unknowns(sys))
    x = Float64.(prob.u0)
    maxchange = Inf
    for n in 1:opts.maxcycles
        reinit!(integ, x; t0 = T, tf = T + ii)
        solve!(integ)
        successful_retcode(integ.sol.retcode) || throw(SteadyStateError(
            "no steady state reached after $n cycles: the solver stopped with retcode $(integ.sol.retcode) (ii = $ii)"))
        xn = Float64.(integ.u)
        all(isfinite, xn) || throw(SteadyStateError(
            "no steady state reached after $n cycles: the state became non-finite (the model grows without bound over a cycle, ii = $ii)"))
        if all(i -> abs(xn[i] - x[i]) <= opts.atol + opts.rtol * abs(xn[i]), eachindex(xn))
            return Float64.(getall(integ))
        end
        maxchange = maximum(abs(xn[i] - x[i]) / max(opts.atol + opts.rtol * abs(xn[i]), floatmin()) for i in eachindex(xn))
        x = xn
    end
    throw(SteadyStateError("no steady state reached after $(opts.maxcycles) cycles (ii = $ii; the largest state change over the last cycle was $(round(maxchange; sigdigits = 3)) times the tolerance atol + rtol |x| = $(opts.atol) + $(opts.rtol) |x|)"))
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

has_steady_state(events) = any(e -> e.kind === :steady_state, events)

# ---- Base problem and per-subject problems ----------------------------------------------

# Build the base ODEProblem once per model (first call wins; guarded by a lock so threads
# starting together do not build it twice).
function base_problem(model::PopulationModel, subject::Subject)
    cache = model.cache
    lock(cache.lock) do
        cache.prob === nothing || return cache.prob
        t0 = model.tspan_default[1]
        c0 = covariates_at(subject, t0)
        η0 = zero_eta(model)
        pp = collect(Pair, model.params(model.θ, η0, c0))
        check_no_route_rates(model, [first(p) for p in pp])
        ii = collect(Pair, model.init(model.θ, η0, c0))
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

# Observation recorder, always LAST in the callback set so that it sees the state after
# the dosing and covariate callbacks that fire at the same time.
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

# Parameter changes at the subject's covariate times, values computed up front.
function covariate_callback(model::PopulationModel, subject::Subject, η, tspan)
    cache = model.cache
    times = [t for t in subject.covariate_times if tspan[1] < t <= tspan[2]]
    isempty(times) && return nothing
    vals = map(times) do t
        pp = model.params(model.θ, η, covariates_at(subject, t))
        check_keys("params", [first(p) for p in pp], cache.param_keys)
        Float64[last(p) for p in pp]
    end
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

# Evaluate params/init for a subject once, in the calling thread, so that mistakes (a
# different set of parameters per subject, a dose into an unknown compartment) surface as
# ordinary errors before the possibly threaded solves start.
function check_subject(model::PopulationModel, s::Subject, η::NamedTuple, tspan)
    cache = model.cache
    c0 = covariates_at(s, tspan[1])
    check_keys("params", [first(p) for p in model.params(model.θ, η, c0)], cache.param_keys)
    check_keys("init", [first(p) for p in model.init(model.θ, η, c0)], cache.init_keys)
    for t in s.covariate_times
        tspan[1] < t <= tspan[2] || continue
        check_keys("params", [first(p) for p in model.params(model.θ, η, covariates_at(s, t))], cache.param_keys)
    end
    events = resolve(s.regimen, model.routes, s.control)
    for e in events
        e.kind === :steady_state && e.time < tspan[1] && throw(ArgumentError(
            "subject $(s.id): steady-state dose at time $(e.time) is before the start of tspan ($(tspan[1])); start the simulation at or before the steady-state dose"))
    end
    return events
end

# Problem for one subject without callbacks (its parameters and initial conditions), and the
# dose events with the steady-state vector filled in when the regimen has a steady-state dose.
# The steady state uses the subject's parameters at the steady-state time (covariates frozen at
# their value there) and starts the cycle iteration from the subject's initial conditions.
function subject_start(model::PopulationModel, base, s::Subject, η::NamedTuple, tspan, events, opts::SteadyStateOptions)
    cache = model.cache
    c0 = covariates_at(s, tspan[1])
    pp = model.params(model.θ, η, c0)
    ii = model.init(model.θ, η, c0)
    check_keys("params", [first(p) for p in pp], cache.param_keys)
    check_keys("init", [first(p) for p in ii], cache.init_keys)
    prob0 = remake(base; u0 = collect(Pair, ii), p = collect(Pair, pp), tspan)
    if has_steady_state(events)
        T = first(e.time for e in events if e.kind === :steady_state)
        if T > tspan[2]
            events = filter(e -> e.kind !== :steady_state, events)      # never reached
        else
            cache.time_dependent && !opts.allow_time_dependent && throw(ArgumentError(
                "subject $(s.id): steady state (ss = true) needs autonomous dynamics, but the model depends explicitly on time; the cycle iteration would repeat the forcing of time T in every cycle. " *
                "If the forcing is periodic with period ii and aligned with T, pass steady_state = SteadyStateOptions(allow_time_dependent = true)"))
            cT = covariates_at(s, T)
            probT = cT == c0 ? prob0 :
                remake(base; u0 = collect(Pair, ii), p = collect(Pair, model.params(model.θ, η, cT)), tspan)
            x = steady_state_state(model, probT, steady_state_cycle(s.regimen, model.routes, s.control), opts)
            events = with_steady_state(events, x)
        end
    end
    return prob0, events
end

# Problem for one subject, plus the buffer its observer fills (rows = unique obstimes,
# columns = endpoints). `prob0` comes from `subject_start`.
function subject_problem(model::PopulationModel, prob0, s::Subject, η::NamedTuple, tspan, events)
    cache = model.cache
    ut = unique(s.obstimes)
    buf = fill(NaN, length(ut), length(cache.getters))

    callbacks = Any[]
    tstops = Float64[]
    dcb, dstops = dosing_callback(model.sys, events, model.routes)
    push!(callbacks, dcb)
    append!(tstops, dstops)
    cov = covariate_callback(model, s, η, tspan)
    if cov !== nothing
        push!(callbacks, cov[1])
        append!(tstops, cov[2])
    end
    if !isempty(ut)
        push!(callbacks, observer_callback(cache.getters, ut, buf))
        append!(tstops, ut)
    end
    tstops = sort!(unique(filter(x -> tspan[1] <= x <= tspan[2], tstops)))

    # Only the final solution matters: observations are recorded by the observer callback.
    prob = remake(prob0; callback = CallbackSet(callbacks...), tstops, save_on = false, save_start = false, save_end = false)
    return prob, buf
end

# Turn a finished solve into obs rows. Residual draws (when `rng` is given) run endpoint by
# endpoint, in observation order, one draw per row.
function observation_rows(model::PopulationModel, s::Subject, replicate::Int, buf::Matrix{Float64}, retcode::Symbol, rng)
    ut = unique(s.obstimes)
    rows = NamedTuple[]
    for (e, name) in enumerate(keys(model.observed))
        err = model.error[name]
        for t in s.obstimes
            ipred = buf[searchsortedfirst(ut, t), e]
            dv = rng === nothing ? ipred : sample_dv(err, ipred, rng)
            push!(rows, (replicate = replicate, id = s.id, time = t, endpoint = name, ipred = ipred, dv = dv, retcode = retcode))
        end
    end
    return rows
end

subject_row(s::Subject, replicate::Int, η::NamedTuple, retcode::Symbol) =
    (replicate = replicate, id = s.id, η = η, covariates = s.covariates, retcode = retcode)

retcode_symbol(sol) = Symbol(sol.retcode)

# ---- Public API -------------------------------------------------------------------------

function check_eta_mode(η)
    η === :sample || η === :zero || η isa NamedTuple || η isa AbstractVector{<:Real} ||
        throw(ArgumentError("η must be :sample, :zero, a NamedTuple or a vector, got $(repr(η))"))
    return η
end
eta_mode(η) = η isa Symbol ? η : :given

function sim_settings(model::PopulationModel; seed, replicates, residual, η, steady_state)
    kw = model.solve_kwargs
    return (; seed, replicates, residual, η = eta_mode(η), steady_state,
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
             steady_state = SteadyStateOptions()) -> SimResult

Simulate one subject. `η = :sample` draws η ~ N(0, Ω) from `rng`, `:zero` uses η = 0 (the
zero-random-effect prediction), or pass a NamedTuple / vector of η values. `residual = false`
returns `dv = ipred`. The residual error uses its own stream, seeded from `rng` right after
the η draw (so η does not depend on `residual`). `tspan` defaults to the model's, extended
to the last observation time if needed.

Steady-state doses (`ss = true`, see [`Regimen`](@ref)) and resets are handled here: for a
steady-state dose at time `T` the state at `T` is found per subject with that subject's
parameters (covariates frozen at their value at `T`), by cycling from the subject's initial
conditions as set by [`SteadyStateOptions`](@ref); the solve then continues from that state with
the infusions and pending doses of the periodic history switched on. `T` must not be before
the start of `tspan` (an error); a steady-state dose after the end is ignored. If there is no
steady state, [`SteadyStateError`](@ref) is thrown by this single-subject method. (In a cohort
the same failure is a per-trajectory outcome, see below.) A reset sets every state to the subject's initial condition at the start of the solve.
"""
function simulate(model::PopulationModel, subject::Subject; η = :sample, residual::Bool = true,
        rng = Xoshiro(0), tspan = nothing, steady_state::SteadyStateOptions = SteadyStateOptions())
    check_eta_mode(η)
    ts = subject_tspan(model, subject, tspan)
    base = base_problem(model, subject)
    ηv = η === :sample ? sample_eta(model, rng) : η === :zero ? zero_eta(model) : eta_value(model, η)
    resid_rng = Xoshiro(rand(rng, UInt64))      # always taken, so rng ends in the same state
    events = check_subject(model, subject, ηv, ts)
    prob0, events = subject_start(model, base, subject, ηv, ts, events, steady_state)
    prob, buf = subject_problem(model, prob0, subject, ηv, ts, events)
    rc = retcode_symbol(solve(prob, model.solver; model.solve_kwargs...))
    rows = observation_rows(model, subject, 1, buf, rc, residual ? resid_rng : nothing)
    warn_failures(1, rc !== :Success ? 1 : 0)
    return SimResult(rows, [subject_row(subject, 1, ηv, rc)],
        sim_settings(model; seed = nothing, replicates = 1, residual, η, steady_state))
end

"""
    simulate(model, subjects; replicates = 1, seed = 0, η = :sample, residual = true,
             ensemble = EnsembleThreads(), tspan = nothing, steady_state = SteadyStateOptions()) -> SimResult

Simulate a cohort. `replicates` is the number of times the whole cohort is re-simulated
(new η and residual draws each time). Subject ids must be unique.

Reproducibility: every (replicate, subject id) has its own random streams, derived from
`seed`, one for η and a separate one for the residual error. The results therefore do not
depend on the number of threads, on `ensemble`, or on the order of `subjects` (the order
of the output rows follows the cohort order).

A trajectory whose steady-state dose has no steady state ([`SteadyStateError`](@ref)) does not
abort the cohort: it is not solved, its subject row has `retcode = :SteadyStateFailed` and an
extra field `message`, its observation rows have `ipred = dv = NaN` and the same retcode, and
the usual single warning counts it. A malformed regimen or model (steady-state dose before
`tspan`, explicit time dependence, unknown compartment, ...) still throws for the whole call.
"""
function simulate(model::PopulationModel, subjects::AbstractVector{<:Subject}; replicates::Integer = 1, seed::Integer = 0,
        η = :sample, residual::Bool = true, ensemble = EnsembleThreads(), tspan = nothing,
        steady_state::SteadyStateOptions = SteadyStateOptions())
    replicates >= 1 || throw(ArgumentError("replicates must be >= 1"))
    check_eta_mode(η)
    ids = [s.id for s in subjects]
    allunique(ids) || throw(ArgumentError("subject ids must be unique, duplicated: $(unique(filter(i -> count(==(i), ids) > 1, ids)))"))
    settings = sim_settings(model; seed, replicates, residual, η, steady_state)
    n = length(subjects)
    n == 0 && return SimResult(NamedTuple[], NamedTuple[], settings)

    spans = [subject_tspan(model, s, tspan) for s in subjects]
    base = base_problem(model, subjects[1])
    N = replicates * n
    rep_of(k) = (k - 1) ÷ n + 1
    sub_of(k) = (k - 1) % n + 1

    # η for every trajectory, drawn up front from its own stream.
    ηs = Vector{NamedTuple}(undef, N)
    given = η isa Symbol ? nothing : eta_value(model, η)
    events = Vector{Any}(undef, n)                  # dose events depend on the subject only
    for k in 1:N
        j = sub_of(k)
        s = subjects[j]
        ηs[k] = η === :sample ? sample_eta(model, stream_rng(seed, rep_of(k), s.id, :eta)) :
                η === :zero ? zero_eta(model) : given
        ev = check_subject(model, s, ηs[k], spans[j])
        k <= n && (events[j] = ev)
    end

    # Steady-state doses: the state at the start of the interval depends on the trajectory's
    # parameters, so it is found per trajectory before the main solves (in threads for
    # EnsembleThreads). A trajectory without a steady state is not solved: it is reported with
    # retcode :SteadyStateFailed and NaN observations, the others go on. Any other error is
    # rethrown here, the first in trajectory order.
    starts = Vector{Any}(undef, N)
    ss_failed = Dict{Int, String}()
    ssks = [k for k in 1:N if has_steady_state(events[sub_of(k)])]
    if !isempty(ssks)
        out = ordered_map(ensemble, ssks) do k
            j = sub_of(k)
            subject_start(model, base, subjects[j], ηs[k], spans[j], events[j], steady_state)
        end
        for (k, (ok, o)) in zip(ssks, out)
            ok ? (starts[k] = o) : (ss_failed[k] = o.msg)
        end
    end
    kruns = [k for k in 1:N if !haskey(ss_failed, k)]

    bufs = Vector{Matrix{Float64}}(undef, N)       # slot k is written by trajectory k only
    prob_func = function (prob, ctx)
        k = kruns[ctx.sim_id]
        j = sub_of(k)
        prob0, evs = isassigned(starts, k) ? starts[k] : subject_start(model, prob, subjects[j], ηs[k], spans[j], events[j], steady_state)
        p, bufs[k] = subject_problem(model, prob0, subjects[j], ηs[k], spans[j], evs)
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
            nan_buf = fill(NaN, length(unique(s.obstimes)), length(model.cache.getters))
            rng = residual ? stream_rng(seed, rep_of(k), s.id, :residual) : nothing
            append!(obs, observation_rows(model, s, rep_of(k), nan_buf, :SteadyStateFailed, rng))
            push!(subs, merge(subject_row(s, rep_of(k), ηs[k], :SteadyStateFailed), (; message = ss_failed[k])))
        else
            append!(obs, results[k].rows)
            push!(subs, subject_row(s, rep_of(k), ηs[k], results[k].retcode))
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

Prediction with all random effects at zero and no residual error. This is the
typical-individual curve for the subject's covariates; it is generally not the
population mean.
"""
zero_re_prediction(model::PopulationModel, subject::Subject; kwargs...) = simulate(model, subject; η = :zero, residual = false, kwargs...)
