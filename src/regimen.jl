# Dosing regimens: what is given (Regimen), how it reaches the model (Route),
# subject-specific F and lag (DoseControl), and the solver callback (dosing_callback).
# Semantics (event ordering, F and lag, infusions, steady state, resets): docs/src/simulation.md.

import ModelingToolkit            # unknowns(sys), for resets and steady-state jumps
using SciMLBase: derivative_discontinuity!
using SymbolicIndexingInterface: getp, setp, getu, setu

abstract type AbstractDoseTarget end

# ---------------------------------------------------------------------------------------------
# Route
# ---------------------------------------------------------------------------------------------

"""
    Route(cmt, state; rate = nothing, scale = 1.0)

Map a dose compartment name `cmt` to the model. `state` is the unknown of the compiled
system that receives bolus amounts (for example `sys.A1`). `rate` is the declared
parameter that receives infusion rates, used in the model as `D(A1) ~ ... + Rin`; leave it
`nothing` for bolus-only routes. `scale` converts dose units to model units:
model amount = `scale` * dose amount.
"""
struct Route <: AbstractDoseTarget
    cmt::Symbol
    state
    rate
    scale::Float64
    function Route(cmt::Symbol, state, rate, scale::Real)
        (isfinite(scale) && scale > 0) || throw(ArgumentError("Route $cmt: scale must be positive and finite, got $scale"))
        return new(cmt, state, rate, Float64(scale))
    end
end
Route(cmt::Symbol, state; rate = nothing, scale = 1.0) = Route(cmt, state, rate, scale)

# ---------------------------------------------------------------------------------------------
# Dose
# ---------------------------------------------------------------------------------------------

"""
    Dose(time, amt; cmt, rate = 0.0, duration = 0.0, ss = false, ii = 0.0)

One dosing event as entered. With `rate = duration = 0` it is a bolus. Give `rate` for a
zero-order infusion of fixed rate (`duration = amt / rate`), or `duration` for an infusion
of fixed duration (`rate = amt / duration`). If both are given they must agree. `amt` must
be non-negative; a zero amount is kept as a zero bolus.

The field `fixed` records which of the two was entered (`:rate`, `:duration`, or `:none` for
a bolus), because bioavailability acts differently on each (see [`DoseControl`](@ref)).

`ss = true` (NONMEM SS=1) says the subject is at steady state on this dose repeated every
`ii` when it is given; see [`Regimen`](@ref). It needs `ii > 0` and `amt > 0`. For a
non-steady-state dose `ii` is ignored and stored as 0.
"""
struct Dose
    time::Float64
    amt::Float64
    cmt::Symbol
    rate::Float64
    duration::Float64
    fixed::Symbol
    ss::Bool
    ii::Float64
    function Dose(time::Real, amt::Real, cmt::Symbol, rate::Real, duration::Real; ss::Bool = false, ii::Real = 0.0)
        isfinite(time) || throw(ArgumentError("dose time must be finite, got $time"))
        (isfinite(amt) && amt >= 0) || throw(ArgumentError("dose amt must be finite and >= 0, got $amt"))
        (isfinite(rate) && rate >= 0) || throw(ArgumentError("dose rate must be finite and >= 0, got $rate"))
        (isfinite(duration) && duration >= 0) || throw(ArgumentError("dose duration must be finite and >= 0, got $duration"))
        if ss
            (isfinite(ii) && ii > 0) || throw(ArgumentError("a steady-state dose (ss = true) needs a finite interval ii > 0, got $ii"))
            amt > 0 || throw(ArgumentError("a steady-state dose needs amt > 0 (constant-infusion steady state, amt = 0 with a rate, is not supported)"))
        end
        iiv = ss ? Float64(ii) : 0.0
        if amt == 0 || (rate == 0 && duration == 0)
            # a bolus (a zero amount never infuses)
            return new(time, amt, cmt, 0.0, 0.0, :none, ss, iiv)
        end
        if rate > 0 && duration > 0
            isapprox(duration, amt / rate; rtol = 1e-8) ||
                throw(ArgumentError("inconsistent infusion: duration $duration != amt / rate = $(amt / rate)"))
            return new(time, amt, cmt, rate, duration, :rate, ss, iiv)
        elseif rate > 0
            return new(time, amt, cmt, rate, amt / rate, :rate, ss, iiv)
        else
            return new(time, amt, cmt, amt / duration, duration, :duration, ss, iiv)
        end
    end
end
Dose(time::Real, amt::Real; cmt::Symbol, rate::Real = 0.0, duration::Real = 0.0, ss::Bool = false, ii::Real = 0.0) =
    Dose(time, amt, cmt, rate, duration; ss, ii)

isinfusion(d::Dose) = d.fixed !== :none

# ---------------------------------------------------------------------------------------------
# Regimen
# ---------------------------------------------------------------------------------------------

"""
    Reset(time)

Reset event (NONMEM EVID=3): at `time` every state goes back to its initial condition, running
infusions are cancelled and doses that are still pending (lagged, not yet in the model) are
cancelled. Combine with doses through [`Regimen`](@ref).
"""
struct Reset
    time::Float64
    function Reset(time::Real)
        isfinite(time) || throw(ArgumentError("reset time must be finite, got $time"))
        return new(time)
    end
end

"""
    Regimen(amt; time = 0.0, cmt, ii = 0.0, addl = 0, rate = 0.0, duration = 0.0, ss = false, reset = false)
    Regimen(r1, r2, ...)       # Regimens and Resets
    vcat(r1, r2, ...)

NONMEM-style dose record: `amt` at `time`, then `addl` more doses every `ii`. `rate` and
`duration` make it an infusion, see [`Dose`](@ref). Several regimens combine into one,
sorted by dose time (doses at equal times keep their order). Get the dose list with
[`doses`](@ref).

  - `ss = true` (NONMEM SS=1, needs `ii > 0`): the FIRST dose of the record is a steady-state
    dose. The subject is taken to have received this dose every `ii` for a very long time (the
    state at `time` is the periodic state of the cycle simulated from the initial conditions),
    including lagged doses and infusions that are still running; the `addl` doses then follow
    as ordinary doses. A regimen can hold several steady-state doses only if they share `time`
    and `ii`. The record resets the state at `time` (see below); everything the regimen did
    before `time` is superseded at `time`.
  - `reset = true` (EVID=4, reset and dose): a [`Reset`](@ref) at `time`, then the doses.
    `Regimen(Reset(t))` is a pure reset (EVID=3), also usable inside `Regimen(r1, Reset(t), r2)`.
  - A reset at the time of a dose applies BEFORE that dose, and the dose (with its lag) is kept.
    A reset cancels doses given before it that have not entered the model yet (lag lands at or
    after the reset time) and stops infusions that are running (the rest of their amount is never
    given). Steady-state doses are not allowed to combine with `reset = true`.
"""
struct Regimen
    doses::Vector{Dose}
    resets::Vector{Float64}
    function Regimen(d::Vector{Dose}, r::Vector{Float64})
        # stable sort by time, so concatenation keeps the order of equal-time doses
        ds = sort(d; by = x -> x.time, alg = Base.Sort.MergeSort)
        ss = [x for x in ds if x.ss]
        if !isempty(ss) && !all(x -> x.time == ss[1].time && x.ii == ss[1].ii, ss)
            throw(ArgumentError("all steady-state doses of a regimen must have the same time and the same ii (one steady-state record per regimen)"))
        end
        return new(ds, sort(unique(r)))
    end
end
Regimen(d::Vector{Dose}) = Regimen(d, Float64[])

function Regimen(amt::Real; time::Real = 0.0, cmt::Symbol, ii::Real = 0.0, addl::Integer = 0, rate::Real = 0.0,
        duration::Real = 0.0, ss::Bool = false, reset::Bool = false)
    addl >= 0 || throw(ArgumentError("addl must be >= 0, got $addl"))
    (isfinite(ii) && ii >= 0) || throw(ArgumentError("ii must be finite and >= 0, got $ii"))
    addl > 0 && ii == 0 && throw(ArgumentError("addl = $addl needs ii > 0"))
    ss && reset && throw(ArgumentError("ss = true already resets the state to the steady state; do not combine it with reset = true"))
    ds = Dose[Dose(time + i * ii, amt, cmt, rate, duration; ss = ss && i == 0, ii) for i in 0:addl]
    return Regimen(ds, reset ? Float64[Reset(time).time] : Float64[])
end
Regimen(r::Reset) = Regimen(Dose[], Float64[r.time])
function Regimen(rs::Union{Regimen, Reset}...)
    ds = Dose[]
    res = Float64[]
    for r in rs
        if r isa Reset
            push!(res, r.time)
        else
            append!(ds, r.doses)
            append!(res, r.resets)
        end
    end
    return Regimen(ds, res)
end
Base.vcat(rs::Regimen...) = Regimen(rs...)

"""
    doses(r::Regimen) -> Vector{Dose}

The expanded dose list of a regimen, sorted by time.
"""
doses(r::Regimen) = r.doses

"""
    reset_times(r::Regimen) -> Vector{Float64}

The reset times of a regimen, sorted.
"""
reset_times(r::Regimen) = r.resets

# ---------------------------------------------------------------------------------------------
# DoseControl
# ---------------------------------------------------------------------------------------------

"""
    DoseControl(; bioav = Dict(), lag = Dict())

Subject-specific dose control, per compartment name (`Dict(:depot => 0.8)`). Missing
entries mean `bioav = 1`, `lag = 0`.

  - `bioav` (F): scales bolus amounts. For an infusion entered with a rate, F scales the
    duration and keeps the rate (NONMEM and Pumas convention). For an infusion entered with
    a duration, F scales the rate and keeps the duration. Either way the infused amount is
    F * amt.
  - `lag`: delays the start of every dose in that compartment (an infusion keeps its length).
"""
struct DoseControl
    bioav::Dict{Symbol, Float64}
    lag::Dict{Symbol, Float64}
    function DoseControl(bioav::Dict{Symbol, Float64}, lag::Dict{Symbol, Float64})
        for (k, v) in bioav
            (isfinite(v) && v >= 0) || throw(ArgumentError("bioav for $k must be finite and >= 0, got $v"))
        end
        for (k, v) in lag
            (isfinite(v) && v >= 0) || throw(ArgumentError("lag for $k must be finite and >= 0, got $v"))
        end
        return new(bioav, lag)
    end
end
DoseControl(; bioav = Dict{Symbol, Float64}(), lag = Dict{Symbol, Float64}()) =
    DoseControl(Dict{Symbol, Float64}(bioav), Dict{Symbol, Float64}(lag))

# ---------------------------------------------------------------------------------------------
# resolve
# ---------------------------------------------------------------------------------------------

"""
    DoseEvent(time, kind, cmt, amount, rate)

One resolved solver event: `kind` is `:bolus`, `:infusion_start` or `:infusion_stop`. Times
are after lag; `amount` (bolus) and `rate` (infusion) are in model units, after F and the
route scale. Infusion events have `amount = 0`, bolus events `rate = 0`.

Two more kinds set the whole state (`cmt` is `:all`, `amount` and `rate` are 0): `:reset`
sets every state back to its value when the solve started (the initial conditions), and
`:steady_state` sets every state to the vector in the sixth field `state` (the steady state;
[`resolve`](@ref) leaves it `nothing` and the simulation fills it in per subject). The state
is set first, then the other events of the same time apply on top of it.
"""
struct DoseEvent
    time::Float64
    kind::Symbol
    cmt::Symbol
    amount::Float64
    rate::Float64
    state::Union{Nothing, Vector{Float64}}
end
DoseEvent(time, kind, cmt, amount, rate) = DoseEvent(time, kind, cmt, amount, rate, nothing)

# Order of coincident events: state sets, then stops, then starts, then boluses.
_rank(kind::Symbol) = (kind === :reset || kind === :steady_state) ? 0 : kind === :infusion_stop ? 1 : kind === :infusion_start ? 2 : 3

function _route_map(routes::Vector{Route})
    m = Dict{Symbol, Route}()
    for r in routes
        haskey(m, r.cmt) && throw(ArgumentError("two routes for compartment $(r.cmt)"))
        m[r.cmt] = r
    end
    return m
end

function _route_of(rmap, d::Dose)
    route = get(rmap, d.cmt, nothing)
    route === nothing && throw(ArgumentError("dose into :$(d.cmt) has no Route (routes: $(collect(keys(rmap))))"))
    return route
end

# A dose given at nominal time `td`, after F and lag: `nothing` (an infusion with F = 0) or
# (kind, start, stop, value), where value is the bolus amount or the infusion rate in model
# units (stop == start for a bolus).
function _landed(d::Dose, td::Real, route::Route, control::DoseControl)
    F = get(control.bioav, d.cmt, 1.0)
    ts = td + get(control.lag, d.cmt, 0.0)
    isinfusion(d) || return (:bolus, ts, ts, route.scale * F * d.amt)
    route.rate === nothing && throw(ArgumentError("infusion into :$(d.cmt) needs a Route with a rate parameter"))
    F == 0 && return nothing
    # fixed-rate: F shortens/lengthens the infusion; fixed-duration: F scales the rate
    rate, dur = d.fixed === :rate ? (d.rate, F * d.duration) : (F * d.rate, d.duration)
    te = ts + dur
    te > ts || throw(ArgumentError(
        "infusion into :$(d.cmt) at time $ts is too short to resolve: duration $dur vanishes in floating point (need t + duration > t)"))
    return (:infusion, ts, te, route.scale * rate)
end

# First barrier time strictly after `t` (Inf if none); `barriers` is sorted.
_cutoff(barriers::Vector{Float64}, t::Real) = (i = searchsortedlast(barriers, t) + 1; i <= length(barriers) ? barriers[i] : Inf)

# Push the events of one landed dose. `cutoff`: a reset (or steady-state jump) at this time
# cancels a dose that lands at or after it and stops an infusion that is still running. `tclip`:
# an infusion that started earlier is switched on at this time instead.
function _emit!(events::Vector{DoseEvent}, cmt::Symbol, l, cutoff::Float64, tclip::Float64)
    l === nothing && return
    kind, ts, te, val = l
    ts >= cutoff && return
    if kind === :bolus
        push!(events, DoseEvent(ts, :bolus, cmt, val, 0.0))
        return
    end
    ts = max(ts, tclip)
    te = min(te, cutoff)
    te > ts || return
    push!(events, DoseEvent(ts, :infusion_start, cmt, 0.0, val))
    push!(events, DoseEvent(te, :infusion_stop, cmt, 0.0, val))
    return
end

# Dose number j (j >= 0) before the steady-state dose `d` (nominal time T - j ii) can still
# matter at/after T until this j; guards against absurd lag / interval ratios.
function _ss_depth(d::Dose, route::Route, control::DoseControl)
    l = _landed(d, d.time, route, control)
    l === nothing && return 0
    span = l[3] - d.time
    n = floor(span / d.ii) + 1
    n <= 1_000_000 || throw(ArgumentError("steady state: lag plus infusion duration ($span) spans more than 10^6 intervals of ii = $(d.ii)"))
    return Int(n)
end

_sortevents!(events) = sort!(events; by = e -> (e.time, _rank(e.kind)), alg = Base.Sort.MergeSort)

"""
    resolve(regimen, routes, control = DoseControl()) -> Vector{DoseEvent}

Turn a regimen into the events the solver needs, for one subject: lag shifts every start,
bioavailability scales amounts (see [`DoseControl`](@ref)), the route scale converts to model
units, and an infusion becomes a start and a stop event. The result is sorted by time. At
the same time infusion stops come first, then infusion starts, then boluses; all boluses at
one time are kept (they add up). An infusion whose stop time equals its start time in floating
point (duration lost next to the dose time) is an error, not silently dropped.

Resets ([`Reset`](@ref), `reset = true`) become `:reset` events, sorted before all other events
of the same time. A dose is cancelled if a reset comes after the time it was GIVEN but not
after it ENTERS the model: a bolus (or infusion start) whose lagged time is at or after the
first reset later than the dose time is dropped; an infusion that is running at that reset is
stopped there. A dose given at the reset time itself is after the reset and kept.

A steady-state dose (`ss = true`) becomes a `:steady_state` event at its time (state still
`nothing`: it depends on the model and the subject, see `simulate`), which acts as a reset at
that time for everything given before. For the doses of the periodic history that are still
pending at that time it adds the events they cause from there on: boluses whose lagged time is
at or after the steady-state time, and infusions that are still running (switched on at the
steady-state time, stopped at their original stop time). Those carried events follow
the same rules for later resets as ordinary doses.
"""
function resolve(regimen::Regimen, routes::Vector{Route}, control::DoseControl = DoseControl())
    rmap = _route_map(routes)
    ssd = [d for d in regimen.doses if d.ss]
    T = isempty(ssd) ? NaN : ssd[1].time
    resets = isempty(ssd) ? regimen.resets : filter(!=(T), regimen.resets)    # a reset at T is superseded by the steady state
    barriers = isempty(ssd) ? resets : sort!(vcat(resets, T))
    events = DoseEvent[]
    for r in resets
        push!(events, DoseEvent(r, :reset, :all, 0.0, 0.0))
    end
    isempty(ssd) || push!(events, DoseEvent(T, :steady_state, :all, 0.0, 0.0))
    for d in regimen.doses
        route = _route_of(rmap, d)
        _emit!(events, d.cmt, _landed(d, d.time, route, control), _cutoff(barriers, d.time), -Inf)
    end
    # doses of the periodic history that are still pending at T (j intervals before)
    cut = isempty(ssd) ? Inf : _cutoff(barriers, T)
    for d in ssd
        route = _route_of(rmap, d)
        for j in 1:_ss_depth(d, route, control)
            l = _landed(d, T - j * d.ii, route, control)
            l === nothing && continue
            (l[1] === :bolus ? l[2] >= T : l[3] > T) || continue
            _emit!(events, d.cmt, l, cut, T)
        end
    end
    return _sortevents!(events)
end

"""
    steady_state_cycle(regimen, routes, control = DoseControl()) -> nothing or (T, ii, events)

The periodic cycle behind a steady-state dose, for the state iteration: `T` and `ii` of the
steady-state record, and the events (boluses, infusion starts and stops, ordered like
[`resolve`](@ref)) of the cycle `[T, T + ii)` of the infinitely repeated record: every dose of
the periodic history whose lagged time falls in the cycle, plus the infusions of earlier doses
that are running at `T` (switched on at `T`). So lag larger than `ii` wraps into the cycle and
an infusion longer than `ii` overlaps itself. Returns `nothing` if the regimen has no
steady-state dose.
"""
function steady_state_cycle(regimen::Regimen, routes::Vector{Route}, control::DoseControl = DoseControl())
    ssd = [d for d in regimen.doses if d.ss]
    isempty(ssd) && return nothing
    rmap = _route_map(routes)
    T, ii = ssd[1].time, ssd[1].ii
    events = DoseEvent[]
    for d in ssd
        route = _route_of(rmap, d)
        for j in 0:_ss_depth(d, route, control)
            l = _landed(d, T - j * ii, route, control)
            l === nothing && continue
            l[2] < T + ii || continue                              # lands in a later cycle
            (l[1] === :bolus ? l[2] >= T : l[3] > T) || continue    # finished before the cycle
            _emit!(events, d.cmt, l, Inf, T)
        end
    end
    return (T, ii, _sortevents!(events))
end

"""
    with_steady_state(events, state) -> Vector{DoseEvent}

The events with `state` (a vector in the order of the system's unknowns) filled into the
`:steady_state` event.
"""
with_steady_state(events::Vector{DoseEvent}, state::Vector{Float64}) =
    [e.kind === :steady_state ? DoseEvent(e.time, e.kind, e.cmt, 0.0, 0.0, state) : e for e in events]

# ---------------------------------------------------------------------------------------------
# dosing_callback
# ---------------------------------------------------------------------------------------------

# One bolus step of a group: add `value` to the state of route number `route`.
struct _Bolus
    route::Int
    value::Float64
end

# Per-integrator data captured when the solve starts: the rate parameters' values and the
# initial state (the target of a reset).
struct _Start
    base::Vector{Float64}
    u0::Vector{Float64}
end

# One infusion: contributes `rate` to its route's rate parameter while start <= t < stop.
struct _Infusion
    start::Float64
    stop::Float64
    rate::Float64
end

"""
    dosing_callback(sys, events, routes) -> (callback, tstops)

Callback for an `ODEProblem` built from the compiled system `sys`, applying the resolved
`events`. A bolus adds its amount to the route's state. For infusions the route's rate
parameter is set at every start or stop to `baseline + sum of the rates of the infusions
running now`, recomputed from scratch, where `baseline` is the parameter value found when the
solve started. So overlapping infusions sum, there is no accumulated rounding, and the
parameter is restored to the baseline exactly when the solve ends (also mid-infusion), which
leaves `sol.prob` as found. The per-solve baseline is kept per integrator, so one callback
object can be used by many solves at once (for example `EnsembleThreads`) and for different
`tspan`s.

`:reset` and `:steady_state` events set all states (in the order of `unknowns(sys)`) before the
other events of their time: a reset to the state the solve started with, a steady state to the
vector stored in the event. They do nothing about infusions or pending doses: those are
cancelled when the events are produced by [`resolve`](@ref).

Events are grouped only when their times are exactly equal; times a rounding error apart stay
separate events. An infusion is on for `start <= t < stop`. Events at the initial time are
applied before the first step; events before it are ignored, except that an infusion already
running at the initial time is on. Events after the end of `tspan` are never reached; an event
exactly at the end fires. Pass `callback` (and `tstops`) to `ODEProblem` or `solve`. `tstops`
lists the event times; the callback schedules them too.
"""
function dosing_callback(sys, events::Vector{DoseEvent}, routes::Vector{Route})
    _route_map(routes)   # duplicate check
    # Rates are computed per route and written to the route's parameter, so two routes
    # writing one parameter would overwrite each other's infusions.
    rates = [r.rate for r in routes if r.rate !== nothing]
    for i in eachindex(rates), j in (i + 1):lastindex(rates)
        isequal(rates[i], rates[j]) && throw(ArgumentError("routes share the rate parameter $(rates[i]); give each route its own"))
    end
    idx = Dict(r.cmt => i for (i, r) in enumerate(routes))
    nr = length(routes)
    getstate = Any[nothing for _ in routes]
    setstate = Any[nothing for _ in routes]
    getrate = Any[nothing for _ in routes]
    setrate = Any[nothing for _ in routes]

    sorted = sort(events; by = e -> (e.time, _rank(e.kind)), alg = Base.Sort.MergeSort)
    # Group by exactly equal time.
    gtimes = Float64[]
    gbolus = Vector{Vector{_Bolus}}()
    grefresh = Vector{Vector{Int}}()                 # routes whose infusion rate changes in the group
    gset = Any[]                                     # nothing, :reset, or the state vector to set
    nstate = length(ModelingToolkit.unknowns(sys))
    needset = any(e -> e.kind === :reset || e.kind === :steady_state, events)
    getall = needset ? getu(sys, ModelingToolkit.unknowns(sys)) : nothing
    setall = needset ? setu(sys, ModelingToolkit.unknowns(sys)) : nothing
    infusions = [_Infusion[] for _ in routes]
    open = [Tuple{Float64, Float64}[] for _ in routes]   # started, not yet stopped: (start, rate)
    for e in sorted
        stateset = e.kind === :reset || e.kind === :steady_state
        stateset || haskey(idx, e.cmt) || throw(ArgumentError("event into :$(e.cmt) has no Route"))
        i = stateset ? 0 : idx[e.cmt]
        if isempty(gtimes) || e.time != gtimes[end]
            push!(gtimes, e.time)
            push!(gbolus, _Bolus[])
            push!(grefresh, Int[])
            push!(gset, nothing)
        end
        if stateset
            gset[end] === nothing || throw(ArgumentError("two state-setting events (reset / steady_state) at time $(e.time)"))
            if e.kind === :reset
                gset[end] = :reset
            else
                e.state === nothing && throw(ArgumentError("steady_state event at time $(e.time) has no state; fill it with with_steady_state (simulate does this per subject)"))
                length(e.state) == nstate || throw(ArgumentError("steady_state event at time $(e.time) has $(length(e.state)) values, the system has $nstate states"))
                gset[end] = e.state
            end
        elseif e.kind === :bolus
            if getstate[i] === nothing
                getstate[i] = getu(sys, routes[i].state)
                setstate[i] = setu(sys, routes[i].state)
            end
            push!(gbolus[end], _Bolus(i, e.amount))
        elseif e.kind === :infusion_start || e.kind === :infusion_stop
            routes[i].rate === nothing && throw(ArgumentError("infusion into :$(routes[i].cmt) needs a Route with a rate parameter"))
            if getrate[i] === nothing
                getrate[i] = getp(sys, routes[i].rate)
                setrate[i] = setp(sys, routes[i].rate)
            end
            i in grefresh[end] || push!(grefresh[end], i)
            if e.kind === :infusion_start
                push!(open[i], (e.time, e.rate))
            else
                # a stop belongs to the earliest open start with the same rate
                k = findfirst(o -> o[2] == e.rate, open[i])
                k === nothing && throw(ArgumentError(
                    "infusion_stop into :$(e.cmt) at time $(e.time) has no earlier infusion_start with rate $(e.rate)"))
                push!(infusions[i], _Infusion(open[i][k][1], e.time, e.rate))
                deleteat!(open[i], k)
            end
        else
            throw(ArgumentError("unknown event kind $(e.kind)"))
        end
    end
    for i in 1:nr, (st, r) in open[i]       # a start without a stop runs to the end
        push!(infusions[i], _Infusion(st, Inf, r))
    end

    # Rate parameter value at time `tau`: baseline plus the infusions running then.
    function rate_at(i, tau, baseline)
        active = 0.0
        on = false
        for inf in infusions[i]
            if inf.start <= tau < inf.stop
                active += inf.rate
                on = true
            end
        end
        return on ? baseline + active : baseline
    end

    # Baselines of the solves in progress, one entry per integrator (so solves do not share).
    baselines = Base.WeakKeyDict{Any, _Start}()

    function affect!(integrator)
        k = searchsortedfirst(gtimes, integrator.t)
        (k <= length(gtimes) && gtimes[k] == integrator.t) || return
        st = gset[k]
        if st !== nothing
            setall(integrator, st === :reset ? baselines[integrator].u0 : st)
        end
        for b in gbolus[k]
            setstate[b.route](integrator, getstate[b.route](integrator) + b.value)
        end
        if !isempty(grefresh[k])
            base = baselines[integrator].base
            for i in grefresh[k]
                setrate[i](integrator, rate_at(i, integrator.t, base[i]))
            end
        end
        derivative_discontinuity!(integrator, true)   # the right-hand side changed: drop the cached derivative
        return
    end

    # Capture the baselines and switch on infusions that are running at the initial time.
    # An integrator that is initialized again without having finished (init followed by reinit!)
    # already holds its baselines; its parameters may be switched on, so do not read them again.
    function init!(c, u, t, integrator)
        old = get(baselines, integrator, nothing)
        base = old === nothing ? zeros(nr) : old.base
        for i in 1:nr
            getrate[i] === nothing && continue
            old === nothing && (base[i] = getrate[i](integrator))
            setrate[i](integrator, rate_at(i, t, base[i]))
        end
        baselines[integrator] = _Start(base, needset ? Float64.(getall(integrator)) : Float64[])
        return
    end

    # Put the rate parameters back exactly as found (the solution's p may be shared).
    function finish!(c, u, t, integrator)
        start = pop!(baselines, integrator, nothing)
        start === nothing && return
        for i in 1:nr
            getrate[i] === nothing || setrate[i](integrator, start.base[i])
        end
        return
    end

    cb = PresetTimeCallback(gtimes, affect!; initialize = init!, finalize = finish!, save_positions = (true, true))
    return cb, copy(gtimes)
end
