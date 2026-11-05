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
be non-negative; a zero amount is kept as a zero bolus (except for a constant infusion, below).

`rate = :model` (NONMEM RATE=-1) takes the infusion rate from the subject's [`DoseControl`](@ref)
(`rate` entry of the compartment), `duration = :model` (RATE=-2) takes the duration from it
(`duration` entry); the dose is then an infusion with that rate or duration. A dose that
needs a modelled value its `DoseControl` does not have is an error when it is resolved.

The field `fixed` records which of the two was entered (`:rate`, `:duration`, `:model_rate`,
`:model_duration`, or `:none` for a bolus), because bioavailability acts differently on each
(see [`DoseControl`](@ref)).

`ss` is the NONMEM SS item: `false`/`0` an ordinary dose, `true`/`1` (SS=1) the subject is at
steady state on this dose repeated every `ii` when it is given (the state is replaced by the
steady state), `2` (SS=2) the steady state of this dose alone is added to the current state.
A steady-state dose needs `ii > 0` and `amt > 0`, except a constant infusion at steady state
(NONMEM AMT=0, RATE>0, SS=1 or 2): `amt = 0` with `rate > 0` (or `rate = :model`) and `ii = 0`, an
infusion that has been running since time minus infinity and ENDS at `time` (NONMEM convention;
stored with `duration = Inf`). For a non-steady-state dose `ii` is ignored and stored as 0. The
field `ss` holds 0, 1 or 2.
"""
struct Dose
    time::Float64
    amt::Float64
    cmt::Symbol
    rate::Float64
    duration::Float64
    fixed::Symbol
    ss::Int
    ii::Float64
    function Dose(time::Real, amt::Real, cmt::Symbol, rate::Union{Real, Symbol}, duration::Union{Real, Symbol};
            ss::Union{Bool, Integer} = false, ii::Real = 0.0)
        ssm = ss_mode(ss)
        isfinite(time) || throw(ArgumentError("dose time must be finite, got $time"))
        (isfinite(amt) && amt >= 0) || throw(ArgumentError("dose amt must be finite and >= 0, got $amt"))
        for (nm, v) in (("rate", rate), ("duration", duration))
            v isa Symbol && v !== :model && throw(ArgumentError("dose $nm must be a number or :model, got $(repr(v))"))
        end
        rate === :model && duration === :model && throw(ArgumentError("rate and duration cannot both be modelled"))
        modelled = rate === :model ? :model_rate : duration === :model ? :model_duration : :no
        if modelled === :model_rate
            duration == 0 || throw(ArgumentError("a modelled rate (rate = :model) cannot be combined with a duration"))
            rate = 0.0
        elseif modelled === :model_duration
            rate == 0 || throw(ArgumentError("a modelled duration (duration = :model) cannot be combined with a rate"))
            duration = 0.0
        end
        (isfinite(rate) && rate >= 0) || throw(ArgumentError("dose rate must be finite and >= 0 (or :model), got $rate"))
        (isfinite(duration) && duration >= 0) || throw(ArgumentError("dose duration must be finite and >= 0 (or :model), got $duration"))
        constant = ssm != 0 && amt == 0 && (modelled === :model_rate || (modelled === :no && rate > 0))
        if constant
            ii == 0 || throw(ArgumentError("a constant infusion at steady state (amt = 0 with a rate, ss = $ssm) needs ii = 0, got $ii"))
            duration == 0 || throw(ArgumentError("a constant infusion at steady state (amt = 0) cannot have a duration"))
            return new(time, 0.0, cmt, modelled === :no ? Float64(rate) : 0.0, Inf, modelled === :no ? :rate : :model_rate, ssm, 0.0)
        end
        if ssm != 0
            (isfinite(ii) && ii > 0) || throw(ArgumentError("a steady-state dose (ss = $ssm) needs a finite interval ii > 0, got $ii"))
            amt > 0 || throw(ArgumentError("a steady-state dose needs amt > 0 (or amt = 0 with a rate and ii = 0 for a constant infusion)"))
        end
        iiv = ssm != 0 ? Float64(ii) : 0.0
        if amt == 0 || (modelled === :no && rate == 0 && duration == 0)
            # a bolus (a zero amount never infuses)
            return new(time, amt, cmt, 0.0, 0.0, :none, ssm, iiv)
        end
        modelled === :no || return new(time, amt, cmt, 0.0, 0.0, modelled, ssm, iiv)
        if rate > 0 && duration > 0
            isapprox(duration, amt / rate; rtol = 1e-8) ||
                throw(ArgumentError("inconsistent infusion: duration $duration != amt / rate = $(amt / rate)"))
            return new(time, amt, cmt, rate, duration, :rate, ssm, iiv)
        elseif rate > 0
            return new(time, amt, cmt, rate, amt / rate, :rate, ssm, iiv)
        else
            return new(time, amt, cmt, amt / duration, duration, :duration, ssm, iiv)
        end
    end
end
Dose(time::Real, amt::Real; cmt::Symbol, rate::Union{Real, Symbol} = 0.0, duration::Union{Real, Symbol} = 0.0,
    ss::Union{Bool, Integer} = false, ii::Real = 0.0) =
    Dose(time, amt, cmt, rate, duration; ss, ii)

# NONMEM SS item: false/true or 0/1/2.
function ss_mode(ss)
    ss isa Bool && return Int(ss)
    ss in (0, 1, 2) || throw(ArgumentError("ss must be false/true or 0, 1, 2 (NONMEM SS); SS=3 and other values are not supported, got $ss"))
    return Int(ss)
end

isinfusion(d::Dose) = d.fixed !== :none
isconstant(d::Dose) = isinf(d.duration)       # constant infusion at steady state (amt = 0)

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
[`doses`](@ref). `rate = :model` / `duration = :model` take the rate / duration from the
subject's [`DoseControl`](@ref) (NONMEM RATE=-1 / -2).

  - `ss = true` or `ss = 1` (NONMEM SS=1, needs `ii > 0`): the FIRST dose of the record is a
    steady-state dose. The subject is taken to have received this dose every `ii` for a very
    long time (the state at `time` is the periodic state of the cycle simulated from the initial
    conditions), including lagged doses and infusions that are still running; the `addl` doses
    then follow as ordinary doses. A regimen can hold several SS=1 doses only if they share
    `time` and `ii` (they form one periodic cycle). The record resets the state at `time` (see
    below); everything the regimen did before `time` is superseded at `time`, and the `addl`
    doses of earlier records that are still to come are cancelled (NONMEM: an SS=1 dose cancels
    pending additional doses), even those whose nominal time is after `time`.
  - `ss = 2` (NONMEM SS=2): the steady state of this dose alone is ADDED to the state at `time`
    instead of replacing it; nothing given before is cancelled. The added amount is the periodic
    state of the dose minus the periodic state without it (both from the initial conditions), so
    a baseline of the model is not counted twice. This is exact for linear models and an
    approximation otherwise (as in NONMEM, which assumes linear kinetics for SS=2). Each SS=2
    record is its own steady state; it may be at any time and have its own `ii`.
  - Constant infusion at steady state (NONMEM AMT=0, RATE>0, SS=1 or 2): `Regimen(0.0; rate = R,
    ss = true)` with `ii = 0` and no `addl`. The state at `time` is the equilibrium under the
    infusion, and the infusion ENDS at `time`, as in NONMEM (the record marks the end of an
    infusion that started at minus infinity). To keep it running, add an ordinary infusion from
    `time` (for example `Regimen(R * 100; time, rate = R, cmt)` for 100 time units).
    Bioavailability does not apply to it (NONMEM) and lag has no effect on it.
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
    record_times::Vector{Float64}     # time of the record each dose comes from (ADDL doses: their record's time)
    function Regimen(d::Vector{Dose}, r::Vector{Float64}, rt::Vector{Float64} = Float64[x.time for x in d])
        length(rt) == length(d) || throw(ArgumentError("one record time per dose needed"))
        # stable sort by time, so concatenation keeps the order of equal-time doses
        perm = sortperm(d; by = x -> x.time, alg = Base.Sort.MergeSort)
        ds = d[perm]
        ss = [x for x in ds if x.ss == 1]
        if !isempty(ss) && !all(x -> x.time == ss[1].time && x.ii == ss[1].ii, ss)
            throw(ArgumentError("all SS=1 steady-state doses of a regimen must have the same time and the same ii (one SS=1 record per regimen; use ss = 2 to add further steady states)"))
        end
        return new(ds, sort(unique(r)), rt[perm])
    end
end
Regimen(d::Vector{Dose}) = Regimen(d, Float64[])

function Regimen(amt::Real; time::Real = 0.0, cmt::Symbol, ii::Real = 0.0, addl::Integer = 0, rate::Union{Real, Symbol} = 0.0,
        duration::Union{Real, Symbol} = 0.0, ss::Union{Bool, Integer} = false, reset::Bool = false)
    ssm = ss_mode(ss)
    addl >= 0 || throw(ArgumentError("addl must be >= 0, got $addl"))
    (isfinite(ii) && ii >= 0) || throw(ArgumentError("ii must be finite and >= 0, got $ii"))
    addl > 0 && ii == 0 && throw(ArgumentError("addl = $addl needs ii > 0"))
    ssm != 0 && reset && throw(ArgumentError("a steady-state dose (ss = $ssm) cannot be combined with reset = true (ss = 1 already resets the state to the steady state)"))
    ds = Dose[Dose(time + i * ii, amt, cmt, rate, duration; ss = i == 0 ? ssm : 0, ii) for i in 0:addl]
    return Regimen(ds, reset ? Float64[Reset(time).time] : Float64[], fill(Float64(time), length(ds)))
end
Regimen(r::Reset) = Regimen(Dose[], Float64[r.time])
function Regimen(rs::Union{Regimen, Reset}...)
    ds = Dose[]
    res = Float64[]
    rt = Float64[]
    for r in rs
        if r isa Reset
            push!(res, r.time)
        else
            append!(ds, r.doses)
            append!(res, r.resets)
            append!(rt, r.record_times)
        end
    end
    return Regimen(ds, res, rt)
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
    DoseControl(; bioav = Dict(), lag = Dict(), rate = Dict(), duration = Dict())

Subject-specific dose control, per compartment name (`Dict(:depot => 0.8)`). Missing
entries mean `bioav = 1`, `lag = 0`. Values that depend on random effects, covariates or the
occasion come from the model instead (`PopulationModel(...; dosecontrol)`, evaluated per dose).

  - `bioav` (F): scales bolus amounts. For an infusion entered with a rate, F scales the
    duration and keeps the rate (NONMEM and Pumas convention). For an infusion entered with
    a duration, F scales the rate and keeps the duration. Either way the infused amount is
    F * amt.
  - `lag`: delays the start of every dose in that compartment (an infusion keeps its length).
  - `rate`, `duration`: the modelled infusion rate (NONMEM R1, R2, ...) and duration (D1, ...)
    used by doses with `rate = :model` and `duration = :model` (NONMEM RATE=-1 and -2). A modelled rate acts like an entered rate
    (F scales the duration), a modelled duration like an entered duration (F scales the rate).
"""
struct DoseControl
    bioav::Dict{Symbol, Float64}
    lag::Dict{Symbol, Float64}
    rate::Dict{Symbol, Float64}
    duration::Dict{Symbol, Float64}
    function DoseControl(bioav::Dict{Symbol, Float64}, lag::Dict{Symbol, Float64},
            rate::Dict{Symbol, Float64} = Dict{Symbol, Float64}(), duration::Dict{Symbol, Float64} = Dict{Symbol, Float64}())
        for (k, v) in bioav
            (isfinite(v) && v >= 0) || throw(ArgumentError("bioav for $k must be finite and >= 0, got $v"))
        end
        for (k, v) in lag
            (isfinite(v) && v >= 0) || throw(ArgumentError("lag for $k must be finite and >= 0, got $v"))
        end
        for (k, v) in rate
            (isfinite(v) && v > 0) || throw(ArgumentError("modelled rate for $k must be finite and > 0, got $v"))
        end
        for (k, v) in duration
            (isfinite(v) && v > 0) || throw(ArgumentError("modelled duration for $k must be finite and > 0, got $v"))
        end
        return new(bioav, lag, rate, duration)
    end
end
DoseControl(; bioav = Dict{Symbol, Float64}(), lag = Dict{Symbol, Float64}(), rate = Dict{Symbol, Float64}(),
    duration = Dict{Symbol, Float64}()) =
    DoseControl(Dict{Symbol, Float64}(bioav), Dict{Symbol, Float64}(lag), Dict{Symbol, Float64}(rate), Dict{Symbol, Float64}(duration))

"""
    merge_control(a::DoseControl, b::DoseControl) -> DoseControl

Union of two dose controls (the subject's own and the one computed by the model's
`dosecontrol`). The same entry (field and compartment) in both is an error.
"""
function merge_control(a::DoseControl, b::DoseControl)
    out = DoseControl()
    for f in (:bioav, :lag, :rate, :duration)
        da, db, d = getfield(a, f), getfield(b, f), getfield(out, f)
        for k in keys(da)
            haskey(db, k) && throw(ArgumentError(
                "$f for :$k is given both by the subject's DoseControl and by the model's dosecontrol; give it in one place"))
        end
        merge!(d, da, db)
    end
    return out
end

# ---------------------------------------------------------------------------------------------
# resolve
# ---------------------------------------------------------------------------------------------

# The DoseControl of dose `d` from the record at time `rt`: `control` itself, or `control(d, rt)`
# for per-dose control (the dose control computed by the model at the record, see PopulationModel).
_ctl(control::DoseControl, d::Dose, rt) = control
function _ctl(control, d::Dose, rt)
    c = control(d, rt)
    c isa DoseControl || throw(ArgumentError("per-dose control must return a DoseControl, got $(typeof(c))"))
    return c
end

"""
    DoseEvent(time, kind, cmt, amount, rate)

One resolved solver event: `kind` is `:bolus`, `:infusion_start` or `:infusion_stop`. Times
are after lag; `amount` (bolus) and `rate` (infusion) are in model units, after F and the
route scale. Infusion events have `amount = 0`, bolus events `rate = 0`.

Two more kinds set the whole state (`cmt` is `:all`, `amount` and `rate` are 0): `:reset`
sets every state back to its value when the solve started (the initial conditions), and
`:steady_state` sets every state to the vector in the sixth field `state` (the steady state;
[`resolve`](@ref) leaves it `nothing` and the simulation fills it in per subject). The state
is set first, then the other events of the same time apply on top of it. A third kind,
`:steady_state_add` (SS=2), adds the vector in `state` to the current state, after any set at
the same time and before the doses.
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
_rank(kind::Symbol) = (kind === :reset || kind === :steady_state || kind === :steady_state_add) ? 0 : kind === :infusion_stop ? 1 : kind === :infusion_start ? 2 : 3

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
    if isconstant(d)
        # running since minus infinity: lag cannot matter, and F does not apply (NONMEM)
        R = d.fixed === :model_rate ? _modelled(control.rate, d, "rate", -1) : d.rate
        return (:infusion, Float64(td), Inf, route.scale * R)
    end
    F == 0 && return nothing
    # fixed-rate: F shortens/lengthens the infusion; fixed-duration: F scales the rate
    rate, dur = if d.fixed === :rate
        (d.rate, F * d.duration)
    elseif d.fixed === :duration
        (F * d.rate, d.duration)
    elseif d.fixed === :model_rate
        R = _modelled(control.rate, d, "rate", -1)
        (R, F * d.amt / R)
    else
        Dm = _modelled(control.duration, d, "duration", -2)
        (F * d.amt / Dm, Dm)
    end
    te = ts + dur
    te > ts || throw(ArgumentError(
        "infusion into :$(d.cmt) at time $ts is too short to resolve: duration $dur vanishes in floating point (need t + duration > t)"))
    return (:infusion, ts, te, route.scale * rate)
end

function _modelled(dict, d::Dose, what, code)
    v = get(dict, d.cmt, nothing)
    v === nothing && throw(ArgumentError("dose into :$(d.cmt) at time $(d.time) has a modelled $what ($what = :model, NONMEM RATE=$code) but the subject's DoseControl gives no $what for :$(d.cmt); pass DoseControl($what = Dict(:$(d.cmt) => ...))"))
    return v
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
    isfinite(te) && push!(events, DoseEvent(te, :infusion_stop, cmt, 0.0, val))    # the cycle of a constant infusion
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
units, and an infusion becomes a start and a stop event. `control` is a `DoseControl` for every
dose, or a function `(dose, record_time) -> DoseControl` giving each dose its own, where
`record_time` is the time of the record the dose comes from (for `addl` doses the time of their
record, so they share its values; a steady-state dose's value applies to its whole periodic
history). The result is sorted by time. At
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
that time for everything given before: doses of earlier records that land at or after it,
including the `addl` doses of an earlier record whose nominal time is at or after it (NONMEM
cancels pending additional doses; the time of the record decides, not the time of each `addl`
dose), are dropped, and running infusions are stopped there. For the doses of the periodic history that are still
pending at that time it adds the events they cause from there on: boluses whose lagged time is
at or after the steady-state time, and infusions that are still running (switched on at the
steady-state time, stopped at their original stop time). Those carried events follow
the same rules for later resets as ordinary doses. An SS=2 dose becomes a `:steady_state_add`
event (one per SS=2 dose, in time order) plus the pending events of its own periodic history; it
cancels nothing. A constant infusion at steady state leaves no event after its time (it ends
there, NONMEM convention); its effect is the steady state alone.
"""
function resolve(regimen::Regimen, routes::Vector{Route}, control = DoseControl())
    rmap = _route_map(routes)
    blocks = ss_blocks(regimen)
    i1 = findfirst(b -> b.mode == 1, blocks)
    T = i1 === nothing ? NaN : blocks[i1].time
    resets = i1 === nothing ? regimen.resets : filter(!=(T), regimen.resets)    # a reset at T is superseded by the steady state
    barriers = i1 === nothing ? resets : sort!(vcat(resets, T))
    events = DoseEvent[]
    for r in resets
        push!(events, DoseEvent(r, :reset, :all, 0.0, 0.0))
    end
    for b in blocks
        push!(events, DoseEvent(b.time, b.mode == 1 ? :steady_state : :steady_state_add, :all, 0.0, 0.0))
    end
    for (d, rt) in zip(regimen.doses, regimen.record_times)
        route = _route_of(rmap, d)
        l = _landed(d, d.time, route, _ctl(control, d, rt))  # also checks a modelled rate is there
        isconstant(d) && continue                      # ends at its time: only its steady state is left
        cut = _cutoff(barriers, d.time)
        # the SS=1 record cancels the pending ADDL doses of every record before it (NONMEM),
        # decided by the time of the record, not by the nominal time of the ADDL dose
        i1 !== nothing && rt < T && (cut = min(cut, T))
        _emit!(events, d.cmt, l, cut, -Inf)
    end
    # doses of each periodic history that are still pending at its time (j intervals before)
    for b in blocks
        cut = _cutoff(barriers, b.time)
        for d in b.doses
            isconstant(d) && continue
            route = _route_of(rmap, d)
            ctl = _ctl(control, d, d.time)         # a steady-state dose is the first dose of its record
            for j in 1:_ss_depth(d, route, ctl)
                l = _landed(d, b.time - j * d.ii, route, ctl)
                l === nothing && continue
                (l[1] === :bolus ? l[2] >= b.time : l[3] > b.time) || continue
                _emit!(events, d.cmt, l, cut, b.time)
            end
        end
    end
    return _sortevents!(events)
end

# One steady state of a regimen: the SS=1 doses (one block, same time and ii) or one SS=2 dose.
# `ii == 0` for a constant infusion.
struct SSBlock
    time::Float64
    mode::Int
    ii::Float64
    doses::Vector{Dose}
end

# The steady-state blocks of a regimen: the SS=1 block first (if any), then the SS=2 doses in
# time order. This is also the order of their events in `resolve` (stable sort).
function ss_blocks(regimen::Regimen)
    blocks = SSBlock[]
    s1 = [d for d in regimen.doses if d.ss == 1]
    isempty(s1) || push!(blocks, SSBlock(s1[1].time, 1, s1[1].ii, s1))
    for d in regimen.doses
        d.ss == 2 && push!(blocks, SSBlock(d.time, 2, d.ii, [d]))
    end
    return blocks
end

# The cycle [T, T + ii) of one block (see steady_state_cycle). For a constant infusion (ii = 0)
# the events are the infusion switched on at T, without a stop.
function block_cycle(b::SSBlock, rmap, control)
    T, ii = b.time, b.ii
    events = DoseEvent[]
    for d in b.doses
        route = _route_of(rmap, d)
        ctl = _ctl(control, d, d.time)
        if isconstant(d)
            _emit!(events, d.cmt, _landed(d, T, route, ctl), Inf, T)
            continue
        end
        for j in 0:_ss_depth(d, route, ctl)
            l = _landed(d, T - j * ii, route, ctl)
            l === nothing && continue
            l[2] < T + ii || continue                              # lands in a later cycle
            (l[1] === :bolus ? l[2] >= T : l[3] > T) || continue    # finished before the cycle
            _emit!(events, d.cmt, l, Inf, T)
        end
    end
    return (T, ii, _sortevents!(events))
end

"""
    steady_state_cycle(regimen, routes, control = DoseControl()) -> nothing or (T, ii, events)

The periodic cycle behind the SS=1 steady-state dose, for the state iteration: `T` and `ii` of
the steady-state record, and the events (boluses, infusion starts and stops, ordered like
[`resolve`](@ref)) of the cycle `[T, T + ii)` of the infinitely repeated record: every dose of
the periodic history whose lagged time falls in the cycle, plus the infusions of earlier doses
that are running at `T` (switched on at `T`). So lag larger than `ii` wraps into the cycle and
an infusion longer than `ii` overlaps itself. For a constant infusion at steady state `ii` is 0
and the events are the infusion switched on at `T`. Returns `nothing` if the regimen has no
SS=1 dose (SS=2 doses are not covered here).
"""
function steady_state_cycle(regimen::Regimen, routes::Vector{Route}, control = DoseControl())
    blocks = ss_blocks(regimen)
    (isempty(blocks) || blocks[1].mode != 1) && return nothing
    return block_cycle(blocks[1], _route_map(routes), control)
end

"""
    with_steady_state(events, state) -> Vector{DoseEvent}

The events with `state` (a vector in the order of the system's unknowns) filled into the
`:steady_state` event.
"""
with_steady_state(events::Vector{DoseEvent}, state::Vector{Float64}) =
    [e.kind === :steady_state ? DoseEvent(e.time, e.kind, e.cmt, 0.0, 0.0, state) : e for e in events]

# Fill the `:steady_state` event with `set` and the k-th `:steady_state_add` event with adds[k].
function fill_steady_states(events::Vector{DoseEvent}, set, adds::Vector)
    k = 0
    return map(events) do e
        e.kind === :steady_state && return DoseEvent(e.time, e.kind, e.cmt, 0.0, 0.0, set)
        e.kind === :steady_state_add || return e
        k += 1
        return DoseEvent(e.time, e.kind, e.cmt, 0.0, 0.0, adds[k])
    end
end

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
vector stored in the event. `:steady_state_add` events (SS=2) then add their vectors. They do nothing about infusions or pending doses: those are
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
    gadd = Any[]                                     # nothing, or the sum of the state vectors to add (SS=2)
    nstate = length(ModelingToolkit.unknowns(sys))
    needset = any(e -> e.kind === :reset || e.kind === :steady_state || e.kind === :steady_state_add, events)
    getall = needset ? getu(sys, ModelingToolkit.unknowns(sys)) : nothing
    setall = needset ? setu(sys, ModelingToolkit.unknowns(sys)) : nothing
    infusions = [_Infusion[] for _ in routes]
    open = [Tuple{Float64, Float64}[] for _ in routes]   # started, not yet stopped: (start, rate)
    for e in sorted
        stateset = e.kind === :reset || e.kind === :steady_state || e.kind === :steady_state_add
        stateset || haskey(idx, e.cmt) || throw(ArgumentError("event into :$(e.cmt) has no Route"))
        i = stateset ? 0 : idx[e.cmt]
        if isempty(gtimes) || e.time != gtimes[end]
            push!(gtimes, e.time)
            push!(gbolus, _Bolus[])
            push!(grefresh, Int[])
            push!(gset, nothing)
            push!(gadd, nothing)
        end
        if e.kind === :steady_state_add
            e.state === nothing && throw(ArgumentError("steady_state_add event at time $(e.time) has no state; simulate fills it per subject"))
            length(e.state) == nstate || throw(ArgumentError("steady_state_add event at time $(e.time) has $(length(e.state)) values, the system has $nstate states"))
            gadd[end] = gadd[end] === nothing ? copy(e.state) : gadd[end] .+ e.state
        elseif stateset
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
        if gadd[k] !== nothing
            setall(integrator, Float64.(getall(integrator)) .+ gadd[k])
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
