# Population model: error models, subjects, population model, eta sampling.
# User guide and semantics: docs/src/simulation.md.
using LinearAlgebra: Symmetric, Diagonal, cholesky, eigen, issuccess, diag
using Random: AbstractRNG, Xoshiro

# ---- Error models -----------------------------------------------------------------------
# σ are standard deviations (not variances). Negative dv is never clipped.

"""
    ErrorModel

Residual error model for one observed endpoint. `sample_dv(err, ipred, rng)` draws one
observation given the individual prediction, always consuming exactly one standard
normal draw from `rng` (also when σ = 0, or `ipred` is NaN), so that the stream stays
aligned. `sd(err, ipred)` is the standard deviation of that observation given `ipred`.

  - `AdditiveError(σ)`:          dv = ipred + σ ε
  - `ProportionalError(σ)`:      dv = ipred (1 + σ ε)
  - `CombinedError(σadd, σprop)`: dv = ipred + sqrt(σadd² + (σprop ipred)²) ε
  - `ExponentialError(σ)`:       dv = ipred exp(σ ε); median-centred (median dv = ipred,
    the mean is ipred exp(σ²/2)), always positive for positive ipred.
"""
abstract type ErrorModel end

function check_sigma(name, x)
    (isfinite(x) && x >= 0) || throw(ArgumentError("$name must be finite and >= 0, got $x"))
    return Float64(x)
end

"""
    AdditiveError(σ)

Additive residual error: `dv = ipred + σ ε`, `ε ~ N(0, 1)`. `σ` is a standard deviation.
"""
struct AdditiveError <: ErrorModel
    σ::Float64
    AdditiveError(σ::Real) = new(check_sigma("σ", σ))
end
"""
    ProportionalError(σ)

Proportional residual error: `dv = ipred (1 + σ ε)`. `σ` is a standard deviation.
"""
struct ProportionalError <: ErrorModel
    σ::Float64
    ProportionalError(σ::Real) = new(check_sigma("σ", σ))
end
"""
    CombinedError(σadd, σprop)

Combined residual error: `dv = ipred + sqrt(σadd^2 + (σprop ipred)^2) ε`.
"""
struct CombinedError <: ErrorModel
    σadd::Float64
    σprop::Float64
    CombinedError(σadd::Real, σprop::Real) = new(check_sigma("σadd", σadd), check_sigma("σprop", σprop))
end
"""
    ExponentialError(σ)

Exponential residual error: `dv = ipred exp(σ ε)`, median-centred.
"""
struct ExponentialError <: ErrorModel
    σ::Float64
    ExponentialError(σ::Real) = new(check_sigma("σ", σ))
end

"""
    sd(err, ipred)

Standard deviation of the observation `dv` given `ipred` under error model `err`.
For `ExponentialError` this is the sd of `dv` itself (not of log dv).
"""
sd(e::AdditiveError, ipred) = e.σ
sd(e::ProportionalError, ipred) = e.σ * abs(ipred)
sd(e::CombinedError, ipred) = sqrt(e.σadd^2 + (e.σprop * ipred)^2)
sd(e::ExponentialError, ipred) = abs(ipred) * sqrt(expm1(e.σ^2) * exp(e.σ^2))

"""
    sample_dv(err, ipred, rng)

One noisy observation for individual prediction `ipred`.
"""
sample_dv(e::AdditiveError, ipred, rng::AbstractRNG) = ipred + e.σ * randn(rng)
sample_dv(e::ProportionalError, ipred, rng::AbstractRNG) = ipred * (1 + e.σ * randn(rng))
sample_dv(e::CombinedError, ipred, rng::AbstractRNG) = ipred + sqrt(e.σadd^2 + (e.σprop * ipred)^2) * randn(rng)
sample_dv(e::ExponentialError, ipred, rng::AbstractRNG) = ipred * exp(e.σ * randn(rng))

# ---- Subject ----------------------------------------------------------------------------

"""
    Subject(id; regimen, covariates = (;), obstimes, control = DoseControl(), timevarying = nothing,
            obsorder = nothing, endpoints = nothing, covariate_rule = :locf, occasions = nothing)

One individual: dosing `regimen` (with subject-specific `control`), baseline `covariates`
(a NamedTuple) and observation times. `timevarying` lists piecewise-constant covariate
changes as `time => NamedTuple` pairs, for example `[10.0 => (WT = 90.0,)]`; covariates not
listed keep their previous value. How a change applies is set by `covariate_rule`:

  - `:locf` (last value carried forward, the default): the values hold from that time onward
    and override the baseline; an observation at that time sees the new values.
  - `:nocb` (next observation carried backward, the NONMEM convention for data records): the
    values hold over the interval ENDING at that time, `(previous change, time]`, so an
    observation at that time still sees them and the next values start just after it. Before
    the first change the first change's values apply; after the last change its values remain.

`occasions` defines occasions for between-occasion variability (IOV, see [`PopulationModel`](@ref)):
a vector of occasion start times (labelled 1, 2, ...) or of `time => label` pairs with integer
labels. The times follow `covariate_rule` like `timevarying`: under `:locf` the occasion starts AT
its time (an observation or a dose at that time is in the new occasion); under `:nocb` (datasets)
the occasion listed at a time holds over the interval ending at that time. Before the first time
the first occasion applies, after the last time the last one. A label that comes back (1, 2, 1)
is the same occasion, with the same random effects. Without `occasions` the whole record is
occasion 1.

`obsorder` tags each observation (same order as `obstimes`) `:pre` or `:post`: a `:pre`
observation sees the state before the dose, reset and steady-state events at its time, a
`:post` observation (the default) the state after them. Covariate values at that time follow
`covariate_rule` either way. `endpoints` restricts each observation to one endpoint (a Symbol,
or `nothing` for every endpoint, the default). Observations are stored sorted by time, `:pre`
before `:post` at equal times, otherwise in the order given.
"""
struct Subject
    id::String
    regimen::Regimen
    covariates::NamedTuple
    covariate_times::Vector{Float64}
    covariate_values::Vector{NamedTuple}
    obstimes::Vector{Float64}
    control::DoseControl
    obsorder::Vector{Symbol}
    obsendpoints::Vector{Union{Nothing, Symbol}}
    covariate_rule::Symbol
    occasion_times::Vector{Float64}
    occasion_ids::Vector{Int}
    function Subject(id::String, regimen::Regimen, covariates::NamedTuple, ct::Vector{Float64}, cv::Vector{NamedTuple},
            obs::Vector{Float64}, control::DoseControl, order::Vector{Symbol}, eps::Vector{Union{Nothing, Symbol}}, rule::Symbol,
            occt::Vector{Float64}, occi::Vector{Int})
        rule in (:locf, :nocb) || throw(ArgumentError("covariate_rule must be :locf or :nocb, got $(repr(rule))"))
        length(order) == length(obs) == length(eps) || throw(ArgumentError("obsorder and endpoints must have one entry per observation time"))
        all(o -> o === :pre || o === :post, order) || throw(ArgumentError("obsorder entries must be :pre or :post"))
        length(occt) == length(occi) || throw(ArgumentError("occasion times and labels must have the same length"))
        all(isfinite, occt) || throw(ArgumentError("occasion times must be finite"))
        all(i -> occt[i] < occt[i + 1], 1:(length(occt) - 1)) || throw(ArgumentError("occasion times must be strictly increasing, got $occt"))
        return new(id, regimen, covariates, ct, cv, obs, control, order, eps, rule, occt, occi)
    end
end
Subject(id::String, regimen::Regimen, covariates::NamedTuple, ct::Vector{Float64}, cv::Vector{NamedTuple},
    obs::Vector{Float64}, control::DoseControl, order::Vector{Symbol}, eps::Vector{Union{Nothing, Symbol}}, rule::Symbol) =
    Subject(id, regimen, covariates, ct, cv, obs, control, order, eps, rule, Float64[], Int[])

# Positional form without the observation tags: every observation `:post`, all endpoints, LOCF.
Subject(id, regimen, covariates, ct, cv, obs, control) =
    Subject(string(id), regimen, covariates, ct, cv, collect(Float64, obs), control, fill(:post, length(obs)),
        Vector{Union{Nothing, Symbol}}(fill(nothing, length(obs))), :locf)

function Subject(id; regimen, covariates = (;), obstimes, control = DoseControl(), timevarying = nothing,
        obsorder = nothing, endpoints = nothing, covariate_rule::Symbol = :locf, occasions = nothing)
    covariates isa NamedTuple || throw(ArgumentError("covariates must be a NamedTuple, got $(typeof(covariates))"))
    t = collect(Float64, obstimes)
    all(isfinite, t) || throw(ArgumentError("obstimes must be finite"))
    n = length(t)
    ord = obsorder === nothing ? fill(:post, n) : collect(Symbol, obsorder)
    eps = endpoints === nothing ? Vector{Union{Nothing, Symbol}}(fill(nothing, n)) : Vector{Union{Nothing, Symbol}}(collect(endpoints))
    length(ord) == n || throw(ArgumentError("obsorder has $(length(ord)) entries for $n observation times"))
    length(eps) == n || throw(ArgumentError("endpoints has $(length(eps)) entries for $n observation times"))
    all(o -> o === :pre || o === :post, ord) || throw(ArgumentError("obsorder entries must be :pre or :post, got $(unique(ord))"))
    perm = sortperm(collect(zip(t, [o === :pre ? 0 : 1 for o in ord])); alg = Base.Sort.MergeSort)
    ct = Float64[]
    cv = NamedTuple[]
    if timevarying !== nothing
        pairs = sort(collect(timevarying); by = first)
        for (time, vals) in pairs
            isfinite(time) || throw(ArgumentError("covariate change time must be finite, got $time"))
            vals isa NamedTuple || throw(ArgumentError("timevarying values must be NamedTuples, got $(typeof(vals))"))
            !isempty(ct) && ct[end] == time && throw(ArgumentError("two covariate changes at time $time"))
            push!(ct, Float64(time))
            push!(cv, vals)
        end
    end
    occt, occi = occasion_list(occasions)
    return Subject(string(id), regimen, covariates, ct, cv, t[perm], control, ord[perm], eps[perm], covariate_rule, occt, occi)
end

# `occasions` keyword -> (times, labels): start times labelled 1, 2, ... or `time => label` pairs.
function occasion_list(occasions)
    occasions === nothing && return Float64[], Int[]
    times = Float64[]
    ids = Int[]
    for (k, o) in enumerate(occasions)
        t, id = o isa Pair ? (first(o), last(o)) : (o, k)
        (t isa Real && isfinite(t)) || throw(ArgumentError("occasion time must be a finite number, got $(repr(t))"))
        (id isa Real && isinteger(id)) || throw(ArgumentError("occasion label must be an integer, got $(repr(id))"))
        push!(times, Float64(t))
        push!(ids, Int(id))
    end
    all(i -> times[i] < times[i + 1], 1:(length(times) - 1)) ||
        throw(ArgumentError("occasion times must be strictly increasing, got $times"))
    return times, ids
end

"""
    with_obstimes(subject, obstimes) -> Subject

The same subject (regimen, covariates and covariate rule, control, occasions) observed at `obstimes`
instead, every observation `:post` and for all endpoints.
"""
with_obstimes(s::Subject, obstimes) =
    Subject(s.id; regimen = s.regimen, covariates = s.covariates, obstimes, control = s.control,
        timevarying = [t => v for (t, v) in zip(s.covariate_times, s.covariate_values)], covariate_rule = s.covariate_rule,
        occasions = [t => i for (t, i) in zip(s.occasion_times, s.occasion_ids)])

"""
    covariates_at(subject, t)

Covariates in force at time `t`. With `covariate_rule = :locf` the baseline values overridden by
every change at or before `t`, in time order. With `:nocb` the baseline overridden by every
change before `t` and then by the first change at or after `t` (the one whose interval contains `t`).
"""
function covariates_at(s::Subject, t::Real)
    c = s.covariates
    for (ti, vi) in zip(s.covariate_times, s.covariate_values)
        if ti <= t
            c = merge(c, vi)
            s.covariate_rule === :nocb && ti == t && break
        else
            s.covariate_rule === :nocb && (c = merge(c, vi))
            break
        end
    end
    return c
end

# Covariates in force just after `t`: under :nocb those of the interval that starts at `t`.
function covariates_after(s::Subject, t::Real)
    s.covariate_rule === :locf && return covariates_at(s, t)
    c = s.covariates
    for (ti, vi) in zip(s.covariate_times, s.covariate_values)
        c = merge(c, vi)
        ti > t && break
    end
    return c
end

"""
    occasion_at(subject, t) -> Int

Label of the occasion in force at time `t`, by the subject's `covariate_rule` (see
[`Subject`](@ref)): under `:locf` the last occasion starting at or before `t`, under `:nocb` the
first occasion listed at or after `t`. Before the first listed time the first occasion, after the
last one the last occasion; 1 if the subject has no occasions.
"""
function occasion_at(s::Subject, t::Real)
    ts, ids = s.occasion_times, s.occasion_ids
    isempty(ts) && return 1
    if s.covariate_rule === :locf
        i = searchsortedlast(ts, t)
        return ids[max(i, 1)]
    end
    i = searchsortedfirst(ts, t)
    return ids[min(i, length(ts))]
end

# Occasion in force just after `t` (under :nocb the one whose interval starts after `t`).
function occasion_after(s::Subject, t::Real)
    s.covariate_rule === :locf && return occasion_at(s, t)
    ts, ids = s.occasion_times, s.occasion_ids
    isempty(ts) && return 1
    i = searchsortedlast(ts, t) + 1
    return ids[min(i, length(ts))]
end

# Occasion labels of a subject, sorted and unique (1 if none).
occasion_labels(s::Subject) = isempty(s.occasion_ids) ? [1] : sort(unique(s.occasion_ids))

# ---- Population model -------------------------------------------------------------------

# Filled in lazily by the first simulation (see simulate.jl): the base ODEProblem,
# built once per model and remade for each subject.
mutable struct SimCache
    lock::ReentrantLock
    prob::Any
    param_keys::Vector{Any}   # symbolic parameters returned by `params`, in order
    init_keys::Vector{Any}    # states returned by `init`, in order
    param_setter::Any         # cached setp(sys, param_keys)
    init_setter::Any          # cached setu(sys, init_keys), or nothing (see subject_remake)
    initial_setter::Any       # cached setp(sys, Initial.(init_keys)), or nothing
    fast_remake::Bool         # per-subject problems by value writes, not symbolic remake
    getters::Vector{Any}      # cached getu(sys, observable), one per endpoint
    time_dependent::Bool      # the compiled system depends explicitly on time
    SimCache() = new(ReentrantLock(), nothing, Any[], Any[], nothing, nothing, nothing, false, Any[], false)
end

"""
    PopulationModel(sys; routes, θ, ηnames, Ω, params, init = (θ, η, c) -> Pair[], observed,
                    error, tspan, solver, solve_kwargs = (;), dosecontrol = nothing,
                    κnames = Symbol[], Ω_IOV = zeros(0, 0))

A structural model `sys` (compiled MTK system) with typical values `θ`, between-subject
random effects `η ~ N(0, Ω)` (names `ηnames`, covariance `Ω` in that order) and residual
error. Individual values come from

  - `params(θ, η, covariates) -> Vector{Pair}` of MTK parameter values (η and covariates
    are NamedTuples),
  - `init(θ, η, covariates) -> Vector{Pair}` of initial conditions.

`observed` and `error` are NamedTuples keyed by endpoint name: symbolic observable (state
or expression) and its [`ErrorModel`](@ref). `routes` map dose compartments to states
(see [`Route`](@ref)); their rate parameters must not be returned by `params`. The final field `cache` holds the lazily built base problem.

Time-varying covariates re-evaluate `params` at each change time and set the parameters
whose value changed. Initial conditions are not re-applied, and only parameters that
`params` returns are updated.

**Dose control from the model.** `dosecontrol(θ, η, covariates)` returns the bioavailability
(F), lag (ALAG), modelled rate (R) and modelled duration (D) per dose compartment, as a
[`DoseControl`](@ref) or a NamedTuple with any of the fields `bioav`, `lag`, `rate`, `duration`
(each a `Dict(compartment => value)`), for example
`(θ, η, c) -> (bioav = Dict(:depot => θ.F * exp(η.F)), lag = Dict(:depot => θ.ALAG * exp(η.ALAG)))`.
It is evaluated once per dose record, at the record's time (before lag): with the covariates
in force at that time (for a dataset dose record, the values on that record) and the random
effects of the trajectory (and of the record's occasion). The values stay with the dose: the
`addl` doses of a record use the record's values, and a lagged dose or a running infusion keeps
them when covariates or the occasion change later (NONMEM: "the value of ALAGn in effect at time
T applies to the dose and to all subsequent additional doses"; mrgsolve does the same for F).
A steady-state dose uses the values of its own record for its whole periodic history. The
subject's own `control` is combined with it; the same entry (field and compartment) in both is
an error.

**Between-occasion variability (IOV).** `κnames` and `Ω_IOV` declare a second random-effect
vector `κ ~ N(0, Ω_IOV)`, drawn independently for every occasion of every subject (occasions:
see [`Subject`](@ref)). Then `params`, `init` and `dosecontrol` are called as
`f(θ, η, covariates, κ)` with the κ of the occasion in force (`init` at the start of the solve,
`dosecontrol` at the time of the dose record); a function with only a three-argument method is called without
κ. Parameters are re-evaluated at every occasion change, like a covariate change.
`Ω_IOV` is checked like `Ω`.
"""
struct PopulationModel{S, P, I, O, C}
    sys::S
    routes::Vector{Route}
    θ::NamedTuple
    ηnames::Vector{Symbol}
    Ω::Matrix{Float64}
    params::P
    init::I
    observed::O
    error::NamedTuple
    tspan_default::Tuple{Float64, Float64}
    solver::Any
    solve_kwargs::NamedTuple
    dosecontrol::Any                # nothing or (θ, η, c[, κ]) -> DoseControl / NamedTuple
    κnames::Vector{Symbol}
    Ω_IOV::Matrix{Float64}
    cache::C
end

function PopulationModel(sys; routes, θ, ηnames, Ω, params, init = (θ, η, c, κ...) -> Pair[], observed, error,
        tspan, solver, solve_kwargs = (;), dosecontrol = nothing, κnames = Symbol[], Ω_IOV = zeros(0, 0))
    θ isa NamedTuple || throw(ArgumentError("θ must be a NamedTuple"))
    observed isa NamedTuple || throw(ArgumentError("observed must be a NamedTuple"))
    error isa NamedTuple || throw(ArgumentError("error must be a NamedTuple"))
    isempty(observed) && throw(ArgumentError("observed is empty"))
    Set(keys(observed)) == Set(keys(error)) ||
        throw(ArgumentError("observed endpoints $(collect(keys(observed))) and error models $(collect(keys(error))) must have the same names"))
    for (k, e) in pairs(error)
        e isa ErrorModel || throw(ArgumentError("error.$k must be an ErrorModel, got $(typeof(e))"))
    end
    names = collect(Symbol, ηnames)
    allunique(names) || throw(ArgumentError("ηnames must be unique: $names"))
    Ωm = validate_omega(Ω, names)
    knames = collect(Symbol, κnames)
    allunique(knames) || throw(ArgumentError("κnames must be unique: $knames"))
    Ωk = validate_omega(Ω_IOV, knames, "Ω_IOV")
    length(tspan) == 2 && tspan[1] < tspan[2] || throw(ArgumentError("tspan must be (t0, tend) with t0 < tend, got $tspan"))
    ids = [r.cmt for r in routes]
    allunique(ids) || throw(ArgumentError("duplicate route compartments: $ids"))
    model = PopulationModel(sys, collect(Route, routes), θ, names, Ωm, params, init, observed, error,
        (Float64(tspan[1]), Float64(tspan[2])), solver, solve_kwargs, dosecontrol, knames, Ωk, SimCache())
    if has_iov(model)
        η0, κ0 = zero_eta(model), zero_kappa(model)
        any(f -> f !== nothing && applicable(f, θ, η0, (;), κ0), (params, dosecontrol)) || throw(ArgumentError(
            "κnames are declared (between-occasion variability) but neither `params` nor `dosecontrol` takes κ: give one of them a method (θ, η, c, κ)"))
    end
    # `params` may need covariates that are not known yet; then the same check runs when the
    # base problem is built (first simulation).
    pp = try
        call_model(params, model, zero_eta(model), (;), zero_kappa(model))
    catch
        nothing
    end
    pp === nothing || check_no_route_rates(model, [first(p) for p in pp])
    return model
end

has_iov(model::PopulationModel) = !isempty(model.κnames)

# `f(θ, η, c)`, or `f(θ, η, c, κ)` when the model has IOV and `f` has a method for it.
function call_model(f, model::PopulationModel, η, c, κ)
    has_iov(model) && applicable(f, model.θ, η, c, κ) && return f(model.θ, η, c, κ)
    return f(model.θ, η, c)
end

# The value returned by a model's `dosecontrol`, as a DoseControl.
function as_dosecontrol(x)
    x isa DoseControl && return x
    x isa NamedTuple || throw(ArgumentError("dosecontrol must return a DoseControl or a NamedTuple with fields bioav, lag, rate, duration; got $(typeof(x))"))
    for k in keys(x)
        k in (:bioav, :lag, :rate, :duration) || throw(ArgumentError("dosecontrol returned the field $k; allowed are bioav, lag, rate, duration"))
    end
    return DoseControl(; x...)
end

# Route rate parameters belong to the dosing callback. A `params` that returned one would
# overwrite the accumulated infusion rate at every covariate change.
function check_no_route_rates(model::PopulationModel, keys)
    for r in model.routes
        r.rate === nothing && continue
        any(k -> isequal(ModelingToolkit.unwrap(k), ModelingToolkit.unwrap(r.rate)), keys) &&
            throw(ArgumentError("`params` returns the infusion rate parameter of route :$(r.cmt) ($(r.rate)); route rates are set by the dosing callback only, remove it from `params`"))
    end
    return nothing
end

"""
    validate_omega(Ω, ηnames, name = "Ω") -> Matrix{Float64}

Check that Ω is square, matches `ηnames`, is symmetric and positive semidefinite
(zero-variance η allowed). Returns the symmetrised matrix. `name` is used in error messages
(`"Ω_IOV"` for the between-occasion covariance).
"""
function validate_omega(Ω, ηnames, name::AbstractString = "Ω")
    n = length(ηnames)
    Ωm = Matrix{Float64}(Ω)
    size(Ωm) == (n, n) ||
        throw(ArgumentError("$name is $(size(Ωm, 1))×$(size(Ωm, 2)) but there are $n $(name == "Ω" ? "ηnames" : "κnames") $ηnames"))
    all(isfinite, Ωm) || throw(ArgumentError("$name contains non-finite entries"))
    n == 0 && return Ωm
    # Tolerances are relative to the largest entry, so a matrix with tiny variances is
    # judged as strictly as one with large variances.
    scale = maximum(abs, Ωm)
    scale == 0 && return Ωm
    maximum(abs, Ωm - Ωm') <= 1e-8 * scale || throw(ArgumentError("$name is not symmetric"))
    Ωm = (Ωm + Ωm') / 2
    all(>=(0), diag(Ωm)) || throw(ArgumentError("$name has a negative variance on its diagonal"))
    λmin = minimum(eigen(Symmetric(Ωm)).values)
    λmin >= -1e-10 * scale ||
        throw(ArgumentError("$name is not positive semidefinite (smallest eigenvalue $λmin, largest entry $scale)"))
    return Ωm
end

# ---- eta sampling -----------------------------------------------------------------------

"""
    eta_factor(Ω) -> Matrix

Matrix `L` with `L L' = Ω`, so that `η = L z` with `z ~ N(0, I)` has covariance Ω.
Zero-variance η get an exactly zero row (their draw is exactly 0). The positive part is
factored by Cholesky; if it is only positive semidefinite (singular, for example two
perfectly correlated η), an eigen-decomposition is used instead.
"""
function eta_factor(Ω::AbstractMatrix)
    n = size(Ω, 1)
    L = zeros(n, n)
    idx = [i for i in 1:n if Ω[i, i] > 0]
    isempty(idx) && return L
    sub = Symmetric(Matrix(Ω[idx, idx]))
    ch = cholesky(sub; check = false)
    if issuccess(ch)
        L[idx, idx] = Matrix(ch.L)
    else
        # exactly singular PSD: eigen-decomposition. Eigenvalues below zero are rounding
        # noise (validate_omega rejects anything material) and are set to 0; the check
        # below makes sure that does not change Ω.
        e = eigen(sub)
        L[idx, idx] = e.vectors * Diagonal(sqrt.(max.(e.values, 0.0)))
    end
    scale = maximum(abs, Ω)
    maximum(abs, L * L' - Ω) <= 1e-8 * scale ||
        throw(ArgumentError("cannot factor Ω: reconstructed covariance differs from Ω"))
    return L
end

"""
    sample_eta(model, rng) -> NamedTuple

Draw η ~ N(0, Ω) as a NamedTuple keyed by `ηnames`. Always consumes `length(ηnames)`
standard normal draws, whatever the structure of Ω.
"""
function sample_eta(model::PopulationModel, rng::AbstractRNG)
    return eta_from_vector(model, eta_factor(model.Ω) * randn(rng, length(model.ηnames)))
end

eta_from_vector(model::PopulationModel, v::AbstractVector) = NamedTuple{Tuple(model.ηnames)}(Tuple(Float64.(v)))

"""
    eta_value(model, η) -> NamedTuple

Normalise a user-given η (NamedTuple with exactly the `ηnames`, or a vector in that order).
"""
function eta_value(model::PopulationModel, η)
    if η isa NamedTuple
        Set(keys(η)) == Set(model.ηnames) ||
            throw(ArgumentError("η names $(collect(keys(η))) do not match model ηnames $(model.ηnames)"))
        return eta_from_vector(model, [η[k] for k in model.ηnames])
    elseif η isa AbstractVector{<:Real}
        length(η) == length(model.ηnames) ||
            throw(ArgumentError("η has length $(length(η)), model has $(length(model.ηnames)) random effects"))
        return eta_from_vector(model, η)
    end
    throw(ArgumentError("η must be :sample, :zero, a NamedTuple or a vector, got $(repr(η))"))
end

zero_eta(model::PopulationModel) = eta_from_vector(model, zeros(length(model.ηnames)))

# ---- kappa (between-occasion random effects) ----------------------------------------------

kappa_from_vector(model::PopulationModel, v::AbstractVector) = NamedTuple{Tuple(model.κnames)}(Tuple(Float64.(v)))
zero_kappa(model::PopulationModel) = kappa_from_vector(model, zeros(length(model.κnames)))

"""
    sample_kappa(model, rng) -> NamedTuple

Draw one occasion's κ ~ N(0, Ω_IOV) as a NamedTuple keyed by `κnames`. Always consumes
`length(κnames)` standard normal draws.
"""
sample_kappa(model::PopulationModel, rng::AbstractRNG) =
    kappa_from_vector(model, eta_factor(model.Ω_IOV) * randn(rng, length(model.κnames)))

"""
    kappa_value(model, κ) -> NamedTuple

Normalise a user-given κ for one occasion (NamedTuple with exactly the `κnames`, or a vector
in that order).
"""
function kappa_value(model::PopulationModel, κ)
    if κ isa NamedTuple
        Set(keys(κ)) == Set(model.κnames) ||
            throw(ArgumentError("κ names $(collect(keys(κ))) do not match model κnames $(model.κnames)"))
        return kappa_from_vector(model, [κ[k] for k in model.κnames])
    elseif κ isa AbstractVector{<:Real}
        length(κ) == length(model.κnames) ||
            throw(ArgumentError("κ has length $(length(κ)), model has $(length(model.κnames)) occasion random effects"))
        return kappa_from_vector(model, κ)
    end
    throw(ArgumentError("κ of an occasion must be a NamedTuple or a vector, got $(repr(κ))"))
end

# The κ of every occasion of a subject, as a Dict label => NamedTuple. `κ` is :sample (each
# occasion from its own stream `kappa_<label>` of `(seed, replicate, id)`), :zero, one value for
# every occasion, or a Dict label => value (every occasion of the subject must be in it).
function subject_kappas(model::PopulationModel, s::Subject, κ, seed, replicate)
    labels = occasion_labels(s)
    out = Dict{Int, NamedTuple}()
    has_iov(model) || return out
    for j in labels
        out[j] = if κ === :sample
            sample_kappa(model, stream_rng(seed, replicate, s.id, Symbol("kappa_", j)))
        elseif κ === :zero
            zero_kappa(model)
        elseif κ isa AbstractDict
            haskey(κ, j) || throw(ArgumentError("subject $(s.id): κ has no value for occasion $j (occasions $labels)"))
            kappa_value(model, κ[j])
        else
            kappa_value(model, κ)
        end
    end
    return out
end

# ---- Reproducible random streams --------------------------------------------------------

# 64-bit FNV-1a over bytes, then a splitmix64 finaliser. Unlike `hash`, this does not
# depend on the Julia version, so seeds stay the same across sessions and releases.
function fnv1a(bytes, h::UInt64 = 0xcbf29ce484222325)
    for b in bytes
        h = (h ⊻ UInt64(b)) * 0x00000100000001b3
    end
    return h
end
function splitmix(z::UInt64)
    z += 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

"""
    stream_seed(seed, replicate, id, stream) -> UInt64

Seed of the random stream `stream` (a Symbol such as `:eta` or `:residual`) of subject
`id` in replicate `replicate`. It depends only on those four values, never on the
subject's position in the cohort, the thread or the order of execution.
"""
function stream_seed(seed::Integer, replicate::Integer, id::AbstractString, stream::Symbol)
    h = splitmix(UInt64(reinterpret(UInt64, Int64(seed))))
    h = splitmix(h ⊻ UInt64(replicate))
    h = splitmix(fnv1a(codeunits(id), h))
    h = splitmix(fnv1a(codeunits(String(stream)), h))
    return h
end
stream_rng(seed, replicate, id, stream) = Xoshiro(stream_seed(seed, replicate, id, stream))

# ---- Virtual populations ----------------------------------------------------------------

"""
    virtual_population(n; regimen, obstimes, covariate_sampler = rng -> (;), seed = 0, control = DoseControl(),
                       occasions = nothing)

`n` subjects with ids `"1"`, `"2"`, ..., the same regimen, observation times and occasions (see
[`Subject`](@ref)), and baseline covariates drawn by `covariate_sampler(rng)` (a NamedTuple).
Subject `i` uses its own stream derived from `seed` and `i`, so adding subjects does not change
earlier ones.
"""
function virtual_population(n::Integer; regimen, obstimes, covariate_sampler = rng -> (;), seed = 0, control = DoseControl(),
        occasions = nothing)
    n >= 0 || throw(ArgumentError("n must be >= 0"))
    return [Subject(string(i); regimen, obstimes, control, occasions,
                covariates = covariate_sampler(stream_rng(seed, 0, string(i), :covariates))) for i in 1:n]
end
