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
    Subject(id; regimen, covariates = (;), obstimes, control = DoseControl(), timevarying = nothing)

One individual: dosing `regimen` (with subject-specific `control`), baseline `covariates`
(a NamedTuple) and observation times. `timevarying` lists piecewise-constant covariate
changes as `time => NamedTuple` pairs, for example `[10.0 => (WT = 90.0,)]`: the values
hold from that time onward and override the baseline (covariates not listed keep their
previous value). `obstimes` are stored sorted.
"""
struct Subject
    id::String
    regimen::Regimen
    covariates::NamedTuple
    covariate_times::Vector{Float64}
    covariate_values::Vector{NamedTuple}
    obstimes::Vector{Float64}
    control::DoseControl
end

function Subject(id; regimen, covariates = (;), obstimes, control = DoseControl(), timevarying = nothing)
    covariates isa NamedTuple || throw(ArgumentError("covariates must be a NamedTuple, got $(typeof(covariates))"))
    obs = sort(collect(Float64, obstimes))
    all(isfinite, obs) || throw(ArgumentError("obstimes must be finite"))
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
    return Subject(string(id), regimen, covariates, ct, cv, obs, control)
end

"""
    covariates_at(subject, t)

Covariates in force at time `t`: the baseline values overridden by every change at or
before `t`, in time order.
"""
function covariates_at(s::Subject, t::Real)
    c = s.covariates
    for (ti, vi) in zip(s.covariate_times, s.covariate_values)
        ti <= t || break
        c = merge(c, vi)
    end
    return c
end

# ---- Population model -------------------------------------------------------------------

# Filled in lazily by the first simulation (see simulate.jl): the base ODEProblem,
# built once per model and remade for each subject.
mutable struct SimCache
    lock::ReentrantLock
    prob::Any
    param_keys::Vector{Any}   # symbolic parameters returned by `params`, in order
    init_keys::Vector{Any}    # states returned by `init`, in order
    param_setter::Any         # cached setp(sys, param_keys)
    getters::Vector{Any}      # cached getu(sys, observable), one per endpoint
    time_dependent::Bool      # the compiled system depends explicitly on time
    SimCache() = new(ReentrantLock(), nothing, Any[], Any[], nothing, Any[], false)
end

"""
    PopulationModel(sys; routes, θ, ηnames, Ω, params, init = (θ, η, c) -> Pair[], observed,
                    error, tspan, solver, solve_kwargs = (;))

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
    cache::C
end

function PopulationModel(sys; routes, θ, ηnames, Ω, params, init = (θ, η, c) -> Pair[], observed, error,
        tspan, solver, solve_kwargs = (;))
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
    length(tspan) == 2 && tspan[1] < tspan[2] || throw(ArgumentError("tspan must be (t0, tend) with t0 < tend, got $tspan"))
    ids = [r.cmt for r in routes]
    allunique(ids) || throw(ArgumentError("duplicate route compartments: $ids"))
    model = PopulationModel(sys, collect(Route, routes), θ, names, Ωm, params, init, observed, error,
        (Float64(tspan[1]), Float64(tspan[2])), solver, solve_kwargs, SimCache())
    # `params` may need covariates that are not known yet; then the same check runs when the
    # base problem is built (first simulation).
    pp = try
        params(θ, zero_eta(model), (;))
    catch
        nothing
    end
    pp === nothing || check_no_route_rates(model, [first(p) for p in pp])
    return model
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
    validate_omega(Ω, ηnames) -> Matrix{Float64}

Check that Ω is square, matches `ηnames`, is symmetric and positive semidefinite
(zero-variance η allowed). Returns the symmetrised matrix.
"""
function validate_omega(Ω, ηnames)
    n = length(ηnames)
    Ωm = Matrix{Float64}(Ω)
    size(Ωm) == (n, n) ||
        throw(ArgumentError("Ω is $(size(Ωm, 1))×$(size(Ωm, 2)) but there are $n ηnames $ηnames"))
    all(isfinite, Ωm) || throw(ArgumentError("Ω contains non-finite entries"))
    n == 0 && return Ωm
    # Tolerances are relative to the largest entry, so a matrix with tiny variances is
    # judged as strictly as one with large variances.
    scale = maximum(abs, Ωm)
    scale == 0 && return Ωm
    maximum(abs, Ωm - Ωm') <= 1e-8 * scale || throw(ArgumentError("Ω is not symmetric"))
    Ωm = (Ωm + Ωm') / 2
    all(>=(0), diag(Ωm)) || throw(ArgumentError("Ω has a negative variance on its diagonal"))
    λmin = minimum(eigen(Symmetric(Ωm)).values)
    λmin >= -1e-10 * scale ||
        throw(ArgumentError("Ω is not positive semidefinite (smallest eigenvalue $λmin, largest entry $scale)"))
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
    virtual_population(n; regimen, obstimes, covariate_sampler = rng -> (;), seed = 0, control = DoseControl())

`n` subjects with ids `"1"`, `"2"`, ..., the same regimen and observation times, and
baseline covariates drawn by `covariate_sampler(rng)` (a NamedTuple). Subject `i` uses its
own stream derived from `seed` and `i`, so adding subjects does not change earlier ones.
"""
function virtual_population(n::Integer; regimen, obstimes, covariate_sampler = rng -> (;), seed = 0, control = DoseControl())
    n >= 0 || throw(ArgumentError("n must be >= 0"))
    return [Subject(string(i); regimen, obstimes, control,
                covariates = covariate_sampler(stream_rng(seed, 0, string(i), :covariates))) for i in 1:n]
end
