# Percentiles and simulation-based VPC. Semantics: docs/src/simulation.md.
using Statistics: std, median!
#
# Conventions used throughout:
#  - quantiles are Hyndman-Fan type 7 (Statistics.quantile default, linear interpolation of
#    the order statistics, what R's quantile() does);
#  - values that are not finite (failed solves, NaN) are never used, but counted (`n_failed`);
#  - a cell with fewer than `min_n` usable values gives NaN quantiles, with `n` reported.

# ---- Binning ---------------------------------------------------------------------------

# Bins are either the distinct observation times (`exact`, lo == hi == the time) or the
# intervals between edges.
struct Binning
    lo::Vector{Float64}
    hi::Vector{Float64}
    exact::Bool
    index::Dict{Float64, Int}
end

function make_binning(bins, times)
    if bins === :obstimes
        ts = sort(unique(Float64.(times)))
        return Binning(ts, ts, true, Dict(t => i for (i, t) in enumerate(ts)))
    elseif bins isa AbstractVector{<:Real}
        e = collect(Float64, bins)
        (length(e) >= 2 && all(isfinite, e) && all(diff(e) .> 0)) ||
            throw(ArgumentError("bin edges must be finite, strictly increasing and at least two"))
        return Binning(e[1:(end - 1)], e[2:end], false, Dict{Float64, Int}())
    end
    throw(ArgumentError("bins must be :obstimes or a vector of edges, got $(repr(bins))"))
end

# Bin of a time, 0 if it is in no bin. With edges the bins are [e_k, e_{k+1}), except that the
# last one also contains its upper edge, so an observation at the final edge is not lost.
function bin_of(b::Binning, time::Real)
    b.exact && return get(b.index, Float64(time), 0)
    (time < b.lo[1] || time > b.hi[end]) && return 0
    return searchsortedlast(b.lo, time)
end

t_mid(b::Binning, i) = b.exact ? b.lo[i] : (b.lo[i] + b.hi[i]) / 2

# ---- BLQ -------------------------------------------------------------------------------

const BLQ_HANDLING = (:censor_to_lloq, :drop, :keep)

# normalise the blq option to (lloq, handling), or nothing
function check_blq(blq)
    blq === nothing && return nothing
    (blq isa NamedTuple && haskey(blq, :lloq)) || throw(ArgumentError("blq must be nothing or (lloq = x, handling = ...)"))
    lloq = Float64(blq.lloq)
    isfinite(lloq) || throw(ArgumentError("blq.lloq must be finite, got $lloq"))
    handling = get(blq, :handling, :censor_to_lloq)
    handling in BLQ_HANDLING || throw(ArgumentError("blq.handling must be one of $BLQ_HANDLING, got $(repr(handling))"))
    return (; lloq, handling)
end

# ---- Cells: quantiles of one set of values ------------------------------------------------

# Quantiles of `vals` (finite values of one cell). Returns (quantiles, n used, n below lloq).
# BLQ means value < lloq. :censor_to_lloq sets such values to lloq, :drop removes them
# (this biases the quantiles upwards), :keep leaves them; n_blq counts them in every case.
function cell_quantiles(vals::Vector{Float64}, qs, min_n, blq)
    nblq = 0
    if blq !== nothing
        nblq = count(<(blq.lloq), vals)
        if blq.handling === :censor_to_lloq
            vals = max.(vals, blq.lloq)
        elseif blq.handling === :drop
            vals = filter(>=(blq.lloq), vals)
        end
    end
    qv = length(vals) >= min_n ? quantile(vals, collect(Float64, qs)) : fill(NaN, length(qs))
    return qv, length(vals), nblq
end

function check_quantiles(q)
    qs = collect(Float64, q)
    (!isempty(qs) && all(x -> 0 <= x <= 1, qs) && allunique(qs)) ||
        throw(ArgumentError("q must be a non-empty collection of distinct probabilities in [0, 1], got $q"))
    return qs
end

# failure-rate gate: fraction of the values of a cell that were not finite (NaN when the cell has none)
rate_of(nfailed, nfinite) = nfailed + nfinite == 0 ? NaN : nfailed / (nfailed + nfinite)

# max_failure_rate = 0 means report only: nothing is marked
is_unreliable(rate, max_failure_rate) = max_failure_rate > 0 && rate > max_failure_rate     # NaN > x is false

function check_max_failure_rate(x)
    (x isa Real && 0 <= x <= 1) || throw(ArgumentError("max_failure_rate must be in [0, 1] (0 = report only), got $(repr(x))"))
    return Float64(x)
end

"""
    failure_rate(rows) -> Float64
    failure_rate(v::VPCResult) -> Float64

Overall fraction of failed (non-finite) values: the sum of `n_failed` over the bins divided by
the sum of `n_failed + n_finite`, each bin counted once (not once per quantile). For
[`percentiles`](@ref) rows pass the vector it returned; `NaN` if there are no values at all.
The rate per bin is `failure_rate` in each row. A failed solve leaves `NaN` for every
observation it did not reach, so this is the fraction of failed subject trajectories, weighted by
how many observations each one contributes.
"""
function failure_rate(rows::AbstractVector)
    isempty(rows) && return NaN
    first_q = rows[1].q
    once = [r for r in rows if r.q == first_q]
    return rate_of(sum(r.n_failed for r in once), sum(r.n_finite for r in once))
end

# Quantile rows for items = iterable of (replicate, time, value). Replicates are pooled
# (replicate 0) unless by_replicate. Rows are ordered replicate, bin, q; the replicate field is
# present only when by_replicate.
function summarize_items(items, binning::Binning, qs, by_replicate, min_n, blq; max_failure_rate::Real = 0.0)
    vals = Dict{Tuple{Int, Int}, Vector{Float64}}()
    failed = Dict{Tuple{Int, Int}, Int}()
    reps = Set{Int}()
    for (rep, time, v) in items
        by_replicate && push!(reps, rep)
        i = bin_of(binning, time)
        i == 0 && continue
        key = (by_replicate ? rep : 0, i)
        if isfinite(v)
            push!(get!(vals, key, Float64[]), v)
        else
            failed[key] = get(failed, key, 0) + 1
        end
    end
    repids = by_replicate ? sort!(collect(reps)) : [0]
    out = NamedTuple[]
    for rep in repids, i in eachindex(binning.lo)
        v = get(vals, (rep, i), Float64[])
        qv, n, nblq = cell_quantiles(v, qs, min_n, blq)
        nf = get(failed, (rep, i), 0)
        nfin = length(v)
        rate = rate_of(nf, nfin)
        for (qi, qq) in enumerate(qs)
            row = (bin = i, bin_lo = binning.lo[i], bin_hi = binning.hi[i], t_mid = t_mid(binning, i), q = qq, value = qv[qi], n, n_failed = nf, n_blq = nblq,
                n_finite = nfin, failure_rate = rate, unreliable = is_unreliable(rate, max_failure_rate))
            push!(out, by_replicate ? merge((; replicate = rep), row) : row)
        end
    end
    return out
end

function check_value(value)
    value in (:dv, :ipred) || throw(ArgumentError("value must be :dv or :ipred, got $(repr(value))"))
    return value
end

# rows of one endpoint of a SimResult as (replicate, time, value) items
function endpoint_items(r::SimResult, endpoint::Symbol, value::Symbol)
    rows = [o for o in r.obs if o.endpoint === endpoint]
    if isempty(rows)
        have = unique(o.endpoint for o in r.obs)
        throw(ArgumentError("no observations of endpoint $(repr(endpoint)) in the result (endpoints: $(collect(have)))"))
    end
    return [(o.replicate, o.time, Float64(getproperty(o, value))) for o in rows], rows
end

function warn_nonfinite(rows, what)
    nbad = count(r -> r.n_failed > 0, rows)
    nbad == 0 && return nothing
    total = sum(r.n_failed for r in rows)
    @warn "$what: $total non-finite values (failed solves or NaN) were excluded from the quantiles; they are counted in n_failed"
    return nothing
end

"""
    percentiles(r::SimResult, endpoint; q = (0.1, 0.5, 0.9), bins = :obstimes, value = :dv,
                by_replicate = false, min_n = 1, blq = nothing, max_failure_rate = 0.0) -> Vector{NamedTuple}

Quantiles of one endpoint of a simulation, over subjects, per time bin.

Rows (one per bin and q, ordered by bin then q; by replicate first when `by_replicate`):
`(bin, bin_lo, bin_hi, t_mid, q, value, n, n_failed, n_blq, n_finite, failure_rate, unreliable)`,
plus `replicate` first when `by_replicate = true`. By default the replicates are pooled, so `n` counts subjects times
replicates.

  - `value`: `:dv` (with residual error) or `:ipred`.
  - `bins = :obstimes`: one bin per distinct observation time, matched exactly
    (`bin_lo == bin_hi == t_mid` = the time). `bins = edges` (a strictly increasing vector):
    bins `[e_k, e_{k+1})`, half-open, except that the last bin also contains its upper edge;
    observations outside the edges are ignored; `t_mid` is the bin centre.
  - Quantile definition: Hyndman-Fan type 7 (`Statistics.quantile`'s default; linear
    interpolation between order statistics).
  - Values that are not finite (failed solves leave `NaN`) are excluded and counted in
    `n_failed`; one warning reports the total. A bin with fewer than `min_n` usable values
    gives `value = NaN` (its `n` is still reported). The default `min_n = 1` only rejects
    empty bins; quantile estimates from a handful of values are poor, raise it to be strict.
  - `blq = (lloq = x, handling = :censor_to_lloq | :drop | :keep)`: values below `lloq` are
    set to `lloq` (default), removed (this biases quantiles upwards), or kept; `n_blq` counts
    the values below `lloq` in every case, `n` the values that were used.
  - Failure gate: `n_finite` is the number of finite values of the bin (before any BLQ
    handling), `failure_rate = n_failed / (n_failed + n_finite)` (`NaN` for a bin without
    values); [`failure_rate`](@ref) gives the overall rate. `max_failure_rate = 0` (the default)
    only reports. With `max_failure_rate > 0`, bins whose rate is above it get `unreliable = true`
    (and one warning); their `value` is not changed or removed, so a model decides what to do.
"""
function percentiles(r::SimResult, endpoint::Symbol; q = (0.1, 0.5, 0.9), bins = :obstimes, value::Symbol = :dv,
        by_replicate::Bool = false, min_n::Integer = 1, blq = nothing, max_failure_rate::Real = 0.0)
    check_value(value)
    qs = check_quantiles(q)
    mfr = check_max_failure_rate(max_failure_rate)
    items, rows = endpoint_items(r, endpoint, value)
    out = summarize_items(items, make_binning(bins, [o.time for o in rows]), qs, by_replicate, max(min_n, 1), check_blq(blq); max_failure_rate = mfr)
    warn_nonfinite(out, "percentiles")
    warn_unreliable(out, mfr, "percentiles")
    return out
end

function warn_unreliable(rows, max_failure_rate, what)
    max_failure_rate > 0 || return nothing
    firstq = isempty(rows) ? nothing : rows[1].q
    nbad = count(r -> r.unreliable && r.q == firstq, rows)
    nbad == 0 || @warn "$what: $nbad bin(s) have a failure rate above max_failure_rate = $max_failure_rate and are marked unreliable"
    return nothing
end

# ---- VPC -------------------------------------------------------------------------------

"""
    VPCResult

Result of [`vpc`](@ref).

  - `rows`: one per (bin, q), ordered by bin then q:
    `(bin, bin_lo, bin_hi, t_mid, q, median, lo, hi, mc_se, n, n_rep, n_failed, n_blq, failure_rate, unreliable)`. `median`
    is the median over replicates of the replicate's q-quantile in the bin, `lo` and `hi`
    the `pi_level` interval of those quantiles across replicates (the prediction
    interval of the quantile, not a confidence interval of any data), `mc_se` the Monte
    Carlo standard error of `median`. `n` is the median number of usable values per
    replicate, `n_rep` the number of replicates with a defined quantile, `n_failed` the
    total number of excluded non-finite values and `n_blq` the total number of values below
    `lloq` over all replicates. `failure_rate` is the fraction of the simulated values of the bin
    (over all replicates) that were not finite, and `unreliable` is true when it is above
    `max_failure_rate` (never when that is 0, the report-only default).
  - `observed`: `nothing`, or rows `(bin, bin_lo, bin_hi, t_mid, q, value, n, n_failed, n_blq)`
    of the observed quantiles (diagnostic only).
  - `endpoint`, `q`, `pi_level`, `replicates`, `seed`, `value`, `blq`, `n_failed` (all
    excluded simulated values), `n_unbinned` (observed values in no bin), `failure_rate` (overall
    fraction of non-finite simulated values, also `failure_rate(v)`), `max_failure_rate`.
    [`mc_check`](@ref) tests the Monte Carlo precision per row.

Iterating, `length` and indexing act on `rows`.
"""
struct VPCResult
    rows::Vector{NamedTuple}
    observed::Union{Nothing, Vector{NamedTuple}}
    endpoint::Symbol
    q::Vector{Float64}
    pi_level::Float64
    replicates::Int
    seed::Union{Nothing, Int}
    value::Symbol
    blq::Union{Nothing, NamedTuple}
    n_failed::Int
    n_unbinned::Int
    failure_rate::Float64
    max_failure_rate::Float64
end

VPCResult(rows, observed, endpoint, q, pi_level, replicates, seed, value, blq, n_failed, n_unbinned) =
    VPCResult(rows, observed, endpoint, q, pi_level, replicates, seed, value, blq, n_failed, n_unbinned, NaN, 0.0)

failure_rate(v::VPCResult) = v.failure_rate

Base.length(v::VPCResult) = length(v.rows)
Base.iterate(v::VPCResult, st...) = iterate(v.rows, st...)
Base.getindex(v::VPCResult, i) = v.rows[i]
function Base.show(io::IO, v::VPCResult)
    print(io, "VPCResult(", v.endpoint, ", ", length(unique(r.bin for r in v.rows)), " bins, q = ", Tuple(v.q), ", ",
        v.replicates, " replicates", v.observed === nothing ? "" : ", with observed", v.n_failed == 0 ? "" : ", $(v.n_failed) values excluded", ")")
end

"""
    mc_se_median(x, rng; B = 500, warn = true) -> Float64

Monte Carlo standard error of the median of `x`, the per-replicate quantiles of one bin
(treated as independent draws), by the bootstrap over replicates: `B` resamples of `x` with
replacement, the standard deviation of their medians. It needs no density estimate and, unlike
rank-interval formulas, it stays honest for tied or discrete values (BLQ censoring puts many
replicate quantiles at exactly `lloq`; a rank interval then has width zero and says "no error").
`rng` is a random number generator, seeded by `vpc` from its `seed`, so the result is
reproducible. `NaN` for fewer than 10 values. A bootstrap SE of a median is itself noisy
(about 10-20 % at K = 100-400) and slightly rough for discrete values; treat it as a noise
level, not an exact number. It is exactly 0 only when all values are equal (no spread at all).

Degenerate case: with heavily tied values (say 98 % of the replicate quantiles at `lloq`) every
resample can have the same median although `x` is not constant, and the standard deviation of the
medians would be 0, "no error", which is false. When `x` is not constant, fewer than 3 distinct
resampled medians occur and fewer than 1 % of the resamples (at least one) differ from the commonest
median, the result is `NaN` and, with `warn = true`, a warning says so. Use more replicates (or fewer
ties) to get a number. Two distinct medians that both occur often (a balanced 0/1 sample of odd
size) are a usable estimate and are not flagged. [`vpc`](@ref) calls this with `warn = false` and reports
the number of such cells in one warning.
"""
function mc_se_median(x::AbstractVector{<:Real}, rng::AbstractRNG; B::Integer = 500, warn::Bool = true)
    se, degenerate = mc_se_median_checked(x, rng; B)
    degenerate && warn && @warn "mc_se_median: the bootstrap medians of $(length(x)) values that are not all equal (heavily tied) hardly vary; the standard error is not estimable and is reported as NaN, not 0"
    return se
end

# (se, degenerate): degenerate is true when the bootstrap cannot say anything (see mc_se_median)
function mc_se_median_checked(x::AbstractVector{<:Real}, rng::AbstractRNG; B::Integer = 500)
    K = length(x)
    K < 10 && return (NaN, false)
    buf = Vector{Float64}(undef, K)
    meds = Vector{Float64}(undef, B)
    for b in 1:B
        for i in 1:K
            buf[i] = x[rand(rng, 1:K)]
        end
        meds[b] = median!(buf)
    end
    # degenerate: the resampled medians (almost) never move although x is not constant
    if length(unique(x)) > 1
        vals = unique(meds)
        if length(vals) < 3
            offmode = B - maximum(count(==(v), meds) for v in vals)
            offmode < max(1, B ÷ 100) && return (NaN, true)
        end
    end
    return (std(meds), false)
end

"""
    mc_check(v::VPCResult, tolerance; fraction = 0.2) -> Vector{NamedTuple}

Is the Monte Carlo precision of a [`vpc`](@ref) good enough for a comparison against `tolerance`?
One row per (bin, q) of `v`: `(bin, t_mid, q, mc_se, limit, ok)` with `limit = fraction * tolerance`
and `ok = mc_se < limit` (false when `mc_se` is `NaN`, i.e. fewer than 10 replicates or too
many ties). `tolerance` is in the units of the value (for example the allowed deviation of the
median in concentration units); pass a function `row -> tolerance` to give each bin or quantile its
own, for instance when the tolerance is a fraction of the value. A model asserts precision with
`all(c.ok for c in mc_check(v, tol))` and raises `replicates` if it fails. Nothing in the atlas
uses this for grading yet.
"""
function mc_check(v::VPCResult, tolerance; fraction::Real = 0.2)
    (0 < fraction <= 1) || throw(ArgumentError("fraction must be in (0, 1], got $fraction"))
    tolf = tolerance isa Function ? tolerance : (_ -> tolerance)
    return map(v.rows) do r
        tol = tolf(r)
        (tol isa Real && tol > 0) || throw(ArgumentError("tolerance must be a positive number, got $(repr(tol))"))
        limit = fraction * tol
        (bin = r.bin, t_mid = r.t_mid, q = r.q, mc_se = r.mc_se, limit = limit, ok = r.mc_se < limit)
    end
end

"""
    vpc(model, subjects; endpoint, q = (0.1, 0.5, 0.9), bins = :obstimes, replicates = 200,
        pi_level = 0.95, seed = 0, value = :dv, blq = nothing, observed = nothing, min_n = 1,
        max_failure_rate = 0.0, simulate_kwargs...) -> VPCResult
    vpc(r::SimResult; endpoint, ...)

Simulation-based visual predictive check. The cohort `subjects` is simulated `replicates`
times (`simulate(model, subjects; replicates, seed, simulate_kwargs...)`, so it is
reproducible from `seed` and independent of threads). In every replicate the `q`-quantile of
`endpoint` is taken per bin over the subjects ([`percentiles`](@ref) with `by_replicate`:
same `bins`, `value`, quantile definition, failed value handling, `min_n` and `blq`). Per
bin and q the rows then give the median of those replicate quantiles, the `pi_level`
interval across replicates (Hyndman-Fan type 7 quantiles of the K replicate values, at
`(1 - pi_level) / 2` and `(1 + pi_level) / 2`) and `mc_se`, the Monte Carlo standard error
of that median (`QSPAtlas.mc_se_median`: bootstrap over the replicates, 500 resamples, generator
seeded from `seed` and the bin, so it is reproducible and handles ties from BLQ censoring). `mc_se` shrinks like
`1 / sqrt(replicates)`; the interval `lo`..`hi` does not shrink, it estimates the spread of
the quantile for a cohort of this size. Only replicates with a defined quantile count (`n_rep`).

`value = :dv` includes the residual error, `:ipred` does not. For a population target
freeze `subjects`, `bins`, `value`, `replicates` and `seed` before comparing with a
reference, and keep `mc_se` small next to the tolerance.

`observed` is optional observed data, a vector of `(id, time, dv)` named tuples. It gives the
observed quantiles per bin (pooled over subjects, same `bins`, `min_n`, `blq`), stored in
`VPCResult.observed`. This is a diagnostic: nothing in the atlas is graded on it. With
`bins = :obstimes`, observed times must equal simulated times exactly to be binned; the rest
is counted in `n_unbinned`.

Failed or NaN simulated values are excluded and counted (`n_failed`, one warning). Each row gives
the `failure_rate` of its bin; `max_failure_rate > 0` marks bins above it `unreliable` (values
unchanged; the default 0 only reports). `mc_se` is `NaN` (with one warning) when the bootstrap
cannot estimate it from heavily tied replicate quantiles, see [`mc_se_median`](@ref QSPAtlas.mc_se_median);
[`mc_check`](@ref) tests `mc_se` against a tolerance.
"""
function vpc(model::PopulationModel, subjects::AbstractVector{<:Subject}; endpoint::Symbol, q = (0.1, 0.5, 0.9), bins = :obstimes,
        replicates::Integer = 200, pi_level::Real = 0.95, seed::Integer = 0, value::Symbol = :dv, blq = nothing,
        observed = nothing, min_n::Integer = 1, max_failure_rate::Real = 0.0, simulate_kwargs...)
    check_max_failure_rate(max_failure_rate)
    replicates >= 1 || throw(ArgumentError("replicates must be >= 1"))
    check_vpc_options(; q, bins, pi_level, value, blq)          # fail before the expensive simulation
    res = simulate(model, subjects; replicates, seed, simulate_kwargs...)
    return vpc(res; endpoint, q, bins, pi_level, seed, value, blq, observed, min_n, max_failure_rate)
end

function check_vpc_options(; q, bins, pi_level, value, blq)
    check_value(value)
    check_quantiles(q)
    check_blq(blq)
    (0 < pi_level < 1) || throw(ArgumentError("pi_level must be in (0, 1), got $pi_level"))
    bins === :obstimes || make_binning(bins, Float64[])
    return nothing
end

function vpc(r::SimResult; endpoint::Symbol, q = (0.1, 0.5, 0.9), bins = :obstimes, pi_level::Real = 0.95,
        seed = nothing, value::Symbol = :dv, blq = nothing, observed = nothing, min_n::Integer = 1, max_failure_rate::Real = 0.0)
    check_vpc_options(; q, bins, pi_level, value, blq)
    mfr = check_max_failure_rate(max_failure_rate)
    qs = check_quantiles(q)
    blq = check_blq(blq)
    min_n = max(min_n, 1)
    items, rows = endpoint_items(r, endpoint, value)
    binning = make_binning(bins, [o.time for o in rows])
    reps = summarize_items(items, binning, qs, true, min_n, blq)
    # gather the replicate quantiles per (bin, q)
    cells = Dict{Tuple{Int, Float64}, Vector{NamedTuple}}()
    for row in reps
        push!(get!(cells, (row.bin, row.q), NamedTuple[]), row)
    end
    plo, phi = (1 - pi_level) / 2, (1 + pi_level) / 2
    out = NamedTuple[]
    ndegenerate = 0
    for i in eachindex(binning.lo), qq in qs
        cell = cells[(i, qq)]
        vals = Float64[c.value for c in cell if isfinite(c.value)]
        md, lo, hi = isempty(vals) ? (NaN, NaN, NaN) : (median(vals), quantile(vals, plo), quantile(vals, phi))
        se, degenerate = mc_se_median_checked(vals, stream_rng(seed === nothing ? 0 : seed, i, string(qq), :mcse))
        ndegenerate += degenerate
        nf, nfin = sum(c.n_failed for c in cell), sum(c.n_finite for c in cell)
        rate = rate_of(nf, nfin)
        push!(out, (bin = i, bin_lo = binning.lo[i], bin_hi = binning.hi[i], t_mid = t_mid(binning, i), q = qq,
            median = md, lo, hi, mc_se = se, n = round(Int, median([c.n for c in cell])),
            n_rep = length(vals), n_failed = nf, n_blq = sum(c.n_blq for c in cell),
            failure_rate = rate, unreliable = is_unreliable(rate, mfr)))
    end
    ndegenerate > 0 && @warn "vpc: mc_se is NaN in $ndegenerate (bin, q) cell(s) where the replicate quantiles are heavily tied (the bootstrap medians hardly vary); it is not estimable there, not 0"
    nfail = sum(row.n_failed for row in out if row.q == qs[1]; init = 0)
    overall = failure_rate([(; q = row.q, n_failed = row.n_failed, n_finite = sum(c.n_finite for c in cells[(row.bin, row.q)])) for row in out])
    nfail > 0 && @warn "vpc: $nfail non-finite simulated values (failed solves or NaN) were excluded; they are counted in n_failed"
    warn_unreliable(out, mfr, "vpc")
    obs_rows, unbinned = nothing, 0
    if observed !== nothing
        oitems = [(0, Float64(o.time), Float64(o.dv)) for o in observed]
        unbinned = count(o -> bin_of(binning, o[2]) == 0, oitems)
        obs_rows = summarize_items(oitems, binning, qs, false, min_n, blq)
    end
    nrep = length(unique(s.replicate for s in r.subjects))
    return VPCResult(out, obs_rows, endpoint, qs, Float64(pi_level), nrep, seed === nothing ? nothing : Int(seed), value, blq, nfail, unbinned, overall, mfr)
end
