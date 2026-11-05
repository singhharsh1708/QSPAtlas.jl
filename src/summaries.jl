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

const BLQ_HANDLING = (:censor_to_lloq, :drop, :keep, :censored)
const WEIGHTINGS = (:observation, :subject)
const SUBJECT_VALUES = (:mean, :first, :last)

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

function check_weighting(weighting, subject_value)
    weighting in WEIGHTINGS || throw(ArgumentError("weighting must be one of $WEIGHTINGS, got $(repr(weighting))"))
    subject_value in SUBJECT_VALUES || throw(ArgumentError("subject_value must be one of $SUBJECT_VALUES, got $(repr(subject_value))"))
    return nothing
end

# ---- Cells: quantiles of one set of values ------------------------------------------------

# Items are (replicate, time, value, censored, id): `censored` marks a value below the limit
# of quantification. `prepare_items` derives the flag (value < lloq, or an explicit flag for
# observed data) and, for :censor_to_lloq, sets the censored values to lloq. With :censored the
# value of a censored item is unknown, so its slot holds its censoring limit instead (lloq here;
# a prediction correction then transforms it like a value, giving that record's own limit on
# the corrected scale). It is done before any prediction correction, so the flag always refers
# to the raw scale.
function prepare_items(raw, ids, blq; flags = nothing)
    out = Vector{Tuple{Int, Float64, Float64, Bool, String}}(undef, length(raw))
    for (k, (rep, time, v)) in enumerate(raw)
        c = blq !== nothing && isfinite(v) && (flags === nothing ? v < blq.lloq : flags[k])
        w = (c && (blq.handling === :censor_to_lloq || blq.handling === :censored)) ? blq.lloq : v
        out[k] = (rep, Float64(time), w, c, ids[k])
    end
    return out
end

# Hyndman-Fan type 7 quantile of n values whose lowest `ncens` are censored (below the limit):
# `sorted` are the uncensored values, ascending. The quantile interpolates the order statistics
# floor(h) and ceil(h), h = (n - 1) p + 1; when floor(h) is censored the quantile cannot be
# computed (it is below the limit, or depends on a value that is only known to be below it):
# returns (NaN, true). Returns (NaN, false) for no values. Without censoring it is `quantile`.
function censored_quantile(sorted::AbstractVector{Float64}, ncens::Integer, p::Real)
    n = ncens + length(sorted)
    n == 0 && return (NaN, false)
    h = (n - 1) * p + 1
    r = round(h)
    abs(h - r) < 1e-9 && (h = r)                       # guard the floor against rounding noise
    lo, hi = floor(Int, h), ceil(Int, h)
    lo <= ncens && return (NaN, true)
    x1, x2 = sorted[lo - ncens], sorted[hi - ncens]
    return (x1 + (h - lo) * (x2 - x1), false)
end

# The same with a censoring limit per censored value: `sorted` are the quantified values,
# ascending, `limits` the limits of the censored ones (each is only known to be below its own
# limit; after a prediction correction or a subject mean the limits differ, and a censored value
# need not lie below every quantified one). Order statistic k is known when it is the same for
# every value the censored ones can take: k > ncens, x = sorted[k - ncens], and at most as many
# limits lie above x as there are further quantified values tied with x. The type 7 quantile is
# a number when both order statistics it interpolates are known; otherwise it is (NaN, true,
# bound) with `bound` the largest value the lower order statistic can take (limits taken as
# values). With every limit equal to lloq and quantified values >= lloq this is
# `censored_quantile`, with bound = lloq. Returns (NaN, false, NaN) for no values.
function censored_quantile_limits(sorted::AbstractVector{Float64}, limits::AbstractVector{Float64}, p::Real)
    ncens = length(limits)
    n = ncens + length(sorted)
    n == 0 && return (NaN, false, NaN)
    h = (n - 1) * p + 1
    r = round(h)
    abs(h - r) < 1e-9 && (h = r)
    lo, hi = floor(Int, h), ceil(Int, h)
    function known(k)
        k <= ncens && return false
        x = sorted[k - ncens]
        above = count(>(x), limits)
        above == 0 && return true
        ties = count(==(x), @view sorted[(k - ncens + 1):end])
        return above <= ties
    end
    if known(lo) && known(hi)
        x1, x2 = sorted[lo - ncens], sorted[hi - ncens]
        return (x1 + (h - lo) * (x2 - x1), false, NaN)
    end
    return (NaN, true, sort!(vcat(sorted, limits))[lo])
end

# Quantiles of the values of one cell (finite, after weighting) with their censoring flags.
# Returns (quantiles, n used, n censored, below flags, limits). BLQ means value < lloq.
# :censor_to_lloq was applied upstream, :drop removes the censored values (this biases the
# quantiles upwards), :keep leaves them, :censored keeps them as censored, each below its own
# limit (the value slot of a censored item holds its limit, see prepare_items), so a quantile
# that depends on one is NaN with its flag set and its bound in `limits` (NaN elsewhere).
# n_blq counts the censored values in every case.
function cell_quantiles(vals::Vector{Float64}, cens::Vector{Bool}, qs, min_n, blq)
    nblq = count(cens)
    handling = blq === nothing ? :keep : blq.handling
    nq = length(qs)
    if handling === :drop
        vals = vals[.!cens]
    end
    n = length(vals)
    if handling === :censored
        if n >= min_n
            unc = sort(vals[.!cens])
            lims = vals[cens]
            res = [censored_quantile_limits(unc, lims, q) for q in qs]
            return Float64[r[1] for r in res], n, nblq, Bool[r[2] for r in res], Float64[r[3] for r in res]
        end
        return fill(NaN, nq), n, nblq, fill(false, nq), fill(NaN, nq)
    end
    qv = n >= min_n ? quantile(vals, collect(Float64, qs)) : fill(NaN, nq)
    return qv, n, nblq, fill(false, nq), fill(NaN, nq)
end

# One value per subject: its mean, first or last value in the cell (weighting = :subject).
# A mean that includes a censored value is itself censored: it is unknown, whatever placeholder
# the censored value has. Its slot is the mean of the slots, which for :censored (slots of
# censored values hold their limits) is the limit of the mean.
function collapse_cell(cell, weighting::Symbol, subject_value::Symbol)
    if weighting === :observation
        return Float64[c[1] for c in cell], Bool[c[3] for c in cell]
    end
    # one value per subject and replicate (replicates pooled in a cell stay separate subjects)
    groups = Dict{Tuple{String, Int}, Vector{eltype(cell)}}()
    order = Tuple{String, Int}[]
    for c in cell
        haskey(groups, (c[4], c[5])) || push!(order, (c[4], c[5]))
        push!(get!(groups, (c[4], c[5]), eltype(cell)[]), c)
    end
    vals, cens = Float64[], Bool[]
    for id in order
        g = groups[id]
        if subject_value === :mean
            push!(vals, sum(c[1] for c in g) / length(g))
            push!(cens, any(c -> c[3], g))
        else
            best = g[1]
            for c in g
                (subject_value === :last ? c[2] >= best[2] : c[2] < best[2]) && (best = c)
            end
            push!(vals, best[1])
            push!(cens, best[3])
        end
    end
    return vals, cens
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

# Quantile rows for items = (replicate, time, value, censored, id) tuples (see prepare_items).
# Replicates are pooled (replicate 0) unless by_replicate. Rows are ordered replicate, bin, q;
# the replicate field is present only when by_replicate; `below_lloq` and `limit` (the bound of
# a quantile below the limit, NaN otherwise) only for blq handling :censored. `stats`, when given, receives (n censored, n values) per (replicate, bin) cell.
function summarize_items(items, binning::Binning, qs, by_replicate, min_n, blq; max_failure_rate::Real = 0.0,
        weighting::Symbol = :observation, subject_value::Symbol = :mean, stats = nothing)
    cells = Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64, Bool, String, Int}}}()
    failed = Dict{Tuple{Int, Int}, Int}()
    reps = Set{Int}()
    for (rep, time, v, c, id) in items
        by_replicate && push!(reps, rep)
        i = bin_of(binning, time)
        i == 0 && continue
        key = (by_replicate ? rep : 0, i)
        if isfinite(v)
            push!(get!(cells, key, Tuple{Float64, Float64, Bool, String, Int}[]), (v, time, c, id, rep))
        else
            failed[key] = get(failed, key, 0) + 1
        end
    end
    repids = by_replicate ? sort!(collect(reps)) : [0]
    censored_mode = blq !== nothing && blq.handling === :censored
    out = NamedTuple[]
    for rep in repids, i in eachindex(binning.lo)
        cell = get(cells, (rep, i), Tuple{Float64, Float64, Bool, String, Int}[])
        vals, cens = collapse_cell(cell, weighting, subject_value)
        stats === nothing || (stats[(rep, i)] = (nblq = count(cens), ntot = length(vals)))
        qv, n, nblq, below, lims = cell_quantiles(vals, cens, qs, min_n, blq)
        nf = get(failed, (rep, i), 0)
        nfin = length(cell)
        rate = rate_of(nf, nfin)
        for (qi, qq) in enumerate(qs)
            row = (bin = i, bin_lo = binning.lo[i], bin_hi = binning.hi[i], t_mid = t_mid(binning, i), q = qq, value = qv[qi], n, n_failed = nf, n_blq = nblq,
                n_finite = nfin, failure_rate = rate, unreliable = is_unreliable(rate, max_failure_rate))
            censored_mode && (row = merge(row, (; below_lloq = below[qi], limit = lims[qi])))
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
                by_replicate = false, min_n = 1, blq = nothing, max_failure_rate = 0.0,
                weighting = :observation, subject_value = :mean) -> Vector{NamedTuple}

Quantiles of one endpoint of a simulation, over subjects, per time bin.

Estimand: by default (`weighting = :observation`) the values of a bin are pooled, so a subject
with several observations in a bin counts several times. This is the usual VPC and what
the observed-data percentiles of a standard VPC do. With `weighting = :subject` each subject
contributes one value per bin (and replicate), its `subject_value` = `:mean`, `:first` or
`:last` (by time) of its values in the bin; `n` is then the number of subjects. With
`bins = :obstimes` and a common sampling schedule the two are identical. A subject mean that
includes any censored value is itself censored (it is unknown, whatever the placeholder of the
censored value): with `:censored` it is ranked as censored below its own limit (the mean with the
censored values at their limits), with `:drop` the subject is left out of that bin, with
`:censor_to_lloq` and `:keep` its number is used and it counts in `n_blq`. With `:first` or
`:last` the subject value is censored when the picked value is.
`n_failed` and `n_finite` always count observations.

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
  - `blq = (lloq = x, handling = :censor_to_lloq | :drop | :keep | :censored)`: values below
    `lloq` are set to `lloq` (default), removed (this biases quantiles upwards), or kept;
    `n_blq` counts the values below `lloq` in every case, `n` the values that were used.
    `:censored` is the censored-VPC convention: the values below `lloq` are kept as censored,
    each known only to lie below its limit (`lloq`, so below every quantified value), and a
    quantile that depends on a censored order statistic (the lower of the two neighbouring
    order statistics of the type 7 interpolation is censored) is `NaN` with `below_lloq = true`
    (an extra field of the rows in this mode), not a number. When the limits differ (subject
    means, or a prediction-corrected VPC) an order statistic is a number only when it is the
    same wherever the censored values lie below their own limits
    (`QSPAtlas.censored_quantile_limits`). The rows then also have `limit`: for a quantile
    below the limit, the largest value its lower order statistic can take (`lloq` with a common
    limit), `NaN` otherwise.
  - Failure gate: `n_finite` is the number of finite values of the bin (before any BLQ
    handling), `failure_rate = n_failed / (n_failed + n_finite)` (`NaN` for a bin without
    values); [`failure_rate`](@ref) gives the overall rate. `max_failure_rate = 0` (the default)
    only reports. With `max_failure_rate > 0`, bins whose rate is above it get `unreliable = true`
    (and one warning); their `value` is not changed or removed, so a model decides what to do.
"""
function percentiles(r::SimResult, endpoint::Symbol; q = (0.1, 0.5, 0.9), bins = :obstimes, value::Symbol = :dv,
        by_replicate::Bool = false, min_n::Integer = 1, blq = nothing, max_failure_rate::Real = 0.0,
        weighting::Symbol = :observation, subject_value::Symbol = :mean)
    check_value(value)
    qs = check_quantiles(q)
    mfr = check_max_failure_rate(max_failure_rate)
    check_weighting(weighting, subject_value)
    blq = check_blq(blq)
    items, rows = endpoint_items(r, endpoint, value)
    pitems = prepare_items(items, [o.id for o in rows], blq)
    out = summarize_items(pitems, make_binning(bins, [o.time for o in rows]), qs, by_replicate, max(min_n, 1), blq; max_failure_rate = mfr, weighting, subject_value)
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
    `max_failure_rate` (never when that is 0, the report-only default). With
    `blq.handling = :censored` four fields are added: `below_lloq`, `lo_below`, `hi_below`
    (the median, `lo` or `hi` is below the limit: its value is `NaN`) and `n_cens` (replicates
    whose quantile was below the limit; they count in `n_rep`).
  - `observed`: `nothing`, or rows `(bin, bin_lo, bin_hi, t_mid, q, value, n, n_failed, n_blq)`
    of the observed quantiles (diagnostic only), plus `below_lloq` and `limit` with `:censored`
    (see [`percentiles`](@ref)).
  - `blq_fraction`: `nothing` unless `blq` was given, else one row per bin
    `(bin, bin_lo, bin_hi, t_mid, median, lo, hi, mc_se, n, n_rep, observed, n_observed, n_observed_blq)`:
    the fraction of values below `lloq`, simulated (median and `pi_level` interval across the replicates,
    `mc_se` of the median, `n` the median number of values per replicate) and observed
    (`NaN` without observed data).
  - `pred_bin`: for a prediction-corrected VPC the median typical population prediction of each
    bin (over the records of the simulated design), else `nothing`;
    `pc` is `:none`, `:pc` or `:pvc`; `n_nopred` counts values left out because their prediction
    was not positive and finite.
  - `weighting`, `subject_value` (see [`percentiles`](@ref)), `stratum` (the key of this result
    inside a [`StratifiedVPC`](@ref), else `nothing`).
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
    blq_fraction::Union{Nothing, Vector{NamedTuple}}
    pred_bin::Union{Nothing, Vector{Float64}}
    pc::Symbol
    weighting::Symbol
    subject_value::Symbol
    n_nopred::Int
    stratum::Any
end

VPCResult(rows, observed, endpoint, q, pi_level, replicates, seed, value, blq, n_failed, n_unbinned, failure_rate, max_failure_rate;
    blq_fraction = nothing, pred_bin = nothing, pc = :none, weighting = :observation, subject_value = :mean, n_nopred = 0, stratum = nothing) =
    VPCResult(rows, observed, endpoint, q, pi_level, replicates, seed, value, blq, n_failed, n_unbinned, failure_rate, max_failure_rate,
        blq_fraction, pred_bin, pc, weighting, subject_value, n_nopred, stratum)

VPCResult(rows, observed, endpoint, q, pi_level, replicates, seed, value, blq, n_failed, n_unbinned) =
    VPCResult(rows, observed, endpoint, q, pi_level, replicates, seed, value, blq, n_failed, n_unbinned, NaN, 0.0)

failure_rate(v::VPCResult) = v.failure_rate

Base.length(v::VPCResult) = length(v.rows)
Base.iterate(v::VPCResult, st...) = iterate(v.rows, st...)
Base.getindex(v::VPCResult, i) = v.rows[i]
function Base.show(io::IO, v::VPCResult)
    print(io, "VPCResult(", v.endpoint, ", ", length(unique(r.bin for r in v.rows)), " bins, q = ", Tuple(v.q), ", ",
        v.replicates, " replicates", v.observed === nothing ? "" : ", with observed", v.n_failed == 0 ? "" : ", $(v.n_failed) values excluded",
        v.pc === :none ? "" : ", prediction-corrected", v.stratum === nothing ? "" : ", stratum $(repr(v.stratum))", ")")
end

"""
    StratifiedVPC

Result of [`vpc`](@ref) with `strata`: one [`VPCResult`](@ref) per stratum, in `results`, with the
stratum keys in `strata` (same order; every result also has its key in `stratum`). Indexing by a
key gives its result, indexing by an integer the i-th; iterating gives `key => result` pairs;
`length` is the number of strata. [`plot_vpc`](@ref) draws one panel per stratum. Use
[`mc_check`](@ref) and [`failure_rate`](@ref) on the individual results.
"""
struct StratifiedVPC
    strata::Vector{Any}
    results::Vector{VPCResult}
end

Base.length(s::StratifiedVPC) = length(s.results)
Base.iterate(s::StratifiedVPC, st = 1) = st > length(s) ? nothing : (s.strata[st] => s.results[st], st + 1)
function Base.getindex(s::StratifiedVPC, key)
    i = findfirst(isequal(key), s.strata)
    i === nothing && (key isa Integer && 1 <= key <= length(s) ? (i = key) : throw(KeyError(key)))
    return s.results[i]
end
Base.keys(s::StratifiedVPC) = s.strata
Base.show(io::IO, s::StratifiedVPC) = print(io, "StratifiedVPC(", length(s), " strata: ", join(repr.(s.strata), ", "), ")")

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
        max_failure_rate = 0.0, weighting = :observation, subject_value = :mean,
        pc = false, variability_correction = false, strata = nothing,
        simulate_kwargs...) -> VPCResult or StratifiedVPC
    vpc(r::SimResult; endpoint, ..., pred = nothing)

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

# Estimand: what is pooled

The binned percentiles pool observations. With `weighting = :observation` (the default) every
value in a bin counts once, so a subject with more observations in a bin weighs more. That is
the usual VPC (observed and simulated alike) and the right reading when the bins are narrow
or the sampling is the same for everyone; with wide bins and unequal sampling (sparse
or informative visit schedules) the percentile is of the pooled observations, not of subjects.
`weighting = :subject` collapses each subject to one value per bin and replicate, its
`subject_value` (`:mean`, `:first` or `:last` by time), so every subject weighs once (`n` counts
subjects). Both observed and simulated values are treated the same way, so the comparison is
like with like. With one observation per subject and bin the two are identical.

# Censored VPC (blq)

Giving `blq = (lloq = x, handling = ...)` also computes, per bin, the fraction of values below
`lloq` (`VPCResult.blq_fraction`): observed (pooled over observed values, or one per subject
with `weighting = :subject`) and simulated (the fraction of each replicate, summarised by its
median and the `pi_level` interval across replicates). With `handling = :censored` the
quantiles follow the standard censored VPC: values below `lloq` are kept as censored, below
every quantified value; a quantile whose interpolation touches a censored order statistic
is not a number but `NaN` with `below_lloq = true`. The replicate quantiles are then
summarised the same way (a replicate quantile below `lloq` is the lowest value): the
median, `lo` or `hi` over replicates is `NaN` with its flag (`below_lloq`, `lo_below`,
`hi_below`) when its rank falls among the censored replicates, and `mc_se` is computed with the
censored replicate quantiles set to `lloq` (`NaN` when the median itself is below `lloq`).
Observed values are censored when `dv < lloq`, or, if the observed rows have a `blq` field,
when it is true (use it when below-limit records are stored at the limit). Censored values are
flagged on the raw scale also in a prediction-corrected VPC; the correction then moves each
record's limit with its value, to `LLOQ * PRED_bin / PRED_ij` (and through the variability
correction), and a censored value is ranked only against its own limit, so a quantile is a
number only when it does not depend on where any censored value lies below its limit (see
[`percentiles`](@ref); the same rule applies to subject means with `weighting = :subject`).
The replicate quantiles below the limit are then summarised with their own bounds. M3 is an
estimation method, not an option of the simulation, and is not part of this function.

# Prediction-corrected VPC (pc = true)

Bergstrand et al. (2011, AAPS J 13:143) remove the part of the variability that comes from
different doses, covariates or sampling times. Every observed and simulated value is
normalised by the typical population prediction of its record (η and κ zero, no residual
error: [`zero_re_prediction`](@ref)) relative to the median prediction of its bin:

    pcY_ij = Y_ij * PRED_bin / PRED_ij

`PRED_ij` is the typical prediction of record j of subject i, `PRED_bin` the median of the
predictions of the simulated design (each record once) in the bin of the observation;
observed and simulated values use the same `PRED_bin`. Records are matched by (id, time, k),
k numbering the records of a subject at one time in their order (`:pre` before `:post`), so a
pre-dose and a post-dose record at a dose time keep their own predictions; observed rows at
one time are numbered the same way, `order = :pre` rows first when they have an `order` field
(as from [`read_nmdata`](@ref)), otherwise in the order given. `vpc(model, subjects; pc = true)`
computes the predictions itself with η = 0 and κ = 0 whatever `η` and `κ` are passed for the
simulation (observed records that are not records of the cohort design are predicted from
their subject observed at those records). `vpc(r::SimResult; pc = true, pred = ...)` needs
`pred`: a `SimResult` from `simulate(model, subjects; η = :zero, κ = :zero, residual = false)`
(checked: nonzero η or κ, or residual error, is an error), or a dictionary keyed by
`(id, time, k)` or `(id, time)` (applies to every record of the subject at that time, so it
cannot tell a pre-dose from a post-dose record); observed rows may carry their own `pred`
field. The correction
is proportional, so it needs positive predictions: values with a prediction that is not positive
and finite (a concentration of zero at the first time, say) cannot be corrected, are left
out and are counted in `n_nopred` (one warning). The lower-bound term and the log-scale
variant of the paper are not implemented. Under `weighting = :subject` the pooled `PRED_bin` is
still the median over observations.

`variability_correction = true` adds the paper's variability correction (pvcVPC): the spread of
the corrected values around `PRED_bin` is equalised across the observations of a bin,

    pvcY_ij = PRED_bin + (pcY_ij - PRED_bin) * sd_bin / sd_ij

with `sd_ij` the standard deviation over the replicates of the simulated `pcY_ij` of that record
and `sd_bin` the median of the `sd_ij` of the bin (finite and positive ones). It needs at
least 2 replicates, and observed records must be in the simulated design (same id, time and k);
an observation with no spread (`sd_ij` not positive) is left uncorrected. `pred_bin` in the result holds `PRED_bin`.

# Stratification

`strata = subject -> key` splits the cohort and the observed data by the key of each subject
(the function gets the [`Subject`](@ref), or for `vpc(r::SimResult)` the subject row, both with
`id` and `covariates`; a dictionary `id => key` also works) and returns a [`StratifiedVPC`](@ref), one
independent VPC per stratum (also the `PRED_bin` of a pcVPC and the bins' `n_rep`). One simulation serves all strata, so the
replicate and seed are the same as without strata. Observed rows must belong to a subject of the cohort (matched by `id`).
"""
function vpc(model::PopulationModel, subjects::AbstractVector{<:Subject}; endpoint::Symbol, q = (0.1, 0.5, 0.9), bins = :obstimes,
        replicates::Integer = 200, pi_level::Real = 0.95, seed::Integer = 0, value::Symbol = :dv, blq = nothing,
        observed = nothing, min_n::Integer = 1, max_failure_rate::Real = 0.0, weighting::Symbol = :observation,
        subject_value::Symbol = :mean, pc::Bool = false, variability_correction::Bool = false, strata = nothing, simulate_kwargs...)
    check_max_failure_rate(max_failure_rate)
    replicates >= 1 || throw(ArgumentError("replicates must be >= 1"))
    check_vpc_options(; q, bins, pi_level, value, blq, weighting, subject_value, pc, variability_correction, replicates)   # fail before the expensive simulation
    res = simulate(model, subjects; replicates, seed, simulate_kwargs...)
    pred, observed = pc ? model_pred_table(model, subjects, endpoint, observed; simulate_kwargs...) : (nothing, observed)
    strata_arg = strata isa Function ? Dict{String, Any}(s.id => strata(s) for s in subjects) : strata
    return vpc(res; endpoint, q, bins, pi_level, seed, value, blq, observed, min_n, max_failure_rate,
        weighting, subject_value, pc, variability_correction, pred, strata = strata_arg)
end

function check_vpc_options(; q, bins, pi_level, value, blq, weighting = :observation, subject_value = :mean, pc = false,
        variability_correction = false, replicates = 2)
    check_value(value)
    check_quantiles(q)
    check_blq(blq)
    check_weighting(weighting, subject_value)
    (0 < pi_level < 1) || throw(ArgumentError("pi_level must be in (0, 1), got $pi_level"))
    bins === :obstimes || make_binning(bins, Float64[])
    (variability_correction && !pc) && throw(ArgumentError("variability_correction needs pc = true"))
    (variability_correction && replicates < 2) && throw(ArgumentError("variability_correction needs at least 2 replicates (the spread of each observation across replicates)"))
    return nothing
end

function vpc(r::SimResult; endpoint::Symbol, q = (0.1, 0.5, 0.9), bins = :obstimes, pi_level::Real = 0.95,
        seed = nothing, value::Symbol = :dv, blq = nothing, observed = nothing, min_n::Integer = 1, max_failure_rate::Real = 0.0,
        weighting::Symbol = :observation, subject_value::Symbol = :mean, pc::Bool = false, variability_correction::Bool = false,
        pred = nothing, strata = nothing)
    nrep = length(unique(s.replicate for s in r.subjects))
    check_vpc_options(; q, bins, pi_level, value, blq, weighting, subject_value, pc, variability_correction, replicates = max(nrep, 2))
    variability_correction && nrep < 2 && throw(ArgumentError("variability_correction needs at least 2 replicates, the result has $nrep"))
    mfr = check_max_failure_rate(max_failure_rate)
    opts = (; endpoint, q, bins, pi_level, seed, value, blq, min_n, mfr, weighting, subject_value, pc, variability_correction,
        pred = pc ? pred_table(pred, endpoint) : nothing)
    strata === nothing && return vpc_single(r, observed; opts..., stratum = nothing)
    idkey = Dict{String, Any}()
    for srow in r.subjects
        haskey(idkey, srow.id) || (idkey[srow.id] = strata isa AbstractDict ? strata[srow.id] : strata(srow))
    end
    keys_ = unique([idkey[srow.id] for srow in r.subjects])      # order of first appearance
    try
        keys_ = sort(keys_)
    catch
    end
    results = VPCResult[]
    for key in keys_
        sobs = filter(o -> isequal(idkey[o.id], key), r.obs)
        ssub = filter(s -> isequal(idkey[s.id], key), r.subjects)
        oobs = nothing
        if observed !== nothing
            for o in observed
                haskey(idkey, string(get(o, :id, ""))) || throw(ArgumentError("observed id $(repr(get(o, :id, ""))) is not among the subjects, so its stratum is unknown"))
            end
            oobs = [o for o in observed if isequal(idkey[string(o.id)], key)]
        end
        push!(results, vpc_single(SimResult(sobs, ssub, r.settings), oobs; opts..., stratum = key))
    end
    return StratifiedVPC(Any[k for k in keys_], results)
end

# ---- records and their predictions -------------------------------------------------------

# A record of one endpoint is (id, time, k): k numbers the records of that subject at that time
# in their order (a pre-dose record before the post-dose one, as `Subject` stores them), so a
# pre-dose and a post-dose record at a dose time are two records with their own prediction.
const RecKey = Tuple{String, Float64, Int}

# Record keys of simulated rows of one endpoint: k counts within (replicate, id, time).
function sim_record_keys(rows)
    seen = Dict{Tuple{Int, String, Float64}, Int}()
    out = Vector{RecKey}(undef, length(rows))
    for (i, o) in enumerate(rows)
        g = (o.replicate, o.id, Float64(o.time))
        k = seen[g] = get(seen, g, 0) + 1
        out[i] = (o.id, Float64(o.time), k)
    end
    return out
end

# Record keys of observed rows: k counts within (id, time), `order = :pre` records first when the
# rows have an `order` field (as `read_nmdata` gives them), otherwise in the order given. Rows
# whose `endpoint` field names another endpoint are numbered on their own.
function obs_record_keys(observed, endpoint::Symbol)
    n = length(observed)
    keys_ = Vector{RecKey}(undef, n)
    groups = Dict{Tuple{String, Float64, Bool}, Vector{Int}}()
    for (i, o) in enumerate(observed)
        ep = hasproperty(o, :endpoint) ? o.endpoint : nothing
        other = !(ep === nothing || ep === endpoint)
        push!(get!(groups, (string(get(o, :id, "")), Float64(o.time), other), Int[]), i)
    end
    for ((id, tm, _), idx) in groups
        pre(i) = hasproperty(observed[i], :order) && observed[i].order === :pre ? 0 : 1
        for (k, i) in enumerate(sort(idx; by = i -> (pre(i), i)))
            keys_[i] = (id, tm, k)
        end
    end
    return keys_
end

# Predictions by record: `rec` (id, time, k) => prediction, and `bytime` (id, time) => prediction,
# which applies to every record of the subject at that time (a dictionary keyed by (id, time)).
struct PredTable
    rec::Dict{RecKey, Float64}
    bytime::Dict{Tuple{String, Float64}, Float64}
end

pred_lookup(t::PredTable, key::RecKey) = get(t.rec, key, get(t.bytime, (key[1], key[2]), nothing))

# A SimResult of predictions must be the typical population prediction: η and κ zero, no residual error.
function check_typical(pred::SimResult)
    iszero_nt(x) = x isa NamedTuple ? all(v -> v == 0, values(x)) : true
    for s in pred.subjects
        ok = iszero_nt(get(s, :η, (;))) && (!haskey(s, :κ) || all(iszero_nt, values(s.κ)))
        ok || throw(ArgumentError("pred must be the typical population prediction, with η = 0 and κ = 0: simulate(model, subjects; η = :zero, κ = :zero, residual = false); subject $(repr(s.id)) has nonzero random effects"))
    end
    get(pred.settings, :residual, false) === true &&
        throw(ArgumentError("pred must be simulated without residual error (residual = false)"))
    return nothing
end

# record => zero-random-effect prediction of one endpoint
function pred_table(pred, endpoint::Symbol)
    pred === nothing && return nothing
    pred isa PredTable && return pred
    if pred isa SimResult
        check_typical(pred)
        rows = [o for o in pred.obs if o.endpoint === endpoint]
        isempty(rows) && return PredTable(Dict{RecKey, Float64}(), Dict{Tuple{String, Float64}, Float64}())
        rep1 = minimum(o.replicate for o in rows)
        rows = [o for o in rows if o.replicate == rep1]
        return PredTable(Dict{RecKey, Float64}(k => Float64(o.ipred) for (k, o) in zip(sim_record_keys(rows), rows)), Dict{Tuple{String, Float64}, Float64}())
    elseif pred isa AbstractDict
        rec = Dict{RecKey, Float64}()
        bytime = Dict{Tuple{String, Float64}, Float64}()
        for (k, v) in pred
            if k isa Tuple && length(k) == 3
                rec[(string(k[1]), Float64(k[2]), Int(k[3]))] = Float64(v)
            elseif k isa Tuple && length(k) == 2
                bytime[(string(k[1]), Float64(k[2]))] = Float64(v)
            else
                throw(ArgumentError("pred keys must be (id, time) or (id, time, k), got $(repr(k))"))
            end
        end
        return PredTable(rec, bytime)
    end
    throw(ArgumentError("pred must be a SimResult from simulate(model, subjects; η = :zero, κ = :zero, residual = false) or a Dict (id, time) or (id, time, k) => prediction, got $(typeof(pred))"))
end

# Typical population prediction (η = 0, κ = 0, no residual error), Bergstrand et al. 2011: the
# simulate keywords that set random effects or replicates are replaced.
typical_simulate(model, subjects, kw) =
    simulate(model, subjects; kw..., η = :zero, κ = :zero, residual = false, replicates = 1, seed = 0)

# The predictions for vpc(model, ...; pc = true): the table of the cohort design, and the observed
# rows with a `pred` field. An observed record that is not in the design (its id, time and k) is
# predicted from its subject observed at the observed records of those times (with their
# :pre/:post order), so a record is never given the prediction of another record.
function model_pred_table(model, subjects, endpoint, observed; kwargs...)
    kw = filter(kv -> kv.first ∉ (:η, :κ, :residual, :replicates, :seed), kwargs)
    tbl = pred_table(typical_simulate(model, subjects, kw), endpoint)
    observed === nothing && return tbl, observed
    okeys = obs_record_keys(observed, endpoint)
    own(o) = !hasproperty(o, :endpoint) || o.endpoint === nothing || o.endpoint === endpoint
    need = Dict{String, Set{Float64}}()
    for (o, key) in zip(observed, okeys)
        (hasproperty(o, :pred) || !own(o)) && continue
        haskey(tbl.rec, key) || push!(get!(need, key[1], Set{Float64}()), key[2])
    end
    extra_tbl = Dict{RecKey, Float64}()
    if !isempty(need)
        byid = Dict(s.id => s for s in subjects)
        extra = Subject[]
        for id in sort!(collect(keys(need)))
            haskey(byid, id) || throw(ArgumentError("observed id $(repr(id)) is not among the subjects: give its rows a `pred` field"))
            sel = [o for o in observed if own(o) && string(get(o, :id, "")) == id && Float64(o.time) in need[id]]
            s = byid[id]
            push!(extra, Subject(s.id; regimen = s.regimen, covariates = s.covariates, control = s.control,
                obstimes = [Float64(o.time) for o in sel],
                obsorder = [hasproperty(o, :order) && o.order === :pre ? :pre : :post for o in sel],
                timevarying = [t => v for (t, v) in zip(s.covariate_times, s.covariate_values)], covariate_rule = s.covariate_rule,
                occasions = [t => i for (t, i) in zip(s.occasion_times, s.occasion_ids)]))
        end
        extra_tbl = pred_table(typical_simulate(model, extra, kw), endpoint).rec
    end
    withpred = map(zip(observed, okeys)) do (o, key)
        (hasproperty(o, :pred) || !own(o)) && return o
        p = (key[1] in keys(need) && key[2] in need[key[1]]) ? get(extra_tbl, key, NaN) : tbl.rec[key]
        return (; o..., pred = p)
    end
    return tbl, withpred
end

# ---- prediction correction ---------------------------------------------------------------

# Apply pcY = Y PRED_bin / PRED (and, optionally, the variability correction) to simulated items
# and observed items, matched to their predictions by record (`skeys`, `okeys`). `opred` are the
# predictions of the observed items. `sspread` are the simulated values whose corrected spread
# gives sd_ij (the item values, except for censored items under :censored, whose slot holds the
# limit: the simulated value itself is used there). A censored item under :censored carries its
# limit in the value slot, so the same transformation gives its limit on the corrected scale,
# LLOQ PRED_bin / PRED_ij (and its pvc transform). Returns (simulated items, observed items,
# PRED_bin per bin, n left out).
function pc_apply(sitems, skeys, sspread, oitems, okeys, opred, binning::Binning, table, variability::Bool)
    table === nothing && throw(ArgumentError("pc = true needs the typical population predictions: pass pred = simulate(model, subjects; η = :zero, κ = :zero, residual = false) (vpc(model, subjects; ...) does it itself)"))
    valid(p) = isfinite(p) && p > 0
    nb = length(binning.lo)
    sp = Vector{Float64}(undef, length(sitems))
    design = Dict{RecKey, Float64}()                      # each record of the design once
    for (k, it) in enumerate(sitems)
        p = pred_lookup(table, skeys[k])
        p === nothing && throw(ArgumentError("no prediction for subject $(repr(it[5])) at time $(it[2]) (record $(skeys[k][3]) at that time) in `pred`"))
        sp[k] = p
        valid(p) && bin_of(binning, it[2]) > 0 && (design[skeys[k]] = p)
    end
    binpreds = [Float64[] for _ in 1:nb]
    for ((_, tm, _), p) in design
        push!(binpreds[bin_of(binning, tm)], p)
    end
    pb = Float64[isempty(v) ? NaN : median(v) for v in binpreds]
    nout = 0
    # corrected items, dropping those that cannot be corrected; returns the kept items, their keys and spread values
    function correct(items, keys_, preds, spread)
        out = Tuple{Int, Float64, Float64, Bool, String}[]
        kout = RecKey[]
        sout = Float64[]
        for (k, it) in enumerate(items)
            i = bin_of(binning, it[2])
            i == 0 && continue
            if !(valid(preds[k]) && isfinite(pb[i]))
                nout += 1
                continue
            end
            f = pb[i] / preds[k]
            push!(out, (it[1], it[2], it[3] * f, it[4], it[5]))
            push!(kout, keys_[k])
            spread === nothing || push!(sout, spread[k] * f)
        end
        return out, kout, sout
    end
    spc, skc, ssc = correct(sitems, skeys, sp, sspread)
    opc, okc = oitems === nothing ? (nothing, nothing) : correct(oitems, okeys, opred, nothing)
    if variability
        groups = Dict{RecKey, Vector{Float64}}()
        for (it, key, sv) in zip(spc, skc, ssc)
            isfinite(sv) && push!(get!(groups, key, Float64[]), sv)
        end
        sdij = Dict(k => std(v) for (k, v) in groups if length(v) >= 2)
        binsd = fill(NaN, nb)
        persd = [Float64[] for _ in 1:nb]
        for ((_, tm, _), sdv) in sdij
            (isfinite(sdv) && sdv > 0) && push!(persd[bin_of(binning, tm)], sdv)
        end
        for i in 1:nb
            isempty(persd[i]) || (binsd[i] = median(persd[i]))
        end
        function vc_item(it, key, strict)
            i = bin_of(binning, it[2])
            sd = get(sdij, key, nothing)
            if sd === nothing
                strict && throw(ArgumentError("variability_correction needs every observed record in the simulated design: no simulated record (id, time, k) = ($(repr(key[1])), $(key[2]), $(key[3]))"))
                return it
            end
            f = (isfinite(sd) && sd > 0 && isfinite(binsd[i])) ? binsd[i] / sd : 1.0
            return (it[1], it[2], pb[i] + (it[3] - pb[i]) * f, it[4], it[5])
        end
        spc = [vc_item(it, key, false) for (it, key) in zip(spc, skc)]
        opc === nothing || (opc = [vc_item(it, key, true) for (it, key) in zip(opc, okc)])
    end
    return spc, opc, pb, nout
end

# ---- one VPC ------------------------------------------------------------------------------

function vpc_single(r::SimResult, observed; endpoint, q, bins, pi_level, seed, value, blq, min_n, mfr, weighting, subject_value,
        pc, variability_correction, pred, stratum)
    qs = check_quantiles(q)
    blq = check_blq(blq)
    min_n = max(min_n, 1)
    items, rows = endpoint_items(r, endpoint, value)
    binning = make_binning(bins, [o.time for o in rows])
    sitems = prepare_items(items, [o.id for o in rows], blq)
    oitems, opred, okeys, unbinned = nothing, nothing, nothing, 0
    if observed !== nothing
        oraw = [(0, Float64(o.time), Float64(o.dv)) for o in observed]
        flags = blq === nothing ? nothing : Bool[hasproperty(o, :blq) ? Bool(o.blq) : Float64(o.dv) < blq.lloq for o in observed]
        oitems = prepare_items(oraw, [string(get(o, :id, "")) for o in observed], blq; flags)
        unbinned = count(o -> bin_of(binning, o[2]) == 0, oraw)
        if pc
            okeys = obs_record_keys(observed, endpoint)
            opred = Float64[]
            for (o, key) in zip(observed, okeys)
                p = hasproperty(o, :pred) ? Float64(o.pred) : (pred === nothing ? nothing : pred_lookup(pred, key))
                push!(opred, p === nothing ? NaN : p)
            end
        end
    end
    pred_bin, nnopred = nothing, 0
    if pc
        # the spread of a censored simulated record (:censored keeps its limit in the slot) uses its simulated value
        sspread = [blq !== nothing && blq.handling === :censored && it[4] ? Float64(raw[3]) : it[3] for (it, raw) in zip(sitems, items)]
        sitems, oitems, pred_bin, nnopred = pc_apply(sitems, sim_record_keys(rows), sspread, oitems, okeys, opred, binning, pred, variability_correction)
        nnopred > 0 && @warn "vpc: $nnopred values have no positive finite prediction and could not be prediction-corrected; they are left out (n_nopred)"
    end
    repstats = Dict{Tuple{Int, Int}, NamedTuple}()
    reps = summarize_items(sitems, binning, qs, true, min_n, blq; weighting, subject_value, stats = repstats)
    # gather the replicate quantiles per (bin, q)
    cells = Dict{Tuple{Int, Float64}, Vector{NamedTuple}}()
    for row in reps
        push!(get!(cells, (row.bin, row.q), NamedTuple[]), row)
    end
    plo, phi = (1 - pi_level) / 2, (1 + pi_level) / 2
    censored_mode = blq !== nothing && blq.handling === :censored
    out = NamedTuple[]
    ndegenerate = 0
    for i in eachindex(binning.lo), qq in qs
        cell = cells[(i, qq)]
        extra = (;)
        if censored_mode
            # a replicate quantile below the limit is known only to be below its bound (`limit`:
            # lloq without prediction correction or subject means, its own bound otherwise)
            nums = sort!(Float64[c.value for c in cell if isfinite(c.value)])
            lims = Float64[c.limit for c in cell if c.below_lloq]
            ncens = length(lims)
            md, mdb, _ = censored_quantile_limits(nums, lims, 0.5)
            lo, lob, _ = censored_quantile_limits(nums, lims, plo)
            hi, hib, _ = censored_quantile_limits(nums, lims, phi)
            nr = length(nums) + ncens
            se, degenerate = mdb ? (NaN, false) : mc_se_median_checked([nums; lims], stream_rng(seed === nothing ? 0 : seed, i, string(qq), :mcse))
            extra = (below_lloq = mdb, lo_below = lob, hi_below = hib, n_cens = ncens)
        else
            vals = Float64[c.value for c in cell if isfinite(c.value)]
            md, lo, hi = isempty(vals) ? (NaN, NaN, NaN) : (median(vals), quantile(vals, plo), quantile(vals, phi))
            se, degenerate = mc_se_median_checked(vals, stream_rng(seed === nothing ? 0 : seed, i, string(qq), :mcse))
            nr = length(vals)
        end
        ndegenerate += degenerate
        nf, nfin = sum(c.n_failed for c in cell), sum(c.n_finite for c in cell)
        rate = rate_of(nf, nfin)
        push!(out, merge((bin = i, bin_lo = binning.lo[i], bin_hi = binning.hi[i], t_mid = t_mid(binning, i), q = qq,
            median = md, lo, hi, mc_se = se, n = round(Int, median([c.n for c in cell])),
            n_rep = nr, n_failed = nf, n_blq = sum(c.n_blq for c in cell),
            failure_rate = rate, unreliable = is_unreliable(rate, mfr)), extra))
    end
    ndegenerate > 0 && @warn "vpc: mc_se is NaN in $ndegenerate (bin, q) cell(s) where the replicate quantiles are heavily tied (the bootstrap medians hardly vary); it is not estimable there, not 0"
    nfail = sum(row.n_failed for row in out if row.q == qs[1]; init = 0)
    overall = failure_rate([(; q = row.q, n_failed = row.n_failed, n_finite = sum(c.n_finite for c in cells[(row.bin, row.q)])) for row in out])
    nfail > 0 && @warn "vpc: $nfail non-finite simulated values (failed solves or NaN) were excluded; they are counted in n_failed"
    warn_unreliable(out, mfr, "vpc")
    obs_rows, ostats = nothing, Dict{Tuple{Int, Int}, NamedTuple}()
    if observed !== nothing
        obs_rows = summarize_items(oitems, binning, qs, false, min_n, blq; weighting, subject_value, stats = ostats)
    end
    blqfrac = blq === nothing ? nothing : blq_fraction_rows(binning, repstats, ostats, plo, phi, seed)
    return VPCResult(out, obs_rows, endpoint, qs, Float64(pi_level), nrep_of(r), seed === nothing ? nothing : Int(seed), value, blq, nfail, unbinned, overall, mfr;
        blq_fraction = blqfrac, pred_bin, pc = pc ? (variability_correction ? :pvc : :pc) : :none, weighting, subject_value,
        n_nopred = nnopred, stratum)
end

nrep_of(r::SimResult) = length(unique(s.replicate for s in r.subjects))

# fraction of values below the limit per bin: observed, and simulated across the replicates
function blq_fraction_rows(binning::Binning, repstats, ostats, plo, phi, seed)
    out = NamedTuple[]
    for i in eachindex(binning.lo)
        fr = Float64[]
        ns = Int[]
        for ((_, bi), st) in repstats
            bi == i || continue
            push!(ns, st.ntot)
            st.ntot > 0 && push!(fr, st.nblq / st.ntot)
        end
        md, lo, hi = isempty(fr) ? (NaN, NaN, NaN) : (median(fr), quantile(fr, plo), quantile(fr, phi))
        se, _ = mc_se_median_checked(sort(fr), stream_rng(seed === nothing ? 0 : seed, i, "blq_fraction", :mcse))
        ost = get(ostats, (0, i), (nblq = 0, ntot = 0))
        push!(out, (bin = i, bin_lo = binning.lo[i], bin_hi = binning.hi[i], t_mid = t_mid(binning, i), median = md, lo, hi, mc_se = se,
            n = isempty(ns) ? 0 : round(Int, median(ns)), n_rep = length(fr),
            observed = ost.ntot > 0 ? ost.nblq / ost.ntot : NaN, n_observed = ost.ntot, n_observed_blq = ost.nblq))
    end
    return out
end
