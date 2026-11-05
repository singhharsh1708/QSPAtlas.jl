"""
    Tolerance(; rms = 0.015, p95 = 0.03)

Pass rule for one curve. Errors are distances in the figure's own axis box, where
both axes are scaled to `[0, 1]` (log axes in decades). `rms` bounds the root mean
square distance over all reference points and `p95` bounds the 95th percentile, so
one badly traced point cannot fail a curve on its own. The maximum is still reported.

The defaults allow for digitising error: the plot panels digitised so far are 85 to
294 px wide or high, so one pixel is 0.3 to 1.2 % of an axis, and the tolerance is
roughly 1.3 to 3.5 px vertically. It is a judgement, not yet a measured error
distribution (see docs/src/method.md and the sensitivity tables). Tolerances are frozen in each
model's `protocol.jl` before its simulations are compared.
"""
Base.@kwdef struct Tolerance
    rms::Float64 = 0.015
    p95::Float64 = 0.03
end

"""
    CurveComparison

Result of comparing one simulated curve with one reference curve: per-point distances
`d` (axis-box units), their `rms`, `p95` and `max`, the `tolerance` used and whether
it `passed`.
"""
struct CurveComparison
    id::String
    d::Vector{Float64}
    rms::Float64
    p95::Float64
    max::Float64
    tolerance::Tolerance
    passed::Bool
end

passes(c::CurveComparison) = c.passed

# A simulated value can be placed on an axis if it is finite and, on a log axis,
# positive. Anything else (NaN, Inf, zero or negative on a log axis) is a gap: it is
# never moved onto the axis edge, which could otherwise create a false match.
usable(a::Axis, v) = isfinite(v) && (a.scale === :linear || v > 0)

# some values are finite but cannot be placed (zero or negative on a log axis): the
# simulation lies below the axis there
below_axis(ref, xs, ys) = any(isfinite(x) && isfinite(y) && !(usable(ref.xaxis, x) && usable(ref.yaxis, y)) for (x, y) in zip(xs, ys))

function all_missed(ref, tolerance)
    d = fill(1.0, length(ref.x))
    return CurveComparison(ref.id, d, 1.0, 1.0, 1.0, tolerance, false)
end

# distance from point p to segment ab in 2-D
function segdist(p, a, b)
    ab = b .- a
    L2 = ab[1]^2 + ab[2]^2
    L2 == 0 && return norm(p .- a)
    t = clamp(((p[1] - a[1]) * ab[1] + (p[2] - a[2]) * ab[2]) / L2, 0.0, 1.0)
    return norm(p .- (a .+ t .* ab))
end

"""
    compare_curve(ref, xs, ys; tolerance = Tolerance()) -> CurveComparison

Compare a simulated curve `(xs, ys)` against `ref`.

For `kind = :curve` references (lines drawn by the authors' own simulation), each
reference point's error is its shortest distance to the simulated polyline in the
normalised axis box. Using distance-to-curve instead of vertical error keeps steep
sections (dose spikes, rapid expansion) from being punished for tiny horizontal
digitising error.

For `kind = :points` references (observed data), the same distance is used; the
tolerance should then be loosened by the caller, since a typical-value curve is not
expected to pass through every observation.

For `kind = :bars` references (bar charts, tornado plots, single reported values),
x is a category position, not a continuous coordinate. The simulation must supply a
value at exactly each reference x (matched within 1e-9), and the error is the vertical
distance only, so a bar can never be matched by a neighbouring bar. A missing or
non-finite simulated bar counts as the maximum error (1.0), never as a pass.

For `kind = :band` references (a shaded interval, e.g. the authors' simulated percentile
band) and for `kind = :curve` references with the header `# metric: vertical` (stored as
`meta["metric"]`; percentile lines of a population simulation) the comparison is VERTICAL
at the reference x, not by nearest-polyline distance. The simulated curve is evaluated at
each reference x: a finite simulated point at exactly that x (within 1e-9 of the normalised
axis) is used as it is, otherwise the value is interpolated between the two neighbouring
finite points, linearly in the displayed coordinates (decades on a log axis, x and y). Only
then is the error taken: for a band 0 inside `[ylo, yhi]` (edges inside) and else the
distance to the nearest edge, for a vertical line the distance to the reference y, in
normalised axis units, capped at 1. A reference x outside the simulated range, or in a gap
of non-finite simulated points (gaps are never bridged), counts as the maximum error 1.0.
Simulated x must be sorted; one finite point is enough when it matches. For a band the
default [`Tolerance`](@ref) is the allowance for digitising its edges.

Each point's error is capped at 1.0 (one full axis box), so a simulation that leaves
the figure is scored as far off rather than infinitely off; the cap is applied to the
distance, never to the simulated points.
"""
function compare_curve(ref::ReferenceCurve, xs::AbstractVector, ys::AbstractVector; tolerance::Tolerance = Tolerance())
    length(xs) == length(ys) || throw(DimensionMismatch("xs and ys differ in length"))
    ref.kind === :bars && return compare_bars(ref, xs, ys, tolerance)
    ref.kind === :band && return compare_vertical(ref, xs, ys, tolerance, (ref.ylo, ref.yhi))
    ref.kind === :curve && get(ref.meta, "metric", "nearest") == "vertical" && return compare_vertical(ref, xs, ys, tolerance, (ref.y, ref.y))
    length(xs) >= 2 || throw(ArgumentError("need at least two simulated points"))
    # Non-finite points (a failed or blown-up solve) leave gaps: segments only join
    # neighbouring finite points, so a gap is never bridged by a straight line.
    # Distances are measured on the true (unclipped) geometry and only then capped at
    # 1.0: clipping the simulated points first would bend a steep segment towards the
    # reference and could report a far-off curve as a perfect match.
    finite = [usable(ref.xaxis, x) && usable(ref.yaxis, y) for (x, y) in zip(xs, ys)]
    sim = [(normalise(ref.xaxis, x), normalise(ref.yaxis, y)) for (x, y) in zip(xs, ys)]
    segs = [i for i in 1:(length(sim) - 1) if finite[i] && finite[i + 1]]
    if isempty(segs)
        # a simulation that lies entirely below a log axis (zero or negative values) is a
        # complete miss, not an error; one with no finite values at all is an error
        below_axis(ref, xs, ys) && return all_missed(ref, tolerance)
        throw(ArgumentError("simulation for $(ref.id) has no two neighbouring finite points"))
    end
    d = map(zip(ref.x, ref.y)) do (x, y)
        p = (normalise(ref.xaxis, x), normalise(ref.yaxis, y))
        min(minimum(segdist(p, sim[i], sim[i + 1]) for i in segs), 1.0)
    end
    r = sqrt(mean(abs2, d))
    q = quantile(d, 0.95)
    return CurveComparison(ref.id, d, r, q, maximum(d), tolerance, r <= tolerance.rms && q <= tolerance.p95)
end

function compare_bars(ref::ReferenceCurve, xs, ys, tolerance::Tolerance)
    d = map(zip(ref.x, ref.y)) do (x, y)
        i = findfirst(v -> isapprox(v, x; atol = 1e-9), xs)
        (i === nothing || !isfinite(ys[i])) && return 1.0
        return min(abs(normalise(ref.yaxis, ys[i]) - normalise(ref.yaxis, y)), 1.0)
    end
    r = sqrt(mean(abs2, d))
    q = quantile(d, 0.95)
    return CurveComparison(ref.id, d, r, q, maximum(d), tolerance, r <= tolerance.rms && q <= tolerance.p95)
end

# Vertical error at each reference x against the interval [lo, hi] (a band, or lo == hi == y for
# a line), see `compare_curve`. The simulated curve is interpolated in normalised coordinates
# FIRST; only the resulting error is capped, never the simulated values.
function compare_vertical(ref::ReferenceCurve, xs, ys, tolerance::Tolerance, (los, his))
    length(xs) >= 1 || throw(ArgumentError("need at least one simulated point"))
    finite = [usable(ref.xaxis, x) && usable(ref.yaxis, y) for (x, y) in zip(xs, ys)]
    if !any(finite)
        below_axis(ref, xs, ys) && return all_missed(ref, tolerance)
        throw(ArgumentError("simulation for $(ref.id) has no finite points"))
    end
    issorted(xs[finite]) || throw(ArgumentError("simulated x must be sorted to be compared at matched x"))
    nx = [usable(ref.xaxis, x) ? normalise(ref.xaxis, x) : NaN for x in xs]
    ny = [usable(ref.yaxis, y) ? normalise(ref.yaxis, y) : NaN for y in ys]
    segs = [i for i in 1:(length(xs) - 1) if finite[i] && finite[i + 1]]
    eps = 1e-9
    d = map(zip(ref.x, los, his)) do (x, lo, hi)
        px = normalise(ref.xaxis, x)
        s = NaN
        i = findfirst(i -> finite[i] && abs(nx[i] - px) <= eps, eachindex(xs))     # exact finite match
        if i !== nothing
            s = ny[i]
        else
            k = findfirst(i -> nx[i] < px < nx[i + 1], segs)                       # strictly inside a segment
            if k !== nothing
                j = segs[k]
                s = ny[j] + (px - nx[j]) / (nx[j + 1] - nx[j]) * (ny[j + 1] - ny[j])
            end
        end
        isnan(s) && return 1.0                                                     # outside the range, or in a gap
        nlo, nhi = normalise(ref.yaxis, lo), normalise(ref.yaxis, hi)
        return min(s < nlo ? nlo - s : s > nhi ? s - nhi : 0.0, 1.0)
    end
    r = sqrt(mean(abs2, d))
    q = quantile(d, 0.95)
    return CurveComparison(ref.id, d, r, q, maximum(d), tolerance, r <= tolerance.rms && q <= tolerance.p95)
end

function Base.show(io::IO, c::CurveComparison)
    print(io, c.passed ? "PASS " : "FAIL ", c.id, @sprintf(" rms=%.4f p95=%.4f max=%.4f (tol rms %.3f, p95 %.3f)", c.rms, c.p95, c.max, c.tolerance.rms, c.tolerance.p95))
end
