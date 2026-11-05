"""
    Tolerance(; rms = 0.015, p95 = 0.03)

Pass rule for one curve. Errors are distances in the figure's own axis box, where
both axes are scaled to `[0, 1]` (log axes in decades). `rms` bounds the root mean
square distance over all reference points and `p95` bounds the 95th percentile, so
one badly traced point cannot fail a curve on its own. The maximum is still reported.

The defaults rest on a measured digitising error study (tools/digitise/error_study/):
twice the worst panel type's 99th percentile of digitising-only error gives RMS 0.015
and P95 0.031 (the atlas keeps 0.03). Panels are 85 to 294 px, so the tolerance is
roughly 1.3 to 3.5 px vertically. See docs/src/method.md and the sensitivity tables. Tolerances are frozen in each
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

# ---- report-only diagnostics ---------------------------------------------------------------
# Checks reported next to the pass/fail of `compare_curve`. They never enter the pass rule or the
# grade; they are stored in score.toml so that a later, separate decision can use them.

"""
    CurveDiagnostics

Report-only checks of one scored curve (none of them changes `passed` or any grade). All
distances are in axis-box units, like [`CurveComparison`](@ref). `NaN` (and `0` for the counts)
means "does not apply" or "nothing to measure".

- `coverage`: span of the digitised reference x values divided by the axis x range (both in
  the displayed coordinate, so decades on a log axis), clamped to `[0, 1]`. A reference that
  covers 40 % of the axis only constrains the simulation there. `n_points` is the number of
  reference points.
- `sim_coverage`: fraction of the axis x range on which the simulated curve is defined, i.e. the
  union of the segments between neighbouring simulated points that are both finite and placeable
  on the axes (a gap of non-finite points is not bridged), clipped to the axis.
- `rev_rms`, `rev_p95`, `rev_n` (`kind = :curve` only): the reverse distance. For every simulated
  point that is placeable and whose x lies inside the reference's x range, the shortest distance
  to the reference polyline (reference points in file order), capped at 1 like the forward error.
  This is the other side of the two-sided check: a simulation that wanders away from the
  reference where the reference has no points shows up here but not in the forward distance. It
  is measured at simulated points, so a sparse reference polyline can inflate it on strongly curved
  sections. `rev_n` is the number of simulated points used.
- `band_centre_dev_rms`, `band_n` (`kind = :band` only): root mean square over the reference x,
  where the simulation is defined, of `|s - centre| / half_width`: how far the simulated value `s`
  (interpolated as for the band score) is from the band centre, in units of the band half-width,
  in the displayed coordinate. 1 is the band edge; the band score is 0 anywhere up to 1, so this
  separates "inside but off-centre" from "on the centre line". Reference x with a zero half-width
  are left out. `band_n` counts the reference x where the simulation is defined.
- `band_edge_rms`, `band_edge_max` (only when the model passed simulated interval edges, see
  [`curve_diagnostics`](@ref)): distance of the simulated lower edge to the reference lower edge and
  of the simulated upper edge to the reference upper edge at every reference x where both are
  defined, capped at 1; root mean square and maximum over both edges.
"""
struct CurveDiagnostics
    coverage::Float64
    n_points::Int
    sim_coverage::Float64
    rev_rms::Float64
    rev_p95::Float64
    rev_n::Int
    band_centre_dev_rms::Float64
    band_n::Int
    band_edge_rms::Float64
    band_edge_max::Float64
end

# span of the reference x values as a fraction of the axis, in displayed coordinates
function ref_coverage(ref::ReferenceCurve)
    nx = [normalise(ref.xaxis, x) for x in ref.x if isfinite(x)]
    isempty(nx) && return 0.0
    return clamp(maximum(nx) - minimum(nx), 0.0, 1.0)
end

# fraction of the axis x range covered by segments joining neighbouring placeable simulated points
function sim_x_coverage(ref::ReferenceCurve, xs, ys)
    ivs = Tuple{Float64, Float64}[]
    ok = [usable(ref.xaxis, x) && usable(ref.yaxis, y) for (x, y) in zip(xs, ys)]
    for i in 1:(length(xs) - 1)
        (ok[i] && ok[i + 1]) || continue
        a, b = minmax(normalise(ref.xaxis, xs[i]), normalise(ref.xaxis, xs[i + 1]))
        a, b = max(a, 0.0), min(b, 1.0)
        b > a && push!(ivs, (a, b))
    end
    isempty(ivs) && return 0.0
    sort!(ivs)
    total, (cur_lo, cur_hi) = 0.0, ivs[1]
    for (a, b) in ivs[2:end]
        if a <= cur_hi
            cur_hi = max(cur_hi, b)
        else
            total += cur_hi - cur_lo
            cur_lo, cur_hi = a, b
        end
    end
    return clamp(total + cur_hi - cur_lo, 0.0, 1.0)
end

# distance of every simulated point inside the reference x range to the reference polyline
function reverse_distances(ref::ReferenceCurve, xs, ys)
    xlo, xhi = extrema(ref.x)
    rp = [(normalise(ref.xaxis, x), normalise(ref.yaxis, y)) for (x, y) in zip(ref.x, ref.y)]
    d = Float64[]
    for (x, y) in zip(xs, ys)
        (usable(ref.xaxis, x) && usable(ref.yaxis, y) && xlo <= x <= xhi) || continue
        p = (normalise(ref.xaxis, x), normalise(ref.yaxis, y))
        dd = length(rp) == 1 ? norm(p .- rp[1]) : minimum(segdist(p, rp[i], rp[i + 1]) for i in 1:(length(rp) - 1))
        push!(d, min(dd, 1.0))
    end
    return d
end

# Simulated value (normalised y) at every reference x, NaN where it is not defined. Same rules as
# `compare_vertical`: an exact finite match first, else linear interpolation inside a segment of two
# neighbouring finite points, never across a gap or outside the simulated range.
function sim_at_reference_x(ref::ReferenceCurve, xs, ys)
    finite = [usable(ref.xaxis, x) && usable(ref.yaxis, y) for (x, y) in zip(xs, ys)]
    nx = [usable(ref.xaxis, x) ? normalise(ref.xaxis, x) : NaN for x in xs]
    ny = [usable(ref.yaxis, y) ? normalise(ref.yaxis, y) : NaN for y in ys]
    segs = [i for i in 1:(length(xs) - 1) if finite[i] && finite[i + 1]]
    return map(ref.x) do x
        px = normalise(ref.xaxis, x)
        i = findfirst(i -> finite[i] && abs(nx[i] - px) <= 1e-9, eachindex(xs))
        i !== nothing && return ny[i]
        k = findfirst(i -> nx[i] < px < nx[i + 1], segs)
        k === nothing && return NaN
        j = segs[k]
        return ny[j] + (px - nx[j]) / (nx[j + 1] - nx[j]) * (ny[j + 1] - ny[j])
    end
end

rms_(v) = isempty(v) ? NaN : sqrt(sum(abs2, v) / length(v))

"""
    curve_diagnostics(ref, xs, ys; sim_edges = nothing) -> CurveDiagnostics

Report-only checks of a simulated curve `(xs, ys)` against `ref` (see [`CurveDiagnostics`](@ref)):
coverage of the reference and of the simulation, the reverse distance for `:curve` references, and
the distance from the band centre for `:band` references. Nothing here is used for pass/fail.

`sim_edges = (lo = ..., hi = ...)` is optional and only read for `:band` references: the model's own
simulated interval edges (for example the `lo` and `hi` of a VPC) at the simulated `xs`, vectors of the
same length as `xs`. They are compared with the reference edges at the reference x (same
interpolation and gap rules as the band score) and reported as `band_edge_rms` and `band_edge_max`.
[`run_curves`](@ref) passes them on when the model's function returns `(x, y, (lo = lo, hi = hi))`.
"""
function curve_diagnostics(ref::ReferenceCurve, xs::AbstractVector, ys::AbstractVector; sim_edges = nothing)
    length(xs) == length(ys) || throw(DimensionMismatch("xs and ys differ in length"))
    cov = ref_coverage(ref)
    scov = sim_x_coverage(ref, xs, ys)
    rev_rms, rev_p95, rev_n = NaN, NaN, 0
    if ref.kind === :curve
        d = reverse_distances(ref, xs, ys)
        if !isempty(d)
            rev_rms, rev_p95, rev_n = rms_(d), quantile(d, 0.95), length(d)
        end
    end
    centre_dev, band_n, edge_rms, edge_max = NaN, 0, NaN, NaN
    if ref.kind === :band
        s = sim_at_reference_x(ref, xs, ys)
        nlo = [normalise(ref.yaxis, v) for v in ref.ylo]
        nhi = [normalise(ref.yaxis, v) for v in ref.yhi]
        band_n = count(!isnan, s)
        dev = [abs(si - (l + h) / 2) / ((h - l) / 2) for (si, l, h) in zip(s, nlo, nhi) if !isnan(si) && h > l]
        centre_dev = rms_(dev)
        if sim_edges !== nothing
            (haskey(sim_edges, :lo) && haskey(sim_edges, :hi)) || throw(ArgumentError("sim_edges must be (lo = ..., hi = ...)"))
            (length(sim_edges.lo) == length(xs) && length(sim_edges.hi) == length(xs)) ||
                throw(DimensionMismatch("sim_edges.lo and sim_edges.hi must have the length of xs"))
            slo = sim_at_reference_x(ref, xs, sim_edges.lo)
            shi = sim_at_reference_x(ref, xs, sim_edges.hi)
            e = Float64[]
            for i in eachindex(slo)
                (isnan(slo[i]) || isnan(shi[i])) && continue
                push!(e, min(abs(slo[i] - nlo[i]), 1.0), min(abs(shi[i] - nhi[i]), 1.0))
            end
            isempty(e) || ((edge_rms, edge_max) = (rms_(e), maximum(e)))
        end
    end
    return CurveDiagnostics(cov, length(ref.x), scov, rev_rms, rev_p95, rev_n, centre_dev, band_n, edge_rms, edge_max)
end
