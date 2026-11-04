"""
    Axis(label, unit, scale, lo, hi)

One plot axis of a published figure. `scale` is `:linear` or `:log`; `lo` and `hi`
are the axis limits as drawn in the figure. They set the scale used for comparison,
so an error of 0.01 means one percent of the drawn axis.
"""
struct Axis
    label::String
    unit::String
    scale::Symbol
    lo::Float64
    hi::Float64
    function Axis(label, unit, scale, lo, hi)
        scale in (:linear, :log) || throw(ArgumentError("axis scale must be :linear or :log, got $scale"))
        hi > lo || throw(ArgumentError("axis upper limit $hi must exceed lower limit $lo"))
        scale === :log && lo <= 0 && throw(ArgumentError("log axis needs a positive lower limit, got $lo"))
        return new(String(label), String(unit), scale, Float64(lo), Float64(hi))
    end
end

"""
    normalise(axis, v)

Map a data value to the fraction of the drawn axis (0 at `lo`, 1 at `hi`). Log axes
are mapped in decades. Values below zero on a log axis map to `-Inf` clamped at -1,
which counts as a large error rather than an exception.
"""
function normalise(a::Axis, v::Real)
    if a.scale === :linear
        return (v - a.lo) / (a.hi - a.lo)
    else
        v <= 0 && return -1.0
        return (log10(v) - log10(a.lo)) / (log10(a.hi) - log10(a.lo))
    end
end

"""
    ReferenceCurve

Points digitised from one curve of one published figure panel.

Fields: `id` (unique within the model, e.g. `fig4_left_vehicle_P`), `figure`
(human description, e.g. "Figure 4 left, vehicle, tumour volume P"), `x`, `y`,
`xaxis`, `yaxis` ([`Axis`](@ref)), `kind` (`:curve` for a line drawn by the authors'
simulation, `:points` for plotted observations, `:bars` for bar heights or reported
single values at category positions, `:band` for a shaded interval such as the authors'
simulated 5th to 95th percentile band), `meta` (all other header keys,
for example `source` and `method`), and `ylo`, `yhi` (lower and upper edge at each `x`;
only for `:band`, empty otherwise). For a band, `y` is the midpoint `(ylo + yhi) / 2`,
which is only used for plotting and the point count.

The eight-argument constructor (without `ylo`, `yhi`) is kept for every other kind.
"""
struct ReferenceCurve
    id::String
    figure::String
    x::Vector{Float64}
    y::Vector{Float64}
    xaxis::Axis
    yaxis::Axis
    kind::Symbol
    meta::Dict{String, String}
    ylo::Vector{Float64}
    yhi::Vector{Float64}
    function ReferenceCurve(id, figure, x, y, xaxis, yaxis, kind, meta, ylo, yhi)
        if kind === :band
            length(ylo) == length(x) == length(yhi) == length(y) ||
                throw(ArgumentError("band reference $id: x, y, ylo and yhi must have the same length"))
            all(ylo .<= yhi) || throw(ArgumentError("band reference $id: every lower edge must be <= its upper edge"))
        else
            (isempty(ylo) && isempty(yhi)) || throw(ArgumentError("reference $id: ylo and yhi are only for kind = :band"))
        end
        return new(id, figure, x, y, xaxis, yaxis, kind, meta, ylo, yhi)
    end
end

ReferenceCurve(id, figure, x, y, xaxis, yaxis, kind, meta) =
    ReferenceCurve(id, figure, x, y, xaxis, yaxis, kind, meta, Float64[], Float64[])

function parse_axis(s::AbstractString)
    # "time [day] linear 0 60"
    m = match(r"^(.*?)\s*\[(.*?)\]\s+(linear|log)\s+(\S+)\s+(\S+)\s*$", strip(s))
    m === nothing && throw(ArgumentError("cannot parse axis line \"$s\"; expected: label [unit] linear|log lo hi"))
    return Axis(m[1], m[2], Symbol(m[3]), parse(Float64, m[4]), parse(Float64, m[5]))
end

format_axis(a::Axis) = string(a.label, " [", a.unit, "] ", a.scale, " ", a.lo, " ", a.hi)

"""
    load_reference(path) -> ReferenceCurve

Read a reference CSV. The file starts with `# key: value` header lines and then has a
two-column `x,y` table. Required keys: `id`, `figure`, `x`, `y`. Optional: `kind`
(`curve`, `points`, `bars` or `band`, default `curve`) and anything else, kept in `meta`.
A `# kind: band` file has three columns `x,lo,hi` instead (a shaded interval, lower and
upper edge at each x); the reference's `y` is then the midpoint of `lo` and `hi`. A `curve`
may carry `# metric: vertical` (default `nearest`): compared by vertical distance at the
reference x, as for bands (see [`compare_curve`](@ref)); use it for simulated percentile lines.

```
# id: fig4_left_vehicle_P
# figure: Figure 4 left, vehicle, tumour volume P
# x: time [day] linear 0 60
# y: tumour volume [mm^3] linear 0 2500
# kind: curve
# source: PSP4-4-660-g004.jpg (PMC4716582)
# method: colour mask, see digitise.toml
x,y
0.0,150.2
...
```
"""
function load_reference(path::AbstractString)
    header = Dict{String, String}()
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        s = strip(line)
        isempty(s) && continue
        if startswith(s, "#")
            kv = split(lstrip(s, ['#', ' ']), ':'; limit = 2)
            length(kv) == 2 && (header[strip(kv[1])] = strip(kv[2]))
        elseif lowercase(s) in ("x,y", "x,lo,hi")
            continue
        else
            push!(rows, parse.(Float64, strip.(split(s, ','))))
        end
    end
    for k in ("id", "figure", "x", "y")
        haskey(header, k) || throw(ArgumentError("$path: missing header key \"$k\""))
    end
    isempty(rows) && throw(ArgumentError("$path: no data rows"))
    kind = Symbol(get(header, "kind", "curve"))
    kind in (:curve, :points, :bars, :band) || throw(ArgumentError("$path: kind must be curve, points, bars or band"))
    ncol = kind === :band ? 3 : 2
    all(r -> length(r) == ncol, rows) ||
        throw(ArgumentError("$path: kind $kind needs $ncol columns (" * (kind === :band ? "x,lo,hi" : "x,y") * ")"))
    meta = Dict(k => v for (k, v) in header if !(k in ("id", "figure", "x", "y", "kind")))
    # `# metric: vertical` (curve references only) scores by vertical distance at the reference x
    # instead of the default nearest-polyline distance, see `compare_curve`
    if haskey(meta, "metric")
        meta["metric"] in ("vertical", "nearest") || throw(ArgumentError("$path: metric must be vertical or nearest"))
        kind === :curve || throw(ArgumentError("$path: a metric is only for kind curve"))
    end
    meta["file"] = basename(path)
    x = [r[1] for r in rows]
    if kind === :band
        lo = [r[2] for r in rows]
        hi = [r[3] for r in rows]
        return ReferenceCurve(header["id"], header["figure"], x, (lo .+ hi) ./ 2,
            parse_axis(header["x"]), parse_axis(header["y"]), kind, meta, lo, hi)
    end
    return ReferenceCurve(header["id"], header["figure"], x, [r[2] for r in rows],
        parse_axis(header["x"]), parse_axis(header["y"]), kind, meta)
end

"""
    load_references(dir) -> Dict{String, ReferenceCurve}

Load every `*.csv` in `dir`, keyed by curve id.
"""
function load_references(dir::AbstractString)
    out = Dict{String, ReferenceCurve}()
    for f in sort(readdir(dir))
        endswith(f, ".csv") || continue
        r = load_reference(joinpath(dir, f))
        haskey(out, r.id) && throw(ArgumentError("duplicate reference id $(r.id) in $dir"))
        out[r.id] = r
    end
    return out
end

"""
    write_reference(path, ref::ReferenceCurve)

Write a reference curve in the format read by [`load_reference`](@ref).
"""
function write_reference(path::AbstractString, r::ReferenceCurve)
    open(path, "w") do io
        println(io, "# id: ", r.id)
        println(io, "# figure: ", r.figure)
        println(io, "# x: ", format_axis(r.xaxis))
        println(io, "# y: ", format_axis(r.yaxis))
        println(io, "# kind: ", r.kind)
        for k in sort(collect(keys(r.meta)))
            k == "file" && continue
            println(io, "# ", k, ": ", r.meta[k])
        end
        if r.kind === :band
            println(io, "x,lo,hi")
            for (x, lo, hi) in zip(r.x, r.ylo, r.yhi)
                println(io, @sprintf("%.6g,%.6g,%.6g", x, lo, hi))
            end
        else
            println(io, "x,y")
            for (x, y) in zip(r.x, r.y)
                println(io, @sprintf("%.6g,%.6g", x, y))
            end
        end
    end
    return path
end
