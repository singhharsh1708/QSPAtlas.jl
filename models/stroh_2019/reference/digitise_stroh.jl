# Colour tracer for stroh_2019 (replaces plain `trace` of tools/digitise/digitise.jl for this model).
#   julia --project=tools/digitise models/stroh_2019/reference/digitise_stroh.jl
# The scans are blurred JPEGs and hollow observation circles share the line colours, so
# a plain colour mask is unreliable. Method per curve (parameters in digitise.toml):
#  1. pixel -> colour class by HUE (blending with white/black keeps hue): a pixel counts for
#     the curve if chroma >= minchroma and its hue is within `hue_tol` degrees of the line colour;
#     pixels inside `exclude_boxes` (legend) are ignored.
#  2. in each column, matching pixels are merged into vertical runs (gap <= 2 px).
#  3. the track starts at the right-most column with exactly one short run and is followed
#     leftwards: next run = the one closest to a linear extrapolation of the last 6 points
#     (must be within `maxdev` px); up to `maxgap` columns without a run are bridged.
#  Columns left of the first matched column / right of the last are simply not part of the target.
using TOML, Statistics
include(joinpath(@__DIR__, "..", "..", "..", "tools", "digitise", "digitise.jl"))

function hue_chroma(p)
    r, g, b = red(p), green(p), blue(p)
    mx, mn = max(r, g, b), min(r, g, b)
    c = mx - mn
    c == 0 && return (0.0, 0.0)
    h = mx == r ? mod((g - b) / c, 6) : mx == g ? (b - r) / c + 2 : (r - g) / c + 4
    return (60 * h, c)
end
hdist(a, b) = (d = abs(mod(a, 360) - mod(b, 360)); min(d, 360 - d))
rgb2hue(c) = hue_chroma(RGB(c...))[1]

function runs_of(rows)
    out = Vector{Vector{Int}}()
    for r in rows
        if !isempty(out) && r - out[end][end] <= 2
            push!(out[end], r)
        else
            push!(out, [r])
        end
    end
    return out
end

function robust_trace(spec, root)
    img = loadimg(joinpath(root, spec["image"]))
    r1, r2, c1, c2 = spec["box"]
    th = get(spec, "hue", rgb2hue(spec["colour"]))
    hue_tol = get(spec, "hue_tol", 14.0); minchroma = get(spec, "minchroma", 0.07)
    maxdev = get(spec, "maxdev", 3.5); maxgap = get(spec, "maxgap", 10); maxrun = get(spec, "maxrun", 40)
    exb = get(spec, "exclude_boxes", [])
    function matched(r, c)
        any(b -> b[1] <= r <= b[2] && b[3] <= c <= b[4], exb) && return false
        h, ch = hue_chroma(img[r, c])
        return ch >= minchroma && hdist(h, th) <= hue_tol
    end
    rowsat = Dict{Int, Vector{Int}}()
    for c in c1:c2
        rs = [r for r in r1:r2 if matched(r, c)]
        isempty(rs) || (rowsat[c] = rs)
    end
    function follow(seedc)
        rr0 = runs_of(rowsat[seedc])
        cols = [seedc]; ys = [mean(rr0[1])]
        gap = 0
        for c in (seedc - 1):-1:c1
            k = min(8, length(cols))
            sl = 0.0
            pred = ys[end]
            if k >= 3
                xs = Float64.(cols[end-k+1:end]); yy = ys[end-k+1:end]
                sl = sum((xs .- mean(xs)) .* (yy .- mean(yy))) / sum((xs .- mean(xs)) .^ 2)
                pred = ys[end] + sl * (c - cols[end])
            end
            w = maxdev + abs(sl) + 0.15 * (cols[end] - c)
            if haskey(rowsat, c)
                near = [r for r in rowsat[c] if abs(r - pred) <= w]
                if !isempty(near)
                    push!(cols, c); push!(ys, mean(near)); gap = 0
                    continue
                end
            end
            gap += 1
            gap > maxgap && break
        end
        return cols, ys
    end
    best = nothing
    for c in c2:-1:c1
        haskey(rowsat, c) || continue
        rr = runs_of(rowsat[c])
        (length(rr) == 1 && length(rr[1]) <= maxrun) || continue
        cols, ys = follow(c)
        if length(cols) >= get(spec, "minlen", 10)
            best = (cols, ys); break
        end
    end
    best === nothing && error("no seed with a long enough track for $(spec["id"])")
    cols, ys = best
    o = sortperm(cols)
    cols = Float64.(cols[o]); ys = ys[o]
    # outlier filter: drop columns that deviate > `outlier` px from the rolling median (+-5 columns),
    # which removes places where a hollow circle or a crossing curve dragged the run centre
    keep = [abs(ys[i] - median(ys[max(1, i - 5):min(end, i + 5)])) <= get(spec, "outlier", 2.0) + 0.5 * abs(ys[min(end, i + 5)] - ys[max(1, i - 5)]) / 10 * 5 for i in eachindex(ys)]
    return cols[keep], ys[keep]
end

function main(tomlpath)
    cfg = TOML.parsefile(tomlpath)
    root = image_root(cfg, tomlpath)
    outdir = dirname(tomlpath)
    for spec in cfg["curve"]
        px, py = robust_trace(spec, root)
        fx = calib(spec["xcal"]..., spec["xscale"]); fy = calib(spec["ycal"]..., spec["yscale"])
        n = spec["n"]
        keep = fx.(px) .>= get(spec, "min_x", -Inf)     # marker-contamination rule, see make_toml.jl
        px, py = px[keep], py[keep]
        idx = unique(round.(Int, range(1, length(px); length = min(n, length(px)))))
        write_csv(joinpath(outdir, spec["id"] * ".csv"), spec, fx.(px[idx]), fy.(py[idx]), length(px))
        println(rpad(spec["id"], 14), length(px), " pixel columns, x ", round(fx(px[1]); digits = 2), "..", round(fx(px[end]); digits = 2))
    end
end
abspath(PROGRAM_FILE) == (@__FILE__) && main(length(ARGS) > 0 ? ARGS[1] : joinpath(@__DIR__, "digitise.toml"))
