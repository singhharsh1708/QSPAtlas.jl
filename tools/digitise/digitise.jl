# Pixel-level figure digitiser. Run with the tools/digitise environment:
#
#   julia --project=tools/digitise tools/digitise/digitise.jl models/<id>/reference/digitise.toml
#
# Figure images are not in the repository: they are read from $QSPATLAS_PAPERS/<id>/supplement/
# (see tools/digitise/README.md). A digitise.toml selects the folder with `paper = "<id>"`.
#
# Each [[curve]] in digitise.toml names an image, a pixel box for one panel, the
# pixel positions of two known ticks on each axis, and a colour to trace. The tool
# keeps pixels within `tol` of that colour inside the box, takes the median row per
# pixel column, maps pixels to data with the tick calibration, thins to `n` points
# and writes reference/<id>.csv in the format read by QSPAtlas.load_reference.
#
# Helpers for finding calibration pixels interactively (include this file in a REPL):
#   axes_lines(img)            long dark horizontal / vertical lines (panel frames, axes)
#   colour_census(img, box)    most common non-grey colours in a box
#   column_profile(img, col)   dark pixel runs down one column (to locate ticks)
#   row_profile(img, row)      dark pixel runs along one row
using FileIO, ImageIO, Colors, TOML, Statistics

include(joinpath(@__DIR__, "papers.jl"))

loadimg(path) = RGB{Float64}.(load(path))

lum(c) = 0.299red(c) + 0.587green(c) + 0.114blue(c)
isdark(c; thr = 0.45) = lum(c) < thr
chroma(c) = maximum((red(c), green(c), blue(c))) - minimum((red(c), green(c), blue(c)))

"""Rows and columns where at least `frac` of pixels are dark: frame and axis lines."""
function axes_lines(img; frac = 0.35, thr = 0.45)
    H, W = size(img)
    rows = [r for r in 1:H if count(c -> isdark(c; thr), @view img[r, :]) >= frac * W]
    cols = [c for c in 1:W if count(x -> isdark(x; thr), @view img[:, c]) >= frac * H]
    return (; rows, cols)
end

"""Same, but within a box (r1, r2, c1, c2) and relative to the box size."""
function axes_lines(img, box; frac = 0.6, thr = 0.45)
    r1, r2, c1, c2 = box
    sub = @view img[r1:r2, c1:c2]
    a = axes_lines(sub; frac, thr)
    return (; rows = a.rows .+ (r1 - 1), cols = a.cols .+ (c1 - 1))
end

runs(v) = isempty(v) ? Tuple{Int, Int}[] : foldl((acc, x) -> (!isempty(acc) && x == acc[end][2] + 1 ? (acc[end] = (acc[end][1], x); acc) : push!(acc, (x, x))), v; init = Tuple{Int, Int}[])

column_profile(img, col; thr = 0.45) = runs([r for r in axes(img, 1) if isdark(img[r, col]; thr)])
row_profile(img, row; thr = 0.45) = runs([c for c in axes(img, 2) if isdark(img[row, c]; thr)])

"""Most common saturated colours in a box, quantised to 1/8 steps."""
function colour_census(img, box; minchroma = 0.25, top = 8)
    r1, r2, c1, c2 = box
    counts = Dict{NTuple{3, Float64}, Int}()
    for c in img[r1:r2, c1:c2]
        chroma(c) < minchroma && continue
        k = (round(red(c) * 8) / 8, round(green(c) * 8) / 8, round(blue(c) * 8) / 8)
        counts[k] = get(counts, k, 0) + 1
    end
    return first(sort(collect(counts); by = last, rev = true), top)
end

colourdist(a, b) = sqrt((red(a) - red(b))^2 + (green(a) - green(b))^2 + (blue(a) - blue(b))^2)

# pixel -> data with a two-point calibration on linear or log axes
function calib(p1, v1, p2, v2, scale)
    if scale == "log"
        l1, l2 = log10(v1), log10(v2)
        return p -> 10^(l1 + (p - p1) * (l2 - l1) / (p2 - p1))
    else
        return p -> v1 + (p - p1) * (v2 - v1) / (p2 - p1)
    end
end

"""
Trace one curve. `spec` keys:
  image, box = [r1, r2, c1, c2], colour = [r, g, b] (0-1) or "dark", tol,
  xcal = [px1, val1, px2, val2], ycal = [py1, val1, py2, val2], xscale, yscale,
  n (points kept), exclude_rows / exclude_cols (pixel ranges to ignore, e.g. legend),
  mode = "median" (default) | "top" | "bottom" (which pixel to keep per column).
"""
function trace(spec, root)
    img = loadimg(joinpath(root, spec["image"]))
    r1, r2, c1, c2 = spec["box"]
    tol = get(spec, "tol", 0.25)
    col = spec["colour"]
    match_px = col == "dark" ? (c -> isdark(c; thr = get(spec, "dark_thr", 0.35)) && chroma(c) < 0.2) :
        (target = RGB(col...); c -> colourdist(c, target) <= tol)
    exr = [a:b for (a, b) in get(spec, "exclude_rows", [])]
    exc = [a:b for (a, b) in get(spec, "exclude_cols", [])]
    mode = get(spec, "mode", "median")
    px = Float64[]; py = Float64[]
    for c in c1:c2
        any(in(c), exc) && continue
        rs = [r for r in r1:r2 if match_px(img[r, c]) && !any(in(r), exr)]
        isempty(rs) && continue
        push!(px, c)
        push!(py, mode == "top" ? minimum(rs) : mode == "bottom" ? maximum(rs) : median(rs))
    end
    isempty(px) && error("no pixels matched for $(spec["id"])")
    fx = calib(spec["xcal"]..., get(spec, "xscale", "linear"))
    fy = calib(spec["ycal"]..., get(spec, "yscale", "linear"))
    n = get(spec, "n", 40)
    idx = unique(round.(Int, range(1, length(px); length = min(n, length(px)))))
    return fx.(px[idx]), fy.(py[idx]), length(px)
end

function write_csv(path, spec, x, y, npx)
    open(path, "w") do io
        println(io, "# id: ", spec["id"])
        println(io, "# figure: ", spec["figure"])
        println(io, "# x: ", spec["xaxis"])
        println(io, "# y: ", spec["yaxis"])
        println(io, "# kind: ", get(spec, "kind", "curve"))
        println(io, "# source: ", spec["source"])
        println(io, "# method: digitise.jl colour trace, ", npx, " matched pixel columns, box ", spec["box"], ", xcal ", spec["xcal"], ", ycal ", spec["ycal"])
        println(io, "x,y")
        for (a, b) in zip(x, y)
            println(io, round(a; sigdigits = 5), ",", round(b; sigdigits = 5))
        end
    end
end

function main(tomlpath)
    cfg = TOML.parsefile(tomlpath)
    root = image_root(cfg, tomlpath)
    outdir = dirname(tomlpath)
    for spec in cfg["curve"]
        x, y, npx = trace(spec, root)
        write_csv(joinpath(outdir, spec["id"] * ".csv"), spec, x, y, npx)
        println("wrote ", spec["id"], ".csv (", length(x), " points from ", npx, " pixel columns)")
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS[1])
