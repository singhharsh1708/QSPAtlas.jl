# Hutchinson 2016 digitiser wrapper. Same as tools/digitise/digitise.jl but each [[curve]] may also have
#   exclude_rects = [[r1, r2, c1, c2], ...]   rectangles (pixel rows / columns) to ignore, e.g. legends and panel letters
# Run from the QSPAtlas.jl folder:
#   julia --project=tools/digitise models/hutchinson_2016/reference/digitise_hutchinson.jl models/hutchinson_2016/reference/digitise.toml
include(joinpath(@__DIR__, "..", "..", "..", "tools", "digitise", "digitise.jl"))

function trace_rects(spec, root)
    img = loadimg(joinpath(root, spec["image"]))
    r1, r2, c1, c2 = spec["box"]
    tol = get(spec, "tol", 0.25)
    col = spec["colour"]
    match_px = col == "dark" ? (c -> isdark(c; thr = get(spec, "dark_thr", 0.35)) && chroma(c) < 0.2) :
        (target = RGB(col...); c -> colourdist(c, target) <= tol)
    rects = get(spec, "exclude_rects", [])
    excluded(r, c) = any(q -> q[1] <= r <= q[2] && q[3] <= c <= q[4], rects)
    px = Float64[]; py = Float64[]
    for c in c1:c2
        rs = [r for r in r1:r2 if match_px(img[r, c]) && !excluded(r, c)]
        length(rs) < get(spec, "min_pixels", 1) && continue
        push!(px, c); push!(py, median(rs))
    end
    isempty(px) && error("no pixels matched for $(spec["id"])")
    # optional clean-up: keep pixel columns inside the plotted time range, and drop isolated outliers
    # (JPEG blends of other colours) that sit more than `outlier_px` rows from the median of their neighbours
    if haskey(spec, "px_range")
        keep = [spec["px_range"][1] <= c <= spec["px_range"][2] for c in px]
        px = px[keep]; py = py[keep]
    end
    if haskey(spec, "outlier_px")
        w = get(spec, "outlier_window", 9)
        keep = trues(length(px))
        for i in eachindex(px)
            nb = [py[j] for j in eachindex(px) if j != i && abs(px[j] - px[i]) <= w]
            length(nb) >= 2 && abs(py[i] - median(nb)) > spec["outlier_px"] && (keep[i] = false)
        end
        px = px[keep]; py = py[keep]
    end
    fx = calib(spec["xcal"]..., get(spec, "xscale", "linear"))
    fy = calib(spec["ycal"]..., get(spec, "yscale", "linear"))
    n = get(spec, "n", 40)
    idx = unique(round.(Int, range(1, length(px); length = min(n, length(px)))))
    return fx.(px[idx]), fy.(py[idx]), length(px), px[idx], py[idx]
end

function main_rects(tomlpath)
    cfg = TOML.parsefile(tomlpath)
    root = image_root(cfg, tomlpath)
    outdir = dirname(tomlpath)
    for spec in cfg["curve"]
        x, y, npx = trace_rects(spec, root)[1:3]
        write_csv(joinpath(outdir, spec["id"] * ".csv"), spec, x, y, npx)
        println("wrote ", spec["id"], ".csv (", length(x), " points from ", npx, " pixel columns)")
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main_rects(ARGS[1])
