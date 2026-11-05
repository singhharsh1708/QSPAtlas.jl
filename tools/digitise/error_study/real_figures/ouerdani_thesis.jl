# Re-digitise Ouerdani 2015 Figure 4 (both panels, P solid and K dashed) from the 300 ppi raster of the
# same figure in the first author's thesis (HAL tel-01496174, PDF page 96, printed page 95). The thesis
# embeds the figure as two strips (top 437 rows, bottom 436 rows) that are stacked here.
#
# Prepare (not part of the repository):
#   pdfimages -png -f 96 -l 96 tel-01496174.pdf <dir>/p96      -> p96-000.png, p96-001.png
# Run:
#   julia --project=tools/digitise tools/digitise/error_study/real_figures/ouerdani_thesis.jl <dir> <output folder>
#
# Calibration (pixel centres, 1-based rows/cols of the stacked image; ticks are the inward tick marks,
# centroid of the dark run in a column/row just inside the frame; least-squares line through all ticks):
#   left  x: ticks 10..70 d at cols 280.5, 382.5, 484.5, 586.5, 689.0, 790.5, 892.5 (row 785)
#   left  y: ticks 50..200 at rows 635.5, 481.5, 327.0, 173.5 (cols 938/940, right frame)
#   right x: ticks 10..80 mo at cols 1122.5, 1215.5, 1309.5, 1402.5, 1495.5, 1588.75, 1682.25, 1775.5
#   right y: ticks 45..65 at rows 662.5, 533.5, 404.5, 275.5, 146.5 (cols 1788/1790)
# The tick rows/columns are re-measured by the script (centroids) and printed with the fit residuals.
#
# Curve rule (the same idea as models/ouerdani_2015/reference/trace_fig4.jl, at 3.3 times the resolution):
# per pixel column from 4 px inside the frame, dark runs (luminance < 0.5; runs that are tick marks are dropped), darkness-weighted centroid per run.
# P (solid, continuous) is followed from the left edge: the run nearest to a linear prediction from the
# last P points. Every other run is K (dashed). Columns where another run lies within MERGE px of the P run,
# or where the P run is much longer than its neighbours (the two lines merged at the crossing), are skipped.
# Isolated K runs (fewer than 3 other K runs within 8 columns and 10 rows) are dropped.
using Statistics, Printf
include(joinpath(@__DIR__, "..", "..", "digitise.jl"))

const SRC = ARGS[1]
const OUT = ARGS[2]
mkpath(OUT)
img = vcat(loadimg(joinpath(SRC, "p96-000.png")), loadimg(joinpath(SRC, "p96-001.png")))

wcentroid(rr, f) = sum(r * (1 - f(r)) for r in rr) / sum(1 - f(r) for r in rr)

# tick centroid along a column (rows) or a row (columns), dark runs of <= 8 px only
function tickpos(kind, line, around; half = 8)
    if kind == :col
        rr = [r for r in (around - half):(around + half) if lum(img[r, line]) < 0.5]
        return wcentroid(rr, r -> lum(img[r, line]))
    else
        cc = [c for c in (around - half):(around + half) if lum(img[line, c]) < 0.5]
        return wcentroid(cc, c -> lum(img[line, c]))
    end
end
function fitline(px, v)
    b = sum((px .- mean(px)) .* (v .- mean(v))) / sum((px .- mean(px)) .^ 2)
    a = mean(v) - b * mean(px)
    return (p -> a + b * p), maximum(abs.(a .+ b .* px .- v))
end

panels = Dict(
    :left => (cols = 183:940, ticks = [(172, 176), (326, 329), (480, 484), (634, 638)], rows = 29:781, framex = (179, 944),
        xt = ([281, 383, 485, 587, 689, 791, 893], 10:10:70, 785), yt = ([636, 481, 327, 174], 50:50:200, (938, 940))),
    :right => (cols = 1034:1790, ticks = [(145, 149), (274, 278), (403, 407), (532, 536), (661, 665)], rows = 29:781, framex = (1030, 1794),
        xt = ([1123, 1216, 1310, 1403, 1496, 1589, 1682, 1776], 10:10:80, 785), yt = ([663, 534, 405, 276, 147], 45:5:65, (1788, 1790))),
)

function calibrate(P, name)
    xs = [tickpos(:row, P.xt[3], c) for c in P.xt[1]]
    ys = [mean(tickpos(:col, c, r) for c in P.yt[3]) for r in P.yt[1]]
    fx, rx = fitline(xs, collect(Float64, P.xt[2]))
    fy, ry = fitline(ys, collect(Float64, P.yt[2]))
    @printf("%s: x ticks at %s, residual %.3f; y ticks at %s, residual %.3f (data units)\n", name, string(round.(xs; digits = 2)), rx, string(round.(ys; digits = 2)), ry)
    return fx, fy
end

const MERGE = 3.0
function trace_panel(P)
    pts_P = Tuple{Float64, Float64}[]; pts_K = Tuple{Float64, Float64}[]
    lenP = Float64[]
    for c in P.cols
        rs = [r for r in P.rows if lum(img[r, c]) < 0.5]
        isempty(rs) && continue
        rr = [(a, b, wcentroid(a:b, r -> lum(img[r, c]))) for (a, b) in runs(rs)]
        # within 12 px of the left or right frame, drop runs that are tick marks (inside a tick row band)
        if c < P.framex[1] + 12 || c > P.framex[2] - 12
            rr = [x for x in rr if !any(t -> t[1] <= x[1] && x[2] <= t[2], P.ticks)]
        end
        isempty(rr) && continue
        if isempty(pts_P)
            length(rr) == 1 || continue
            push!(pts_P, (c, rr[1][3])); push!(lenP, rr[1][2] - rr[1][1] + 1); continue
        end
        k = min(5, length(pts_P))
        last = pts_P[end-k+1:end]
        pred = last[end][2]
        xs = [p[1] for p in last]; ys = [p[2] for p in last]
        sl = k >= 2 ? sum((xs .- mean(xs)) .* (ys .- mean(ys))) / max(sum((xs .- mean(xs)) .^ 2), 1e-9) : 0.0
        pred = ys[end] + sl * (c - xs[end])
        i = argmin([abs(x[3] - pred) for x in rr])
        L = rr[i][2] - rr[i][1] + 1
        others = [x for (j, x) in enumerate(rr) if j != i]
        merged = any(abs(x[3] - rr[i][3]) < MERGE + (L + x[2] - x[1] + 1) / 2 for x in others) ||
                 (L > 1.6 * median(lenP[max(1, end - 9):end]) + 4 && !isempty(pts_K))
        if merged
            continue
        end
        push!(pts_P, (c, rr[i][3])); push!(lenP, L)
        for x in others
            push!(pts_K, (c, x[3]))
        end
    end
    # K clean-up: a dash pixel run has neighbours on the dashed line; drop isolated runs (pieces of the
    # steep P line split by anti-aliasing at the left edge), i.e. fewer than 3 other K runs within 8 columns
    # and 10 rows
    keep = [count(q -> q !== p && abs(q[1] - p[1]) <= 8 && abs(q[2] - p[2]) <= 10, pts_K) >= 3 for p in pts_K]
    println("  K runs dropped as isolated: ", [(Int(p[1]), round(p[2]; digits = 1)) for (p, k) in zip(pts_K, keep) if !k])
    return pts_P, pts_K[keep]
end

ids = Dict(:left => ("fig4_left_P", "fig4_left_K"), :right => ("fig4_right_P", "fig4_right_K"))
for name in (:left, :right)
    P = panels[name]
    fx, fy = calibrate(P, name)
    pP, pK = trace_panel(P)
    for (id, pts) in zip(ids[name], (pP, pK))
        pts = sort(pts)
        open(joinpath(OUT, id * ".csv"), "w") do io
            println(io, "x,y")
            for (c, r) in pts
                @printf(io, "%.6g,%.6g\n", fx(c), fy(r))
            end
        end
        println(id, ": ", length(pts), " points, x ", round(fx(pts[1][1]); digits = 2), "..", round(fx(pts[end][1]); digits = 2))
    end
end
