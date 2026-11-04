# Digitise the shaded simulated-percentile bands of Figure 3a (PSP4-8-285-g003.jpg, PMC6539725).
# Run:  julia --startup-file=no --project=tools/digitise models/stein_2019/reference/trace_fig3a.jl
#
# Bands are light-blue (10th and 90th percentile) and light-pink (50th percentile) areas.
# For each pixel column inside the panel, runs of band-coloured pixels are found (gaps of up to
# GAP rows are bridged: data dots, the empirical percentile line and red outlier marks cross the
# bands). The upper blue run is the 90th-percentile band, the lower blue run the 10th, the pink
# run the 50th. Edges are the outer pixel boundaries of the run. A 5-column running median
# removes single-column spikes. y calibration uses the printed tick labels (centres of the
# labelled rows), x calibration the printed day labels / vertical grid lines, see digitise.toml.
# Output is thinned to one point per STEP pixel columns.
include(joinpath(@__DIR__, "..", "..", "..", "tools", "digitise", "digitise.jl"))
using Statistics, FileIO

# The image is not in the repository: it is read from $QSPATLAS_PAPERS/stein_2019/supplement/ (tools/digitise/README.md).
img = loadimg(paper_image("stein_2019", "PSP4-8-285-g003.jpg"))

isblue(c) = blue(c) > 0.93 && 0.6 < red(c) < 0.86 && 0.74 < green(c) < 0.92
ispink(c) = red(c) > 0.93 && 0.68 < green(c) < 0.88 && 0.68 < blue(c) < 0.9 && abs(green(c) - blue(c)) < 0.08

function runs_col(pred, c, r1, r2; gap = 10, minlen = 3)
    rs = [r for r in r1:r2 if pred(img[r, c])]
    isempty(rs) && return Tuple{Int, Int}[]
    out = Tuple{Int, Int}[]
    s = p = rs[1]
    for r in rs[2:end]
        if r - p - 1 <= gap
            p = r
        else
            push!(out, (s, p)); s = r; p = r
        end
    end
    push!(out, (s, p))
    return filter(x -> x[2] - x[1] + 1 >= minlen, out)
end

function runmedian(v, w = 5)
    return map(eachindex(v)) do i
        win = filter(!isnan, v[max(1, i - w ÷ 2):min(end, i + w ÷ 2)])
        (isnan(v[i]) || isempty(win)) ? NaN : median(win)
    end
end

# panel geometry: pixel columns of the plotting area, rows, calibrations (see digitise.toml)
panels = Dict(
    "left" => (cols = 53:290, rows = (18, 220),
        xcal = (82.5, 50.0, 275.5, 350.0), ycal = (90.5, 1e4, 192.5, 1e2)),
    "right" => (cols = 393:637, rows = (18, 214),
        xcal = (406.0, 2.0, 634.0, 30.0), ycal = (89.0, 1e4, 185.0, 1e2)),
)
ydata(cal, row) = 10^(log10(cal[2]) + (row - cal[1]) * (log10(cal[4]) - log10(cal[2])) / (cal[3] - cal[1]))
xdata(cal, col) = cal[2] + (col - cal[1]) * (cal[4] - cal[2]) / (cal[3] - cal[1])

# per column: (top, bottom) pixel rows of each band (NaN if not found)
function trace_panel(name)
    P = panels[name]
    r1, r2 = P.rows
    cols = collect(P.cols)
    out = Dict(k => fill(NaN, length(cols), 2) for k in ("p90", "p50", "p10"))
    for (j, c) in enumerate(cols)
        B = runs_col(isblue, c, r1, r2; minlen = 5)
        K = runs_col(ispink, c, r1, r2; minlen = 6)
        if !isempty(B)
            # 90th band: topmost run of at least 12 px; 10th: lowest run (at least 5 px) below it
            i90 = findfirst(x -> x[2] - x[1] + 1 >= 12, B)
            i90 !== nothing && (out["p90"][j, :] .= (B[i90][1] - 0.5, B[i90][2] + 0.5))
            lo = [b for (k, b) in enumerate(B) if i90 === nothing || k > i90]
            if !isempty(lo)
                out["p10"][j, :] .= (lo[end][1] - 0.5, lo[end][2] + 0.5)
            end
        end
        if !isempty(K)
            k = argmax([x[2] - x[1] for x in K])
            out["p50"][j, :] .= (K[k][1] - 0.5, K[k][2] + 0.5)
        end
    end
    for k in keys(out), s in 1:2
        out[k][:, s] .= runmedian(out[k][:, s])
    end
    return cols, out
end

# ---- targets: which part of each band is scored, and why --------------------------------
# left (0-380 d): x in [30, 185] only. Left of 30 d the band lies under the dense data cloud and
#   the right panel covers it; right of ~187 d the band is a flat extension of the last bin
#   (the empirical percentile lines stop there), not a model percentile at that time.
# right (0-30 d): p90 from 4.5 d (first bin extension is flat before the first node at ~2.8 d,
#   and 3-4.5 d is crossed by the p50 band); p50 from 3.4 d except 4.4-6.4 d where red outlier
#   marks (dark red areas) cover the band; p10 from 6.4 d (lower edge leaves the axis box before
#   that and red outlier marks cover the band at 5-6 d). All to 29.5 d (frame edge).
const TARGETS = [
    # id, panel, band, xlo, xhi, excluded intervals, step in pixel columns
    ("fig3a_right_p90", "right", "p90", 4.5, 29.5, [], 3),
    ("fig3a_right_p50", "right", "p50", 3.4, 29.5, [(4.4, 6.4)], 3),
    ("fig3a_right_p10", "right", "p10", 6.4, 29.5, [], 3),
    ("fig3a_left_p90", "left", "p90", 30.0, 185.0, [], 6),
    ("fig3a_left_p50", "left", "p50", 30.0, 185.0, [], 6),
    ("fig3a_left_p10", "left", "p10", 30.0, 185.0, [], 6),
]
const PANEL_AXES = Dict(
    "left" => (x = "time since infusion [day] linear -18 383", y = "tisagenlecleucel transgene [copies/ug] log 27.6 300000"),
    "right" => (x = "time since infusion [day] linear 0.28 31", y = "tisagenlecleucel transgene [copies/ug] log 22.6 420000"),
)
const BANDNAME = Dict("p90" => "90th percentile (blue band)", "p50" => "50th percentile (pink band)", "p10" => "10th percentile (blue band)")

function write_band(path, id, panel, band, x, lo, hi, npx)
    open(path, "w") do io
        println(io, "# id: ", id)
        println(io, "# figure: Figure 3a ", panel, " panel, shaded simulated interval of the ", BANDNAME[band])
        println(io, "# x: ", PANEL_AXES[panel].x)
        println(io, "# y: ", PANEL_AXES[panel].y)
        println(io, "# kind: band")
        println(io, "# source: PSP4-8-285-g003.jpg (PMC6539725), authors' simulation (Monolix VPC prediction interval), NOT the empirical line or dots")
        println(io, "# method: trace_fig3a.jl colour-run trace of band-coloured pixels, ", npx, " columns kept, calibration in digitise.toml; band edges = outer pixel boundaries")
        println(io, "x,lo,hi")
        for (a, b, c) in zip(x, lo, hi)
            println(io, round(a; sigdigits = 5), ",", round(b; sigdigits = 4), ",", round(c; sigdigits = 4))
        end
    end
end

# ADDED AFTER REVIEW (2026-11-04, after the first comparison): a band-edge point is dropped when its whole
# interval lies away from its neighbours' interval, i.e. hi < 1/FACTOR times the median lo of the points within
# 2 % of the x range plus 0.5 d, or lo > FACTOR times their median hi (a traced run that belongs to another band, not a feature of
# the band; a steep rise or fall changes lo and hi together and is not caught). The rule is stated once and
# applied to all six targets; it removes only the two points of the right p90 band at x = 8.8772 and 9.9825 d
# (hi 2005 and 4117 where the neighbouring lo are about 47,000).
function outlier_filter(xs, los, his; span = 0.02, factor = 2.0)
    n = length(los)
    keep = trues(n)
    for i in 1:n
        nb = [j for j in 1:n if j != i && abs(xs[j] - xs[i]) <= span * (xs[end] - xs[1]) + 0.5]
        isempty(nb) && continue
        (his[i] < median(los[nb]) / factor || los[i] > median(his[nb]) * factor) && (keep[i] = false)
    end
    return keep
end

function main(outdir)
    traced = Dict(n => trace_panel(n) for n in ("left", "right"))
    for (id, panel, band, xlo, xhi, excl, step) in TARGETS
        cols, out = traced[panel]
        P = panels[panel]
        w = panel == "left" ? 11 : 5
        top = runmedian(out[band][:, 1], w); bot = runmedian(out[band][:, 2], w)
        xs = Float64[]; los = Float64[]; his = Float64[]
        for (j, c) in enumerate(cols)
            x = xdata(P.xcal, c)
            (x < xlo || x > xhi || any(e -> e[1] <= x <= e[2], excl)) && continue
            (isnan(top[j]) || isnan(bot[j])) && continue
            (c - first(cols)) % step == 0 || continue
            push!(xs, x); push!(his, ydata(P.ycal, top[j])); push!(los, ydata(P.ycal, bot[j]))
        end
        keep = outlier_filter(xs, los, his)
        for i in findall(!, keep)
            println("  dropped (post-review outlier rule) ", id, " x = ", round(xs[i]; digits = 4), " lo = ", los[i], " hi = ", his[i])
        end
        xs, los, his = xs[keep], los[keep], his[keep]
        write_band(joinpath(outdir, id * ".csv"), id, panel, band, xs, los, his, length(xs))
        println("wrote ", id, ".csv (", length(xs), " points)")
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(@__DIR__)
